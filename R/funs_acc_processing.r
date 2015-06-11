#' Import data from Wildlife Computers ".tab" text files
#' 
#' @param x filename to be imported or a TDR dataset to be formated.
#' @param dt if TRUE function will return a data.table object, else, a data.frame.
#' @param ... Arguments to be passed to \code{\link{file}} such as \code{encoding}.
#' @details Wildlife Computers > Instrument helper > Save instrument readings > R format
#' @export
#' @keywords raw_processing
#' @import sqldf data.table
read.wcih <- function(x, dt = TRUE, ...) {
  stopifnot(require("data.table"))
  stopifnot(require("sqldf"))
  if (is.character(x)) {
    f <- file(x, ...)
    nms <- unlist(read.table(x, skip = 3, as.is = TRUE, nrows = 1))
    fmt <- list(skip = 4, header = FALSE, row.names = FALSE, sep = " ")
    x <- sqldf("select * from f", dbname = tempfile(), file.format = fmt)[ , -1]
  } else {
    nms <- names(x)
  }
  
  new_nms <- c(
    "Time" = "time", "Depth" = "depth", 
    "External Temperature" = "temp", "Light Level" = "light", 
    "int aX" = "ax", "int aY" = "ay", "int aZ" = "az", 
    "int mX" = "mx", "int mY" = "my", "int mZ" = "mz", 
    "Velocity" = "spd"
  )
  
  x <- setnames(as.data.table(x), new_nms[nms])
  x <- x[ , lapply(.SD, as.numeric)]
  x <- x[ , time := as.POSIXct(floor(time), origin = "1970-01-01", tz = 'UTC')]
  
  x <- if (dt) 
    data.table(x, key = "time")
  else 
    as.data.frame(x)
}

#' Identify Prey Catch attempts
#' 
#' @param x 3 axes acceleration table with time in the first column and acceleration 
#' axes in the following columns. Variables must be entitled "time" for time, 
#' "ax", "ay", and "az" for x, y and z accelerometer axes.
#' @param fs sampling frequency of the input data (Hz).
#' @param fc Cut-off frequency for the butterworth high pass filter (Hz)
#' @return returns a logical vector of prey catch attempts at 1 Hz frequency. 
#' Value is TRUE if the record belong to prey catch attempt FALSE otherwise.
#' @import data.table signal RcppRoll
#' @export
#' @keywords raw_processing
prey_catch_attempts <- function(x, fs = 16, fc = 2.64) {
  stopifnot(require("data.table"))
  stopifnot(require("signal"))
  stopifnot(require("RcppRoll"))
  # Generate Butterworth filter 
  # Critical frequencies of the filter: f_cutoff / (f_sampling/2)
  bf_pca  <- butter(3, W = fc / (0.5*fs), type = 'high')
  
  # Apply filter
  if (!is.data.table(x)) x <- data.table(x, key = "time")
  .f <- function(x) as.numeric(filter(bf_pca, x))
  x <- x[ , 2:4 := lapply(.SD, .f), .SDcols = 2:4]
  gc()
  
  # 1 s fixed window standard deviation + aggregate data to 1 Hz
  x <- x[ , lapply(.SD, sd, na.rm = TRUE), by = time]
  # In case of NAs set ACC to zero
  nas <- x[ , 2:4 := lapply(.SD, is.na), by = time, .SDcols = 2:4]
  nas_vector <- Reduce("|", nas[ , time := NULL])
  if (any(nas_vector)) {
    warning("NAs found and replaced by 0. Na proportion:", mean(nas_vector))
    x$ax[nas$ax] <- 0
    x$ay[nas$ay] <- 0
    x$az[nas$az] <- 0
  }
  gc()
  
  # 5 s moving window standard deviation
  .f <- function(x) c(0,0,roll_sd(x, 5),0,0)
  x <- x[ , 2:4 := lapply(.SD, .f), .SDcols = 2:4]
  gc()
  
  # kmean clutering: "high" = TRUE vs "low" = FALSE
  .f <- function(x) { 
    km_mod <- kmeans(x, 2)
    high_state <- which.max(km_mod$centers)
    as.logical(km_mod$cluster == high_state)
  }
  x <- x[ , 2:4 := lapply(.SD, .f), .SDcols = 2:4]
  
  # Aggregate and return to data.frame
  # records classified as PCA if the three axis are simultaneously in high state
  Reduce("&", x[ , time := NULL])
}

#' Compute swimming effort
#' 
#' @param fc Cut-off frequencies for the butterworth band pass filter (Hz)
#' @inheritParams prey_catch_attempts
#' @return returns a vector of swimming effort values at 1 Hz.
#' @details Only Y accelerometer axe is used to compute swimming effort.
#' @import data.table signal
#' @export
#' @keywords raw_processing
swimming_effort <- function(x, fs = 16, fc = c(0.4416, 1.0176)) {
  stopifnot(require("data.table"))
  stopifnot(require("signal"))
  # Generate a Butterworth filter 
  # Critical frequencies of the filter: f_filter / (f_sampling/2)
  bf_swm  <-  butter(3, W = fc / (0.5*fs), type = 'pass')
  
  # Apply filter
  if (!is.data.table(x)) x <- data.table(x, key = "time")
  x <- x[ , ay := abs(as.numeric(filter(bf_swm, ay)))]
  x <- x[ , c(2, 4) := NULL, with = FALSE] # remove unused "ax" & "az" columns
  
  # 1 s fixed window average + aggregate data to 1 Hz
  x <- x[ , lapply(.SD, mean, na.rm = TRUE), by = time]
  x <- x$ay
}

#' Static acceleration
#' 
#' The raw acceleration is first filtered using a low pass butterworth filter. 
#' The extracted signal is then scaled so that the norm of the the vector G is 1 at 
#' each second.
#' 
#' @param fc Cut-off frequency for the butterworth low pass filter (Hz)
#' @inheritParams prey_catch_attempts
#' @return returns a data.frame with time, and X, Y and Z static accelearyion at 1 Hz.
#' @details This filtered acceleration can be used to compute pitch and roll angles
#' @import data.table signal
#' @keywords raw_processing
#' @export
static_acc <- function(x, fs = 16, fc = 0.01) {
  stopifnot(require("data.table"))
  stopifnot(require("signal"))
  # Generate a Butterworth filter 
  # Critical frequencies of the filter: f_filter / (f_sampling/2)
  bf_grav  <-  butter(3, W = fc / (0.5*fs), type = 'low')
  
  # Apply filter
  if (!is.data.table(x)) x <- data.table(x, key = "time")
  x <- x[ , 2:4 := lapply(.SD, function(x) as.numeric(filter(bf_grav, x))), 
          .SDcols = 2:4]
  
  # 1 s fixed window average + aggregate data to 1 Hz
  x <- x[ , lapply(.SD, mean, na.rm = TRUE), by = time]
  x <- setnames(x, c('time', 'axG', 'ayG', 'azG'))
  
  # Scale axis
  Gnorm <- sqrt(x$axG^2 + x$ayG^2 + x$azG^2)
  x <- x[ , 2:4 := lapply(.SD, function(x) x / Gnorm), .SDcols = 2:4]

  as.data.frame(x)
}

#' Attitude angles from static accelation
#' 
#' @param object A data frame or TDR table including static acceleration variables 
#' entitled "axG", "ayG", and "azG" for X, Y, and Z axes of the accelerometer.
#' @export
#' @keywords raw_processing
pitch <- function(object) {
  -atan(object$axG/sqrt(object$ayG^2 + object$azG^2))
}

#' @rdname pitch
#' @export
#' @details For roll angle, the x axe is not necessary.
#' @return A vector of pitch/roll of the same length as \code{object}.
#' @keywords raw_processing
roll <- function(object) {
  atan2(object$ayG^2, object$azG)
}
