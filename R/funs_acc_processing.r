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
    "Velocity" = "spd", "Internal Temperature" = "itemp"
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
#' Method based on Viviant et al. (2010) (see references) which was originally 
#' implemented on a jaw accelerometer.
#' 
#' @param x 3 axes acceleration table with time in the first column and acceleration 
#' axes in the following columns. Variables must be entitled \code{"time"} for time, 
#' \code{"ax"}, \code{"ay"}, and \code{"az"} for x, y and z accelerometer axes. 
#' @param fs sampling frequency of the input data (Hz).
#' @param fc Cut-off frequency for the butterworth high pass filter (Hz). Frequency 
#' above which signals across x, y and z accelerometer axes are retained. 
#' @param w window width in seconds to use when applying rolling standard deviation.
#' @param grp A vector (with same number of observation than \code{x}) 
#' which provides a group identifier for observations to be treated by rolling 
#' standard deviation separately. For example a dive number 
#' or a \code{\link{brokenstick}}  segment number. 
#' Leave as \code{NULL} to process all data as a single block.
#' @param n_days time period (in days) over which kmeans clustering of negative 
#' and positive prey catch attempt behaviours should be performed. 
#' Leave as \code{NULL} to have no time grouping in kmeans generation.
#' @return returns a logical vector of prey catch attempts at 1 Hz frequency. 
#' Value is TRUE if the record belong to prey catch attempt FALSE otherwise.
#' @import data.table signal RcppRoll
#' @author Yves Le Bras, Samantha Cox
#' @references Viviant, M., Trites, A. W., Rosen, D. A. S., Monestiez, P. and Guinet, C. 
#' #(2010). Prey capture attempts can be detected in Steller sea lions 
#' and other marine predators using accelerometers. Polar Biol 33, 713–719.
#' @export
#' @keywords raw_processing
prey_catch_attempts <- function(x, fs = 16, fc = 2.64, w = 1.5, grp = NULL, n_days = NULL) {
  stopifnot(require("data.table"))
  stopifnot(require("signal"))
  stopifnot(require("RcppRoll"))
  if (!is.data.table(x)) x <- data.table(x, key = "time")
  
  # Cut time into groups of "n_days" each
  time.bak <- copy(x$time)
  if (!is.null(n_days)) {
    days <- as.numeric(time.bak - time.bak[1]) / (24*3600)
    n_day <- max(days)
    n_bins <- n_day %/% n_days + ((n_day %% n_days) != 0)
    brks <- seq(0, n_bins) * n_days
    time_grp <- .bincode(days, breaks = brks, right = FALSE, include.lowest = TRUE)
    rm(days, n_days, n_bins, brks)
  }
  
  # Butterworth filter: Critical frequencies of the filter: f_cutoff / (f_sampling/2)
  bf_pca <- signal::butter(3, W = fc / (0.5 * fs), type = 'high')
  # Apply filter
  .f <- function(x) as.numeric(signal::filtfilt(bf_pca, x))
  x <- x[ , 2:4 := lapply(.SD, .f), .SDcols = 2:4]
  gc()
  
  # Handle missing values
  nas <- lapply(x[, 2:4, with = FALSE], is.na)
  nas_vector <- Reduce("|", nas)
  if (any(nas_vector)) {
    warning("NAs found and replaced by 0. NA proportion:", 
            round(mean(nas_vector), digits = 4))
    x$ax[nas$ax] <- x$ay[nas$ay] <- x$az[nas$az] <- 0
  }
  gc()
 
  # Apply rolling standard deviation
  .f <- function(x) roll_sd(x, fs * w, fill = 0)
  if (is.null(grp)) {
    x <- x[, `:=`(2:4, lapply(.SD, .f)), .SDcols = 2:4]
  } else {
    x <- x[, `:=`(2:4, lapply(.SD, .f)), by = grp, .SDcols = 2:4]
  }
  gc()
  
  # Apply 2-means clustering to the standard deviation signal 
  .f <- function(x) { 
    km_mod <- kmeans(x, 2)
    high_state <- which.max(km_mod$centers)
    as.logical(km_mod$cluster == high_state)
  }
  if (is.null(n_days)) {
    x <- x[, `:=`(2:4, lapply(.SD, .f)), .SDcols = 2:4]
  } else {
    x <- x[, `:=`(2:4, lapply(.SD, .f)), by = time_grp, .SDcols = 2:4]
  }
  
  # Aggregate to 1 Hz. Rule: classify a second as a PCA if at least one TRUE
  x <- x[, lapply(.SD, any), by = time, .SDcols = 2:4]
  
  # Classify obs as PCA if the three axis are simultaneously in high state
  Reduce("&", x[, `:=`(time, NULL)])
}

#' Compute swimming effort
#' 
#' @param fc Cut-off frequencies for the butterworth band pass filter (Hz)
#' @inheritParams prey_catch_attempts
#' @param rms Should the root mean square be used (instead of mean of absolute values) 
#' when averaging the acceleration to 1 Hz ?
#' @return returns a vector of swimming effort values at 1 Hz.
#' @details Only Y accelerometer axe is used to compute swimming effort.
#' @import data.table signal
#' @export
#' @keywords raw_processing
swimming_effort <- function(x, fs = 16, fc = c(0.4416, 1.0176), rms = FALSE) {
  stopifnot(require("data.table"))
  stopifnot(require("signal"))
  # Generate a Butterworth filter 
  # Critical frequencies of the filter: f_filter / (f_sampling/2)
  bf_swm <- signal::butter(3, W = fc / (0.5*fs), type = 'pass')
  
  # Apply filter
  if (!is.data.table(x)) x <- data.table(x, key = "time")
  x <- x[ , ay := abs(as.numeric(signal::filtfilt(bf_swm, ay)))]
  x <- x[ , c(2, 4) := NULL, with = FALSE] # remove unused "ax" & "az" columns
  
  # 1 s fixed window average + aggregate data to 1 Hz
  if (!rms) {
    x <- x[ , lapply(.SD, function(x) mean(abs(x), na.rm = TRUE)), by = time]
  } else {
    x <- x[ , lapply(.SD, function(x) sqrt(mean(x^2, na.rm = TRUE))), by = time]
  }
  x <- x$ay
}
globalVariables("ay")

#' Static acceleration
#' 
#' The raw acceleration is first filtered using a low pass butterworth filter. 
#' Then , the extracted signal can be scaled so that the norm of the the vector 
#' G is 1 at each second.
#' 
#' @param fc Cut-off frequency for the butterworth low pass filter (Hz)
#' @param Gscale Should the values be scaled by the norm of the static 
#' acceleration vector ?
#' @param agg_1hz Should the input be aggregated to 1 Hz ?
#' @inheritParams prey_catch_attempts
#' @return returns a data.frame with time, and X, Y and Z static accelearyion at 1 Hz.
#' @details This filtered acceleration can be used to compute pitch and roll angles
#' @import data.table signal
#' @keywords raw_processing
#' @export
static_acc <- function(x, fs = 16, fc = 0.20, Gscale = TRUE, agg_1hz = TRUE) {
  stopifnot(require("data.table"))
  stopifnot(require("signal"))
  # Generate a Butterworth filter 
  # Critical frequencies of the filter: f_filter / (f_sampling/2)
  bf_grav <- signal::butter(3, W = fc / (0.5*fs), type = 'low')
  
  # Apply filter
  if (!is.data.table(x)) x <- data.table(x, key = "time")
  x <- x[ , 2:4 := lapply(.SD, function(x) as.numeric(signal::filtfilt(bf_grav, x))), 
          .SDcols = 2:4]
  
  # 1 s fixed window average + aggregate data to 1 Hz
  if (agg_1hz) {
    x <- x[ , lapply(.SD, mean, na.rm = TRUE), by = time]
  }
  x <- setnames(x, c('time', 'axG', 'ayG', 'azG'))
  
  # Scale axis
  if (Gscale) {
    Gnorm <- sqrt(x$axG^2 + x$ayG^2 + x$azG^2)
    x <- x[ , 2:4 := lapply(.SD, function(x) x / Gnorm), .SDcols = 2:4]
  }
  
  as.data.frame(x)
}

#' Dynamic (Body) acceleration DBA
#' 
#' DBA is calculated by smoothing data for each axis to calculate the static 
#' acceleration (\code{\link{static_acc}}), and then subtracting it from 
#' the raw acceleration.
#' 
#' @param ... Parameters to be passed to \code{\link{static_acc}} (e.g \code{fc}).
#' @inheritParams static_acc
#' @return returns a data.frame with time, and X, Y and Z static accelearyion at 1 Hz.
#' @details This filtered acceleration can be used to compute ODBA and VeDBA.
#' @import data.table signal
#' @export
#' @keywords raw_processing
dynamic_acc <- function(x, fs = 16, agg_1hz = TRUE, ...) {
  static <- static_acc(copy(x), fs = fs, Gscale = FALSE, agg_1hz = FALSE, ...)
  x <- x[ , `:=`(2:4, Map("-", x[ , 2:4, with = FALSE], static[ , 2:4])), with = FALSE]
  rm(list = "static") ; gc()
  if (agg_1hz) {
    x <- x[ , lapply(.SD, mean, na.rm = TRUE), by = time]
  }
  as.data.frame(setnames(x, c("time", "axD", "ayD", "azD")))
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

#' Overall Dynamic Body Acceleration (ODBA)
#' 
#' @param object A data frame or TDR table including dynamic acceleration variables 
#' entitled "axD", "ayD", and "azD" for X, Y, and Z axes of the accelerometer.
#' @export
#' @return A vector of ODBA of the same length as \code{object}.
#' @keywords raw_processing
overall_DBA <- function(object) {
  object <- as.data.table(object)
  object <- object[ , tot := abs(axD) + abs(ayD) + abs(azD)]
  object$tot
}

#' Vectorial Dynamic Body Acceleration (VeDBA)
#' 
#' @param object A data frame or TDR table including dynamic acceleration variables 
#' entitled "axD", "ayD", and "azD" for X, Y, and Z axes of the accelerometer.
#' @export
#' @return A vector of VeDBA of the same length as \code{object}.
#' @keywords raw_processing
vectorial_DBA <- function(object) {
  object <- as.data.table(object)
  object <- object[ , tot := sqrt(axD^2 + ayD^2 + azD^2)]
  object$tot
}
