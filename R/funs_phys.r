#' Compute temperature profile from TDR recordings
#' 
#' @param x the \code{x} variable.
#' @param y The \code{y} variable.
#' @param by To specify how to bin the x data.
#' @param na.rm Should the NAs be removed from the output ?
#' @export
#' @examples
#' data(exses)
#' temp_profiles <- tdrply(prof, c("depth", "temp"), by = 10, obj = exses)
#' plot(temp_profiles[["!_/#33"]], type = 'b')
prof <- function(x, y = NULL, by = 5, na.rm = FALSE) {
  if (is.recursive(x)) 
    nms <- names(x)
  else 
    nms <- c('x', 'y')
  xy <- as.data.frame(xy.coords(x, y))
  brk <- seq(min(xy$x, na.rm = TRUE), max(xy$x, na.rm = TRUE), by)
  out <- data.frame(rollapply(brk, mean, 2, nas = FALSE), 
                    tapply(xy$y, cut(xy$x, breaks = brk), mean, na.rm = TRUE))
  setNames(`if`(na.rm, na.omit(out), out), nms)
}

#' Use time and location to find if events occured during the day or the night
#' 
#' @inheritParams sun_position
#' @param elevlim Sun elevation thresholds to distinguish between day and night
#' @param type Should the output type be logical (\code{TRUE} for days, \code{FALSE} 
#' for night and \code{NA} otherwise) OR character (\code{"Day"}, \code{"Night"}, 
#' \code{"Transition"} and \code{NA} for missing locations).
#' @seealso \code{\link{sun_position}}
#' @export
#' @examples
#' data(exses)
#' do.call(is_day, exses$stat[ , c("time", "lat", "lon")])
is_day <- function (time, lat, lon, elevlim = c(-18, 18), type = c('character', 'logical')) {
  all(!is.na(time)) || stop('NA not allowed in "time" argument')
  
  loc <- data.frame(time, lat, lon)
  locNA <- apply(is.na(loc[ , -1]), 1, sum) != 0
  sunAngle <- rep(NA, nrow(loc))
  sunAngle[!locNA] <- do.call(sun_position, loc[!locNA, ])$el
  
  isday <- ifelse(sunAngle < min(elevlim), FALSE, NA)
  isday[sunAngle > max(elevlim)] <- TRUE
  if (match.arg(type) == 'character'){
    isday <- sapply(as.character(isday), 
                    function(x) switch(x, 'TRUE' = 'Day', 'FALSE' = 'Night', 'Transition'), 
                    USE.NAMES = FALSE)
    isday[locNA] <- NA
  }
  return(isday)
}

#' Compute sun azimuth and elevation given a location and a date/time
#' 
#' @param time a POSIXlt, POSIXct or data.frame. If a data frame order of 
#' columns must be from Year to Second (6 columns).
#' @param lat numeric vector of latitudes.
#' @param lon numeric vector of longitudes.
#' 
#' @details Based on: 
#' \url{http://stackoverflow.com/questions/8708048/position-of-the-sun-given-time-of-Day-latitude-and-longitude.}. 
#' Update version: Jan 6 '12 at 21:40 by "Josh O'Brien".
#' 
#' @seealso \code{\link{is_day}}
#' @keywords internal
#' @export
sun_position <- function(time, lat, lon) {
  !missing(lat) && !missing(lon) || stop('Please provide location information')
  
  if (inherits(time, "POSIXt")) {
    timelt <- as.POSIXlt(time)
    time <- data.frame(Day = timelt$mday, Month = timelt$mon + 1, Year = timelt$year + 1900, 
                       Hour = timelt$hour, Minute = timelt$min, Second = timelt$sec)
  } else {
    if (ncol(time) != 6) stop('"time" must have 6 columns.')
    time <- replaceMissing(time, na.0 = NA, 0)
    names(time) <- c("Year", "Month", "Day", "Hour", "Minute", "Second")
  }
  return(with(time, sunPos(Year, Month, Day, Hour, Minute, Second, lat, lon)))
}

#' sunPos
#' 
#' Compute sun azimuth and elevation given a location and a date/time. Algorithm taken from: 
#' \url{http://stackoverflow.com/questions/8708048/position-of-the-sun-given-time-of-Day-latitude-and-longitude.}. 
#' Update version: Jan 6 '12 at 21:40 by "Josh O'Brien".
#' 
#' @inheritParams sun_position
#' @seealso \code{\link{sun_position}}
#' @export
#' @keywords internal
sunPos <- function(Year, Month, Day, Hour, Minute, Second, Lat, Lon) {
  
  # From: http://stackoverflow.com/questions/8708048/position-of-the-sun-given-time-of-Day-latitude-and-longitude
  # Update: the Jan 6 '12 at 21:40 by "Josh O'Brien"
  
  twopi <- 2 * pi
  deg2rad <- pi / 180
  
  # Get Day of the Year, e.g. Feb 1 = 32, Mar 1 = 61 on leap Years
  Month.Days <- c(0,31,28,31,30,31,30,31,31,30,31,30)
  Day <- Day + cumsum(Month.Days)[Month]
  leapDays <- Year %% 4 == 0 & (Year %% 400 == 0 | Year %% 100 != 0) & 
    Day >= 60 & !(Month==2 & Day==60)
  Day[leapDays] <- Day[leapDays] + 1
  
  # Get Julian date - 2400000
  Hour <- Hour + Minute / 60 + Second / 3600 # Hour plus fraction
  delta <- Year - 1949
  leap <- trunc(delta / 4) # former leapYears
  jd <- 32916.5 + delta * 365 + leap + Day + Hour / 24
  
  # The input to the Atronomer's almanach is the difference between
  # the Julian date and JD 2451545.0 (noon, 1 January 2000)
  time <- jd - 51545.
  
  # Ecliptic coordinates
  
  # Mean longitude
  mnlong <- 280.460 + .9856474 * time
  mnlong <- mnlong %% 360
  mnlong[mnlong < 0] <- mnlong[mnlong < 0] + 360
  
  # Mean anomaly
  mnanom <- 357.528 + .9856003 * time
  mnanom <- mnanom %% 360
  mnanom[mnanom < 0] <- mnanom[mnanom < 0] + 360
  mnanom <- mnanom * deg2rad
  
  # Ecliptic longitude and obliquity of ecliptic
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.020 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong[eclong < 0] <- eclong[eclong < 0] + 360
  oblqec <- 23.439 - 0.0000004 * time
  eclong <- eclong * deg2rad
  oblqec <- oblqec * deg2rad
  
  # Celestial coordinates
  # Right ascension and declination
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num / den)
  ra[den < 0] <- ra[den < 0] + pi
  ra[den >= 0 & num < 0] <- ra[den >= 0 & num < 0] + twopi
  dec <- asin(sin(oblqec) * sin(eclong))
  
  # Local coordinates
  # Greenwich mean sidereal time
  gmst <- 6.697375 + .0657098242 * time + Hour
  gmst <- gmst %% 24
  gmst[gmst < 0] <- gmst[gmst < 0] + 24.
  
  # Local mean sidereal time
  lmst <- gmst + Lon / 15.
  lmst <- lmst %% 24.
  lmst[lmst < 0] <- lmst[lmst < 0] + 24.
  lmst <- lmst * 15. * deg2rad
  
  # Hour angle
  ha <- lmst - ra
  ha[ha < -pi] <- ha[ha < -pi] + twopi
  ha[ha > pi] <- ha[ha > pi] - twopi
  
  # Latitude to radians
  Lat <- Lat * deg2rad
  
  # Azimuth and elevation
  el <- asin(sin(dec) * sin(Lat) + cos(dec) * cos(Lat) * cos(ha))
  az <- asin(-cos(dec) * sin(ha) / cos(el))
  
  # For logic and names, see Spencer, J.W. 1989. Solar Energy. 42(4):353
  cosAzPos <- (0 <= sin(dec) - sin(el) * sin(Lat))
  sinAzNeg <- (sin(az) < 0)
  az[cosAzPos & sinAzNeg] <- az[cosAzPos & sinAzNeg] + twopi
  az[!cosAzPos] <- pi - az[!cosAzPos]
  
  el <- el / deg2rad
  az <- az / deg2rad
  Lat <- Lat / deg2rad
  
  return(list(elevation=el, azimuth=az))
}

#' Add averaged oceonagraphic fronts to an existing plot
#' 
#' @param name The name of the front(s) to add. To choose 
#' in \code{c('SACCF', 'PF', 'SAF', 'SSTF', 'NSTF')}
#' @param col The colors associated with \code{names}.
#' @param lwd The width of the line.
#' @param rescaleLon If not \code{NULL} then \code{\link{rescale}} is called
#' with these provided arguments.
#' @param rescaleLat If not \code{NULL} then \code{\link{rescale}} is called
#' with these provided arguments.
#' @param ... Other atguments to be passed to \code{\link{lines}}.
#' @details See the documentation of the \code{\link{sofronts}} dataset for more 
#' information about the data.
#' @export
#' @examples
#' plot(c(40, 120), c(-80,-30), type = 'n', xlab = '', ylab = '')
#' front()
front <- function(name = c('SACCF', 'PF', 'SAF', 'SSTF'), 
                  col = 1:4, 
                  lwd = 2, rescaleLon = NULL, rescaleLat = NULL, ...) {
  data(sofronts, envir = environment())
  names(col) <- name
  for (ft in name) {
    df <- sofronts[as.character(sofronts$Name) == ft & sofronts$Lon %bw% par('usr')[1:2], ]
    df$Name <- as.character(df$Name)
    if (!is.null(rescaleLon))
      df$Lon <- do.call(rescale, c(list(x = df$Lon, rescaleLon)))
    if (!is.null(rescaleLat))
      df$Lat <- do.call(rescale, c(list(x = df$Lat, rescaleLat)))
    lines(Lat ~ Lon, data = df, col = col[ft], lwd = lwd, ...)
  }
}
globalVariables("sofronts")

#' Add isobaths to an existing plot (Kerguelen area)
#' 
#' @param depths The isobath(s) which is(are) desired.
#' @param lty The line type.
#' @param col The line color.
#' @param plot Should the lines be drawn ?
#' @param ... Other arguments to be passed to \code{\link{lines}}.
#' @details Isobath at -1000 m corresponds approximatively to the Kerguelen shell 
#' boarder. See the documentation of the \code{\link{kerbathy}} dataset for more
#' information about the data.
#' @export
#' @examples
#' \dontrun{
#' # Require package 'sp' 
#' plot(c(68, 80), c(-54,-48), type = 'n', xlab = '', ylab = '')
#' isobath()
#' }
isobath <- function(depths  = -1000, lty = 2, col = 'gray', plot = TRUE, ...) {
  data(kerbathy, envir = environment())
  tab <- kerbathy[kerbathy$CONTOUR %in% depths, ]
  if (plot) lines(tab, lty = lty, col = col, ...)
  invisible(tab)
}
globalVariables("kerbathy")