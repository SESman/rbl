% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_phys.r
\name{is_day}
\alias{is_day}
\title{Use time and location to find if events occured during the day or the night}
\usage{
is_day(time, lat, lon, elevlim = c(-18, 18), type = c("character",
  "logical"))
}
\arguments{
\item{time}{a POSIXlt, POSIXct or data.frame. If a data frame order of 
columns must be from Year to Second (6 columns).}

\item{lat}{numeric vector of latitudes.}

\item{lon}{numeric vector of longitudes.}

\item{elevlim}{Sun elevation thresholds to distinguish between day and night}

\item{type}{Should the output type be logical (\code{TRUE} for days, \code{FALSE} 
for night and \code{NA} otherwise) OR character (\code{"Day"}, \code{"Night"}, 
\code{"Transition"} and \code{NA} for missing locations).}
}
\description{
Use time and location to find if events occured during the day or the night
}
\examples{
\dontrun{
do.call(is_day, exses$stat[ , c("time", "lat", "lon")])
}
}
\seealso{
\code{\link{sun_position}}
}

