% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_acc_processing.r
\name{dynamic_acc}
\alias{dynamic_acc}
\title{Dynamic (Body) acceleration DBA}
\usage{
dynamic_acc(x, fs = 16, agg_1hz = TRUE, ...)
}
\arguments{
\item{x}{3 axes acceleration table with time in the first column and acceleration 
axes in the following columns. Variables must be entitled \code{"time"} for time, 
\code{"ax"}, \code{"ay"}, and \code{"az"} for x, y and z accelerometer axes.}

\item{fs}{sampling frequency of the input data (Hz).}

\item{agg_1hz}{Should the input be aggregated to 1 Hz ?}

\item{...}{Parameters to be passed to \code{\link{static_acc}} (e.g \code{fc}).}
}
\value{
returns a data.frame with time, and X, Y and Z static accelearyion at 1 Hz.
}
\description{
DBA is calculated by smoothing data for each axis to calculate the static 
acceleration (\code{\link{static_acc}}), and then subtracting it from 
the raw acceleration.
}
\details{
This filtered acceleration can be used to compute ODBA and VeDBA.
}
\keyword{raw_processing}

