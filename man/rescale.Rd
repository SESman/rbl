% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_pkgtools.r
\name{rescale}
\alias{rescale}
\title{Scale a series between two values}
\usage{
rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE))
}
\arguments{
\item{x}{Numeric vector.}

\item{to}{Output range.}

\item{from}{Input range to be rescaled to \code{to}. Default is the range of \code{x}.}
}
\description{
\code{rescale} is a utility to resize the range of values while keeping 
the original spacing between values.
}
\examples{
x <- -10:10
rescale(x)
rescale(x, to = c(-1, 3))
rescale(x, from = c(5, max(x)), to = c(0, 10))
}
\keyword{internal}

