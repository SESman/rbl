% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functionals.r
\name{curry}
\alias{curry}
\title{Curry}
\usage{
curry(FUN, ...)
}
\arguments{
\item{FUN}{a function}

\item{...}{arguments to be set}
}
\description{
Curry
}
\examples{
f <- curry(mean, na.rm = TRUE)
f(c(1, NA))
}
\keyword{functional}

