% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_correction.r
\name{correct_duplicated_time}
\alias{correct_duplicated_time}
\title{Fix duplicates in TDR time stamps}
\usage{
correct_duplicated_time(obj, time_seq = 1, verbose = FALSE)
}
\arguments{
\item{obj}{The TDR dataset}

\item{time_seq}{A variable to use in order to subset \code{obj} and extract 
the time column. The time variable has to be in \code{POSIXct} format.}

\item{verbose}{A logical indicating if the function should be verbose.}
}
\description{
Fix duplicates in TDR time stamps
}
\details{
Table with time stamp and number of replicate of duplicates is 
returned in attributes.
}
\examples{
data(exses)
tmp <- exses$tdr[c(1,1,1,2:30), ]
tmp <- correct_duplicated_time(tmp)
attr(tmp, "correct_time")
}
\author{
Simon Wotherspoon, Yves Le Bras
}
\keyword{correction}
\keyword{internal}

