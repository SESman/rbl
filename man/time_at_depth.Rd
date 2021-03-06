% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_behavior.r
\name{time_at_depth}
\alias{time_at_depth}
\alias{time_at_depth.bsm}
\alias{time_at_depth.default}
\alias{time_at_depth.ses}
\alias{time_at_depth.tdr}
\title{Compute the Time Allocation at Depth (TAD) index of a dive}
\usage{
time_at_depth(x, vs = Inf, fs = NULL, ...)

\method{time_at_depth}{default}(x, vs = Inf, fs = NULL, ...)

\method{time_at_depth}{tdr}(x, vs = Inf, ...)

\method{time_at_depth}{bsm}(x, vs = Inf, ...)

\method{time_at_depth}{ses}(x = ind(), vs = Inf, ...)
}
\arguments{
\item{x}{input data. Can be a numeric vector of depth records 
(then sampling frequency should be provided), a data frame with time and depth 
as columns 1 & 2 or some \code{bsm}/\code{tdr}/\code{ses} objects.}

\item{vs}{Maximum vertical speed achievable.}

\item{fs}{Sampling frequency (in Hz). Optional if a full time sequence is provided.}

\item{...}{\code{na.rm} or other arguments (mainly for S3 methods compatibility).}
}
\description{
Compute the Time Allocation at Depth (TAD) index of a dive
}
\details{
The index takes values from 0 for a dive where the maximum of time was spent
 at the minimum depth, to 0.5 for \code{"V"} shaped dives and 1 
for \code{"U"} shaped dives.
\code{"ses"} and \code{"tdr"} methods assume that time and depth information 
respectively occupy the first and second columns of the \code{tdr} tables.
}
\examples{
data(exses)
bsm_6pts <- tdrply(brokenstick, 1:2, obj = exses)
# These 4 lines return the same result
tad_highres <- tdrply(time_at_depth, 2, obj = exses, vs = 2, fs = 1)
tad_highres <- tdrply(time_at_depth, 1:2, obj = exses, vs = 2)
tad_highres <- sapply(bsm_6pts, time_at_depth, vs = 2) # "data" slot is used
exses$stat$tad <- tad_highres <- time_at_depth(exses, vs = 2)
# When the "data" slot is not available
tad_lowres <- sapply(eco.mem(bsm_6pts), time_at_depth, vs = 2) # "data" slot is not used

plot(tad_highres, tad_lowres) ; abline(0, 1, col = "red", lwd = 3)
plot(tad ~ time, exses$stat)
}
\references{
Fedak, M. A., Lovell, P. and Grant, S. M. (2001). Two Approaches 
to Compressing and Interpreting Time-Depth Information as as Collected by 
Time-Depth Recorders and Satellite-Linked Data Recorders. 
Marine Mammal Science 17, 94--110.
}

