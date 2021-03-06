% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_delim.r
\name{bottom_delim_halsey}
\alias{bottom_delim_halsey}
\title{Delimitate bottom phase of dive using wiggles and steps}
\usage{
bottom_delim_halsey(time, depth = NULL, ledge = 0.5, vert_vel = 0.35)
}
\arguments{
\item{time}{time readings, sorted in chronological order.}

\item{depth}{depth readings, sorted in chronological order.}

\item{ledge}{depth threshold, specified as a percentage which is compared 
to the maximum depth. 0.75 has been used for king pengins. For 
elephant seals: about 95 percent of prey catch attemps occur at depth greater 
than 50 percent of the maximum depth (default value 0.50).}

\item{vert_vel}{A threshold used to define the steps. 0.35 m/s has been used for 
king pengins.}
}
\description{
Improves \code{\link{bottom_delim_std}} by defining a lower ledge (which implies 
longer bottoms) but reducing the bottom phase to a more relevant period 
(where specific behavior seems to indicate intensive foraging search). This 
is achieved by selecting the period between the first and last step/wiggle 
(see function \code{\link{wiggles}}) deeper than the ledge.
}
\examples{
data(exses)
ind(exses)

n <- 65
idx <- tdrply(bottom_delim_halsey, 1:2, "!_/", no = n)[[1]]
tdrply(plot, 1:2, "!_/", no = n, main = n)
tdrply(function(x, st, ed, ...) points(x[st:ed, ]), 1:2, no = n, la = idx)
}
\references{
Halsey, L.G., Bost, C.-A. & Handrich, Y. (2007) A thorough and 
quantified method for classifying seabird diving behaviour. 
Polar Biology, 30, 991-1004.
}
\seealso{
\code{\link{wiggles}}
}
\keyword{internal}

