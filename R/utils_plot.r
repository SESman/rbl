#' @title Color helpers
#'
#' @param col A vector of colors to transform.
#' @param factor A value ranging between 0 and 1 decribing the strenght of
#' the effet.
#' @return \code{transparent} returns a transparent version of the input colors
#' with an alpha equal to the factor argument.
#' @keywords internal
#' @export
transp <- function (col, factor = 0.5) {
  c <- col2rgb(col)
  rgb(c['red', ], c['green', ], c['blue', ], factor * 255, maxColorValue = 255)
}

#' @rdname transp
#' @return \code{lighter} returns a brightened version of the input colors.
#' @export
lighter <- function (col, factor = 0.5) {
  c <- col2rgb(col)
  l <- function (c) 255 * factor + c * (1 - factor)
  rgb(l(c['red', ]), l(c['green', ]), l(c['blue', ]), maxColorValue = 255)
}

#' @rdname transp
#' @return \code{darker} returns a darkened version of the input colors.
#' @export
darker <- function (col, factor = 0.5) {
  c <- col2rgb(col)
  d <- function (c) c * (1 - factor)
  rgb(d(c['red', ]), d(c['green', ]), d(c['blue', ]), maxColorValue = 255)
}

#' Draw a correct axis for POSIXct objects
#' @param x The POSIXct data.
#' @param format The format for ticks label.
#' @inheritParams graphics::axis
#' @keywords internal
#' @export
axisPOSIXct <-function(x, side = 1, format = '%m/%d') {
  xax <- floorPOSIXct(x, 'days')
  xax <- unique(xax)
  xaxlab <- format(xax, format = format)
  xaxlab[!(seq_along(xax) %% 3) == 1] <- ''
  cond <- as.numeric(xax) > par('usr')[1]
  axis(side, at = xax[cond & xaxlab != ''], labels = xaxlab[cond & xaxlab != ''])
  
  opar <- par('xpd')
  par(xpd = FALSE) ; on.exit(par(opar))
  abline(v = xax[cond], lty = 3, lwd = .5)
}

#' plot for TDR data
#' 
#' Reverse Y axe when depth column in inputs
#' 
#' @param x TDR data
#' @param type as in \code{\link{plot}}
#' @param ... as in \code{\link{plot}}
#' @export
#' @keywords internal
plot.tdr <- function(x, type = "l", ...) {
  if (names(x)[1] == "depth") x <- x[ , 2:1]
  if (names(x)[2] == "depth") yl <- rev(range(x$depth)) else range(x[ , 2])
  plot(as.data.frame(x), type = type, ylim = list(...)[["ylim"]] %else% yl, ...)
}
