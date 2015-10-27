#' @title Color helpers
#' @description Set alpha transparency parameter of input colors.
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
#' When depth in input set in as Y and reverse axe. Use third column as a color 
#' variable.
#' 
#' @param x TDR data
#' @param type as in \code{\link{plot}}
#' @param ... as in \code{\link{plot}}
#' @export
#' @keywords internal
#' @examples 
#' data(exses)
#' tdrply(plot, c("time", "depth", "is_pca"), no = 65, obj = exses)
#' 
#' exses$tdr$no_btt <- tdrexpand(exses$stat$no_dive, ty = "_", obj = exses)
#' exses$tdr$is_btt <- !is.na(exses$tdr$no_btt)
#' tdrply(plot, c("time", "depth", "is_btt"), no = 65, obj = exses)
#' tdrply(plot, c("depth", "temp", "is_btt"), no = 65, obj = exses)
plot.tdr <- function(x, type = "l", ...) {
  # If depth is involved set it to Y and inverse axe
  if (names(x)[1] == "depth") x <- cbind(x[ , 2:1], x[ , -(1:2)])
  if (names(x)[2] == "depth") yl <- rev(range(x$depth, na.rm = TRUE))
  if (exists("yl") && all(is.na(yl))) stop("depth variable contains only NA.")
  dots <- list(...)
  if (exists("yl") && !grepl("ylim", names(dots) %else% "")) dots[["ylim"]] <- yl
  
  # If 3 variables in input and no color specification use tird column
  if (ncol(x) == 3 && !grepl("col", names(dots) %else% "")) {
    f <- x[ , 3]
    if (is.logical(f)) {col <- f + 1} 
    else if (is.factor(f) || is.character(f)) {col <- as.integer(as.factor(f))} 
    else if (is.numeric(f)) {col <- grey(rescale(f))}
    dots[["type"]] <- "p"
    dots[["col"]] <- col
    do.call(plot, c(list(x = as.data.frame(x[ , 1:2])), dots))
    } else {
    do.call(plot, c(list(x = as.data.frame(x), type = type), dots))
  }
  invisible(NULL)
}
