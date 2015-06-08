#' Count the number of wiggle in a dive
#' 
#' \code{wiggles} ...
#' 
#' @param x The x data. A subset of a TDR dataset time sequence. Can also 
#' be a list of \code{x} and \code{y} (processed by \code{\link{xy.coords}}). 
#' @param y The y data. A subset of a TDR dataset depth sequence, processed by 
#' \code{\link{xy.coords}}).
#' @param thres.depth minimum depth difference.
#' @param thres.dur minimum time difference.
#' @param plt Should graphics about processing be plotted ?
#' @export
#' @keywords behavior
#' @examples
#' data(exses)
#' sunflowerplot(tdrply(wiggles, ty = '_', obj = exses), exses$stat$pca)
wiggles <- function(x, y = NULL, thres.depth = 2, thres.dur = c(5, 300), plt = FALSE) {
  xy <- as.data.frame(xy.coords(x, y)[1:2])
  bsm <- optBrokenstick(xy, threshold = thres.depth, cost = max_dist_cost)
  
  slp <- per(sign(coef(bsm)$slope))
  if (slp$value[nrow(slp)] == -1) slp <- slp[-nrow(slp), ]
  if (nrow(slp) <= 1) return(0)
  if (slp$value[1] == 1) slp <- slp[-1, ]
  if (nrow(slp) <= 1) return(0)
  slp$no <- cumsum(seq_along(slp$st_idx) %% 2)
  
  wgl_idx <- by(slp, slp$no, function(x) seq(min(x$st_idx), max(x$ed_idx)), simplify = FALSE)
  pts_stck <- which.stick(bsm, bsm$pts, type = 'i')
  wgl_st <- sapply(wgl_idx, function(x) min(bsm$pts[pts_stck %in% x]))
  wgl_ed <- sapply(wgl_idx, function(x) max(bsm$pts[pts_stck %in% x]))
  dt <- mapply(function(st, ed) diff(xy$x[c(st, ed)]), wgl_st, wgl_ed)
  dp <- mapply(function(st, ed) diff(range(xy$y[st:ed])), wgl_st, wgl_ed)
  cond_time <- dt %bw% (thres.dur %else% c(0, Inf))
  if (length(thres.depth) != 2)
    cond_depth <- dp >= (thres.depth %else% 0)
  else 
    cond_depth <- dp %bw% thres.depth
  slp <- slp[rep(cond_time & cond_depth, each = 2), ]
  
  if (plt) {
    plot(xy, type = 'l', ylim = rev(range(xy$y)))
    plot(bsm, add = TRUE)
    segments(x0 = xy$x[wgl_st], x1 = xy$x[wgl_ed], y0 = max(xy$y), lwd = 2)
    segments(x0 = xy$x[wgl_st], y0 = max(xy$y), y1 = xy$y[wgl_st], lty = 2, col = 'gray')
    segments(x0 = xy$x[wgl_ed], y0 = max(xy$y), y1 = xy$y[wgl_ed], lty = 2, col = 'gray')
  }
  
  return(nUN(slp$no))
}

#' Compute straightness index
#' 
#' \code{straightness} compute a straightness index between 0 and 1 by making the 
#' ration \code{L / l} where \code{L} is the cumulated distance between 
#' brokenstick points and and where \code{l} the cumulated distance between 
#' each \code{y} data points.
#' 
#' @param x The x data (time).
#' @param y The y data (depth).
#' @param npts The number of points to use. 2 is the minimum (from the start 
#' to the end of data). Each new point adds a step (using the \code{\link{brokenstick}} 
#' algorithm) which is taken into account when computing \code{L}.
#' @return \code{straightness} returns a number between 0 (maximum sinuosity) 
#' to 1 (maximum straightness). \code{sinuosity} is equivalent to 
#' \code{1 / straightness}.
#' @export
#' @keywords behavior
#' @examples
#' data(exses)
#' ind(exses)
#' sunflowerplot(tdrply(straightness, cl = 1:2, ty = '_'), exses$stat$pca)
straightness <- function(x, y = NULL, npts = 3) {
  bsm <- brokenstick(x, y, npts)
  L <- sum(abs(diff(bsm$data[bsm$pts, 2])))
  l <- sum(abs(diff(bsm$data[ , 2])))
  L / l
}

#' @rdname straightness
#' @export
sinuosity <- function(x, y = NULL, pts = 3) {
  1 / straightness(x, y, pts)
}

#' Shannon entropy index on time at depth proportions
#' 
#' @param x a character vector naming the depth layers.
#' @param base base argument passed to \code{\link{log}}.
#' @param scale if TRUE the output is divided by the \code{log(N)} where N is the 
#' number of layers so that outpout lies between 0 and 1.\code{log(N)} is the 
#' minimum entropy for sample with N layers where each observation of sample is 
#' unique i.e all probability are equal.
#' @export
#' @keywords behavior
#' @examples 
#' data(exses)
#' ind(exses)
#' brks <- do.call(seq, as.list(c(range(exses$tdr$depth), by = 2)))
#' exses$tdr$depth_cat <- as.character(cut(exses$tdr$depth, brks))
#' plot(tdrply(entropy, "depth_cat", ty = "!_/"), exses$stat$pca)

#' H <- tdrply(entropy, "depth_cat", ty = "_")
#' scaled_H <- tdrply(entropy, "depth_cat", ty = "_", scale = TRUE)
#' plot(H, exses$stat$pca)
#' plot(H, scaled_H)
entropy <- function(x, base = 2, scale = FALSE) {
  pi <- tapply(x, x, length) / length(x)
  H <- -sum(pi * log(pi, base))
  "if"(scale, H / log(nUN(x), base), H)
}
