#' Fitter function for brokenstick models
#' 
#' Basic computing engines of the brokenstick algorithm called by 
#' \code{\link{brokenstick}}.
#' 
#' @param xy Data, as returned by \code{\link{xy.coords}}. The \code{x} values 
#' have to be sorted.
#' @param pts The points to use.
#' @param eco.mem A value between 0 and 6 to control the memory size of the 
#' output. \code{0} = no memory savings.
#' @seealso \code{\link{brokenstick}} and \code{\link{optBrokenstick}} which should 
#' be used to fit brokenstick models.
#' @keywords internal
#' @export
bsmfit <- function(xy, pts, eco.mem = 0L) 
{ 
	if (eco.mem > 6) stop('"eco.mem" > 6. Not output !')
	n <- diff(pts)
	n[1L] <- n[1L] + 1 # First segment needs 1 more prediction for the first point
	dy <- diff(xy[pts, 2L]) ; dx <- diff(xy[pts, 1L])
	a <- dy / dx
	b <- xy[pts[-length(pts)], 2L] - (a * xy[pts[-length(pts)], 1L])
	yfit <- rep(a, n) * xy[ , 1L] + rep(b, n)
	
	out <- list(pts.x = xy[pts, 1L], slope = a, intercept = b, 
				pts = pts, residuals = xy[, 2] - yfit, 
				fitted.values = yfit, data = xy)
	`if`(eco.mem != 0L, out[-unique(seq(7L - eco.mem + 1L, 7L))], out)
}

#' Fitting brokenstick models
#' 
#' \code{brokenstick} is used to fit brokenstick models on two-dimensional data 
#' such as Time-Depth, Depth-Temperature or Depth-Light profiles. Brokenstick models 
#' are useful for compressing or extracting the shape of high resolution profiles.
#' 
#' @param x The x data. Note that the \code{x} values have to be sorted. Can also 
#' be a list of \code{x} and \code{y} (processed by \code{\link{xy.coords}}). 
#' Alternatively, it can also be an object of class "\code{\link{formula}}".
#' @param y The y data. (processed by \code{\link{xy.coords}}). Alternatively, 
#' if \code{x} is a \code{'formula'}, \code{y} can also be an optional data frame, 
#' list or environment (or object coercible by \code{\link{as.data.frame}} to a 
#' data frame) containing the variables in the model. If not found in data, the 
#' variables are taken from environment (\code{formula}), typically the environment 
#' from which \code{brokenstick} is called.
#' @param npts The number of points for the brokenstick model to fit. See 
#' \code{\link{optBrokenstick}} for a version of \code{brokenstick} which figures 
#' out this number of points automatically.
#' @param start Some starting points to start the algorithm with.
#' @param na.action A function which indicates what should happen when the data 
#' contain \code{NAs}. The default is set by the \code{na.action} setting of 
#' \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The 
#' "factory-fresh" default is \code{\link{na.omit}}. Another possible value 
#' is \code{NULL}, no action. Value \code{\link{na.exclude}} can be useful.
#' @param ... Further arguments to be passed to \code{\link{bsmfit}} 
#' such as \code{eco.mem}.
#' @return A \code{bsm} object with (depending on \code{eco.mem}) contains:
#' \itemize{
#'   \item pts.x (\code{eco.mem < 7}) The x values of the brokenstick points.
#'   \item slope (\code{eco.mem < 6}) The slopes of each "stick" of the model.
#'   \item intercept (\code{eco.mem < 5}) The intercept of each "stick" of the model.
#'   \item pts (\code{eco.mem < 4}) The index of the brokenstick points.
#'   \item na.action Information from the action which was applied to object if \code{NAs} 
#'   were handled specially. (\code{eco.mem} inefficient)
#'   \item residuals (\code{eco.mem < 3}) The model's residuals.
#'   \item fitted.values (\code{eco.mem < 2}) The fitted values.
#'   \item data (\code{eco.mem < 1}) The input data used for fitting.
#'   \item pts.no The iteration number of points. (\code{eco.mem} inefficient)
#' }
#' @details See \code{\link{bsmfit}}, the function called by \code{brokenstick} 
#' on each iteration to fit the model with specified points.
#' @seealso \code{\link{optBrokenstick}} and \code{\link{predict.bsm}}, \code{\link{residuals.bsm}}, 
#' \code{\link{update.bsm}}, \code{\link{summary.bsm}}, 
#' \code{\link{coef.bsm}}, \code{\link{plot.bsm}}, \code{\link{as.data.frame.bsm}}
#' for other functions with a S3 method for \code{bsm} objects. 
#' @export
#' @keywords brokenstick
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' 
#' # Syntax
#' bsm <- brokenstick(dv$time, dv$depth)
#' bsm <- brokenstick(dv) # if two columns
#' bsm <- brokenstick(depth ~ time, dv) 
#' bsm <- with(dv, brokenstick(depth ~ time))
brokenstick <- function(x, y, npts = 6, start  = NULL, na.action, ...) {
	UseMethod("brokenstick")
}

#' @inheritParams brokenstick
#' @export
brokenstick.default <- function(x, y = NULL, npts = 6, start = NULL, 
                                na.action, ...) {
  if (is.recursive(x)) 
    nms <- names(x)
  else 
    nms <- c(deparse(substitute(x)), deparse(substitute(y)))
  if (any(grepl('\\$', nms)))
    nms <- gsub('(.*\\$)(.*$)', '\\2', nms)
  else if (any(sapply(nms, nchar) > 10))
    nms <- c('x', 'y')
	xy <- setNames(as.data.frame(xy.coords(x, y)[1:2]), nms)
	if (missing(na.action)) na.action <- options("na.action")[[1]]
	xy <- do.call(na.action, list(xy))
	pts <- `if`(is.null(start) || length(start) < 2, c(1, length(xy[ , 1L])), start)
  np <- length(pts)
	pts.no <- rep(1L, np)
	while (np < npts) {
	  brkstk <- bsmfit(xy, pts)
	  absRes <- abs(brkstk$residuals)
	  pts <- c(pts, which.max(absRes))
	  pts.no <- c(pts.no, max(pts.no) + 1L)
    dup <- duplicated(pts)
	  pts <- pts[!dup] ; pts.no <- pts.no[!dup]
	  ord <- order(pts) ; pts <- pts[ord] ; pts.no <- pts.no[ord]
	  if (any(dup)) {
	    warning('Duplicated points found. "npts" may be too high.')
	    break
	  }
	  np <- length(pts)
	}
	out <- bsmfit(xy, pts, ...)
  out$pts.no <- pts.no
	out$na.action <- attr(xy, "na.action")
	class(out) <- c('bsm', 'list')
	out
}

#' @inheritParams brokenstick
#' @keywords internal
#' @export
brokenstick.formula <- function(x, y = NULL, npts = 6, start = NULL, 
                                na.action, ...) {
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("x", "y"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf[[1L]] <- quote(stats::model.frame)
  names(mf) <- c('', 'formula', 'data')[seq_along(mf)]
	mf <- eval(mf, parent.frame())
	if (missing(na.action)) na.action <- options("na.action")[[1]]
	brokenstick.default(mf[2:1], NULL, npts = npts, start = start, na.action, ...)
}

#' Update and Re-fit a brokenstick model
#' 
#' @param object Object of class inheriting from "\code{bsm}".
#' @inheritParams brokenstick
#' @param allow.dup Should the duplicated points be added anyway.
#' @export
#' @keywords internal brokenstick
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' bsm <- brokenstick(dv) 
#' length(bsm$pts)
#' length(update(bsm, 5)$pts)
#' length(update(bsm, 7)$pts)
update.bsm <- function(object, npts, allow.dup = FALSE, ...) 
{
  np <- length(object$pts)
  if (np == npts) return(object)
  if (!'data' %in% names(object)) stop('Please provide data ("eco.mem = 0").')
  if (np < npts) {
    pts <- object$pts ; pts.no <- object$pts.no
    while (np < npts) {
      brkstk <- bsmfit(object$data, pts)
      absRes <- abs(brkstk$residuals)
      pts <- c(pts, which.max(absRes))
      pts.no <- c(pts.no, max(pts.no) + 1L)
      pts <- pts[ord <- order(pts)] ; pts.no <- pts.no[ord]
      dup <- duplicated(pts)
      if (any(dup)) {
        warning('Duplicated points found. "npts" may be too high.')
        if (!allow.dup) {
          pts <- pts[!dup] ; pts.no <- pts.no[!dup]
          break
        }
      }
      np <- length(pts)
    }
  } else {
    cond <- (object$pts.no <= npts)
    k <- 1
    while (sum(cond) > npts) {
      cond <- (object$pts.no <= npts - k)
      k <- k + 1
    }
    pts <- object$pts[cond] ; pts.no <- object$pts.no[cond]
  }
  x <- bsmfit(object$data, pts, ...)
  x$pts.no <- pts.no
  x$na.action <- `if`('na.action' %in% names(object), object$na.action, NULL)
  class(x) <- c('bsm', 'list')
  x
}

#' Brokenstick model predictions
#' 
#' @param object Object of class inheriting from "\code{bsm}".
#' @param newdata An optional data frame in which to look for variables with 
#' which to predict. If omitted, the fitted values are used.
#' @param ... other arguments.
#' @export
#' @details Require \code{eco.mem <= 4} (see \code{\link{bsmfit}}). When 
#' using \code{newdata} argument, the function returns \code{NAs} for values that 
#' do not belong to a stick.
#' @seealso \code{\link{brokenstick}}, \code{\link{optBrokenstick}}, 
#' \code{\link{which.stick}}
#' @keywords internal brokenstick
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' bsm <- brokenstick(dv) 
#' plot(depth ~ time, dv, ylim = rev(range(dv$depth)), type = 'l')
#' plot(bsm, add = TRUE, enumerate = TRUE)
#' 
#' ypred <- predict(bsm, newdata = xpred <- sample(dv$time, 10))
#' points(xpred, ypred, col = 3, pch = 19)
predict.bsm <- function(object, newdata, ...) {
	if (missing(newdata)) {
		if (is.null(fitted(object))) {
			if (!'data' %in% names(object)) {
				stop('Impossible without x values. ', 
					 'Consider add theses or set "eco.mem"  less than 2.')
			} else {
				# This case should not occur as 'data' is missing befor 'fitted' when
				# using a 'eco.mem' value  greater than zero.
				n <- diff(object$pts)
				n[1] <- n[1] + 1 
				return(rep(object$slope, n) * object$data$x + rep(object$intercept, n))
			}
		} else {
			return(fitted(object))
		}
	} else {
		stks <- which.stick(object, newdata)
		if (is(newdata, 'POSIXct')){newdata <- as.numeric(newdata)}
		return(object$slope[stks] * newdata + object$intercept[stks])
	}
	warning('Unexpected case occured in "predict.brokenstick".')
	NA
}

#' To which brokenstick segment a point belongs to ?
#' 
#' Given a brokenstick model and a set of points, the function determines on which 
#' stick the points are located.
#' 
#' @param object Object of class inheriting from "\code{bsm}".
#' @param pts The set of the points to match aginst sticks.
#' @param type The type of values provided in \code{pts}. To choose in 
#' \code{c('x', 'i')} where \code{'x'} stands for x values and \code{'i'} stands 
#' for the index of values.
#' @export
#' @seealso \code{\link{predict.bsm}}
#' @keywords internal brokenstick
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' bsm <- brokenstick(dv) 
#' (pts <- sample(1:nrow(dv), 5))
#' which.stick(bsm, pts, type = 'i')
#' which.stick(bsm, dv[pts, 1], type = 'x')
#' 
#' \dontrun{
#' # For the actual points of the model the result does not matter so much
#' # since both of previous and next segment is valid for prediction.
#' which.stick(bsm, bsm$pts, 'i')
#' }
which.stick <- function(object, pts, type = c('x', 'i')) {
	bsm.pts <- switch(match.arg(type), x = object$pts.x, i = object$pts)
	.f <- function(eql, grt, lst) {
		if (any(eql)) {ifelse(which(eql) == length(eql), which(eql) - 1, which(eql))}
		else if (is.na(grt) || is.na(lst)) {NA}
		else {grt}
	}
	eql <- lapply(pts, function(x) x == bsm.pts)
	grt <- lapply(pts, function(x) max(which(x >= bsm.pts) %else% NA))
	lst <- lapply(pts, function(x) min(which(x <= bsm.pts) %else% NA))
	mapply(.f, eql, grt, lst)
}

#' Extract brokenstick models coefficients
#' 
#' @param object \code{bsm} object, typically result from \code{\link{brokenstick}} 
#' or \code{\link{optBrokenstick}}.
#' @param ... other arguments.
#' @return A data frame with slope and intercepts of each stick.
#' @seealso \code{\link{brokenstick}}, \code{\link{optBrokenstick}} for model fitting.
#' @export
#' @keywords internal brokenstick
#' @details Require \code{eco.mem <= 4} (see \code{\link{bsmfit}}).
#' @examples
#' data(ses)
#' dv <- tdrply(identity, 1:2, no = 400, obj= exses)[[1]]
#' coef(brokenstick(dv))
coef.bsm <- function(object, ...) {
	out <- data.frame(intercept = object$intercept, slope = object$slope)
	row.names(out) <- paste0('seg', seq(nrow(out)))
	out
}

#' Plot brokenstick models
#' 
#' @param x \code{bsm} object to plot, typically result from \code{\link{brokenstick}} 
#' or \code{\link{optBrokenstick}}.
#' @param add If true add the plot to the already existing plot.
#' @param type Character indicating the type of plotting; actually any of the 
#' types as in \code{\link{plot.default}}.
#' @param lwd The line width, a positive number, defaulting to 2.
#' @param ylim the y limits (y1, y2) of the plot. 
#' Here y1 > y2 and leads to a "reversed axis".
#' @param col A specification for the default plotting color.
#' @param enumerate A switch to indicate if the iteration number of points should 
#' be added to the plot.
#' @param ... Further graphical parameters (see \code{\link{par}}) may also be 
#' supplied as arguments.
#' @details Require \code{eco.mem <= 4} (see \code{\link{bsmfit}}).
#' @seealso \code{\link{brokenstick}}, \code{\link{optBrokenstick}} 
#' @keywords internal brokenstick
#' @export
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj= exses)[[1]]
#' bsm <- brokenstick(dv) 
#' plot(depth ~ time, dv, ylim = rev(range(dv$depth)), type = 'l')
#' plot(bsm, add = TRUE, enumerate = TRUE)
plot.bsm <- function(x, type = 'b', lwd = 2, ylim = rev(range(y)), 
                     add = FALSE, col = add+1, enumerate = FALSE, ...) {
  y <- predict(x, newdata = x$pts.x)
  xy <- xy.coords(x$pts.x, y, 'x', 'y')
  if (add) lines(xy, type = type, lwd = lwd, col = col, ...)
  else plot(xy, type = type, lwd = lwd, ylim = ylim, col = col, ...)
  if (enumerate) text(xy, labels = x$pts.no, adj = c(1.5, 1.5), col = col, cex = .8)
}

#' Extract brokenstick model residuals
#' 
#' @param object \code{bsm} object, typically result from \code{\link{brokenstick}} 
#' or \code{\link{optBrokenstick}}.
#' @param type To choose in \code{c('normal', 'absolute')}. The second choice 
#' returning the absolute value of the first choice output.
#' @param newdata \code{x} and \code{y} data to used when interested in specific 
#' residual values (\code{x} in first column, \code{y} in second column).
#' @param ... Other arguments.
#' @return Returns model residuals if \code{newdata} is omited. Returns residuals 
#' computed from \code{newdata} otherwise. In this case, returns \code{NAs} for 
#' values that do not belong to a stick .
#' @details Require \code{eco.mem <= 4} if using \code{newdata} argument but requires 
#' \code{eco.mem <= 2} otherwise (see \code{\link{bsmfit}}).
#' @seealso \code{\link{brokenstick}}, \code{\link{optBrokenstick}}
#' @export
#' @keywords internal brokenstick
#' @examples
#' data(exses)
#' dv <- dv <- tdrply(identity, 1:2, no = 400, obj= exses)[[1]]
#' bsm <- brokenstick(dv)
#' plot(residuals(bsm)) ; abline(v = bsm$pts, h = 0)
residuals.bsm <- function(object, type = c('normal', 'absolute'), newdata, ...) {
	if (!missing(newdata)) {
		if (length(newdata) < 2) {stop('Please provide "x" and "y" data in "newdata" argument.')}
		ypred <-  predict.bsm(object, newdata = newdata[ , 1])
		out <- newdata[ , 2] - ypred
	} else {
		if (all(!'residuals' %in% names(object))) {
			stop('Cant return residuals without "residuals" slot or "newdata" argument.')
		} else {
			out <- object$residuals
		}
	}
	switch(match.arg(type), normal = out, absolute = abs(out))
}


#' Fitting automatic brokenstick models
#' 
#' \code{optBrokenstick} is similar to \code{\link{brokenstick}} except that a 
#' cost function can be used to determine the optimal number of points.
#' 
#' @inheritParams brokenstick
#' @param x The x data. Note that the \code{x} values have to be sorted. Can also 
#' be a list of \code{x} and \code{y} data as returned by \code{\link{xy.coords}}.
#' @param y The y data.
#' @param threshold A threshold value for the cost function to be used instead of 
#' the minimum. If provided the search of a local minimum in the cost function is 
#' abandoned.
#' @param cost The cost function to use. Two included in the package 
#' \code{\link{dist_per_pt_cost}} and \code{\link{max_dist_cost}}. Feel free to use 
#' a custom one.
#' @param npmin Minimun number of points.
#' @param npmax Maximum number of points.
#' @return Same as \code{\link{brokenstick}} with the value of the cost function.
#' @seealso \code{\link{brokenstick}} and \code{\link{predict.bsm}}, \code{\link{residuals.bsm}}, 
#' \code{\link{update.bsm}}, \code{\link{summary.bsm}}, 
#' \code{\link{coef.bsm}}, \code{\link{plot.bsm}}, \code{\link{as.data.frame.bsm}}
#' for other functions with a S3 method for \code{bsm} objects. 
#' @export
#' @keywords brokenstick
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' plot(depth ~ time, dv, ylim = rev(range(dv$depth)), type = 'l')
#' bsm <- optBrokenstick(dv)
#' plot(bsm, add = TRUE)
#' bsm <- optBrokenstick(dv, threshold = 20, cost = max_dist_cost) 
#' plot(bsm, add = TRUE, col = 'blue', enumerate = TRUE)
optBrokenstick <- function(x, y = NULL, threshold, cost = max_dist_cost, 
						   npmin = 2, npmax = Inf, start = NULL, na.action, ...) {
	if (npmin < length(start)) 
		stop('"npmin" must be equal or greater than the length of "start".')
	if (npmax <= length(start)) 
		stop('"npmax" must be greater than the length of "start".')
	if (npmax <= npmin) 
		stop('"npmax" must be greater than "npmin".')
	if (missing(na.action)) na.action <- options("na.action")[[1]]
	if (missing(threshold)) threshold <- -Inf
	bsm0 <- brokenstick(x, y, npmin, start, na.action, ...)
	S0 <- cost(bsm0)
	repeat {
		if (S0 <= threshold || length(bsm0$pts) == npmax) break
		bsm <- update(bsm0, npts = length(bsm0$pts) + 1)
		S <- cost(bsm)
		if (is.infinite(threshold) && S >= S0) break
		bsm0 <- bsm ; S0 <- S
	}
	bsm <- update(bsm0, npts = length(bsm0$pts), ...)
	bsm$cost <- cost(bsm)
	bsm
}

#' Cost functions for automatic brokenstick models.
#' 
#' Given a model and data the function returns a single value which is used to 
#' determine if the number of points is optimized: the cost function value being 
#' minimized by the optimal set of parameters.
#'
#' @param object Object of class inheriting from "\code{bsm}".
#' @return A statistic to minimize. 
#' @details \code{max_dist_cost} In this function the statistic is the maximun distance between an 
#' observation and its fitted value. Hence the function constantly decrease with 
#' increasing number of point in the model and a threshold must be provided along 
#' with this function to avoid infinite looping. Yet, the pros of this cost function 
#' is a threshold that is simple to determine and to interpret.
#' @seealso \code{\link{optBrokenstick}}
#' @keywords internal brokenstick
#' @export
max_dist_cost <- function(object) {
	max(residuals.bsm(object, type = 'absolute'))
}

#' @rdname max_dist_cost
#' @inheritParams max_dist_cost
#' @details \code{dist_per_pt_cost} In this function the statistic is the 
#' average distance between observations and fitted values divided by the number of
#' points used by the model. This function has does rech a minimun value but 
#' generally for high numbers of points.
#' @keywords internal brokenstick
dist_per_pt_cost <- function(object) {
	res <- residuals.bsm(object, type = 'absolute')
	mean(res, na.rm = TRUE) / length(object$pts)
}

#' Coerce brokenstick model to data.frame
#' 
#' @param x a brokenstick model.
#' @inheritParams base::as.data.frame
#' @export
#' @keywords internal brokenstick
#' @examples 
#' data(exses)
#' bsm <- tdrply(brokenstick, 1:2, no = 50:53, obj = exses)
#' lapply(bsm, as.data.frame)
as.data.frame.bsm <- function(x, row.names = NULL, ...) {
  n <- length(x$pts)
  df <- data.frame(st_tm = x$pts.x[-n], ed_tm = x$pts.x[-1], no_seg = seq(1, n-1),
             bsm_slope = x$slope, intercept = x$intercept)
  as.data.frame(df, row.names = row.names, ...)
}

#' Print summary of a brokenstick model
#' 
#' @param object a brokenstick model
#' @param ... for method consistency
#' @export
#' @keywords internal brokenstick
#' @examples 
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 53, obj = exses)[[1]]
#' summary(brokenstick(dv))
#' summary(brokenstick(dv, npts = 12))
summary.bsm <- function(object, ...) {
  df <- as.data.frame(object)
  df$dur <- apply(df[ , 1:2], 1, function(object) diff(as.numeric(object)))
  print(df[ , -(1:2)])
  r2 <- 1 - (var(resid(object)) / var(object$data[ , 2]))
  cat("\nR-squared =", r2, "\nMax residual =", maxr <- max(resid(object)), 
      "\nMean squared resiluals =", meanr <- mean(resid(object)^2))
  invisible(list(df = df[ , -(1:2)], r2 = r2, max_res = maxr, mean_res = meanr))
}
