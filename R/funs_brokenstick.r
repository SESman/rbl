#' Fitter function for brokenstick models
#' 
#' Basic computing engines of the brokenstick algorithm called by 
#' \code{\link{brokenstick}}.
#' 
#' @param xy Data, as returned by \code{\link{xy.coords}}. The \code{x} values 
#' have to be sorted.
#' @param pts The points to use.
#' @param eco.mem An integer value between 0 and 7 to control the memory size of the 
#' output. \code{0} = no memory savings.
#' @seealso \code{\link{brokenstick}} and \code{\link{optBrokenstick}} which should 
#' be used to fit brokenstick models.
#' @keywords internal
#' @export
bsmfit <- function(xy, pts, eco.mem = 0L) { 
  if (eco.mem > 4) stop('"eco.mem" must be between 0 and 4.')
  
  # Compute stick slopes and coefficients
  n <- diff(pts)
  n[1L] <- n[1L] + 1 # First segment needs 1 more prediction for the first point
  dy <- diff(xy[pts, 2L]) ; dx <- diff(xy[pts, 1L])
  a <- dy / dx
  b <- xy[pts[-length(pts)], 2L] - (a * xy[pts[-length(pts)], 1L])
  
  # Compute fitted values and residuals
  fit <- rep(a, n) * xy[ , 1L] + rep(b, n)
  res <- xy[, 2] - fit
  
  # Format output
  out <- list(pts.x = xy[pts, 1L], pts.y = xy[pts, 2L], slope = a, intercept = b, 
              pts = pts, residuals = res, fitted.values = fit, data = xy)
  "if"(eco.mem != 0L, out[-unique(seq(8L - eco.mem + 1L, 8L))], out)
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
#' @return A \code{bsm} object with (depending on \code{eco.mem}):
#' \itemize{
#'   \item pts.x The x values of the brokenstick points (\code{eco.mem} inefficient).
#'   \item pts.y The y values of the brokenstick points (\code{eco.mem} inefficient).
#'   \item slope The slopes of each "stick" of the model (\code{eco.mem} inefficient).
#'   \item intercept The intercept of each "stick" of the model (\code{eco.mem} inefficient).
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
brokenstick.default <- function(x, y = NULL, npts = 6, start = NULL, na.action, ...) {
  # Format input data
  nms <- if (is.recursive(x)) {names(x)} else {c(deparse(substitute(x)), deparse(substitute(y)))}
  if (any(grepl('\\$', nms))) nms <- gsub('(.*\\$)(.*$)', '\\2', nms)
  else if (any(sapply(nms, nchar) > 10)) nms <- c('x', 'y')
  xy <- setNames(as.data.frame(xy.coords(x, y)[1:2]), nms) # formated data
  
  # Apply na.action
  if (missing(na.action)) na.action <- options("na.action")[[1]]
  xy <- do.call(na.action, list(xy))
  
  # Broken sticks algorithm
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
  
  # Format output
  out$na.action <- attr(xy, "na.action")
  class(out) <- c('bsm', 'list')
  out
}

#' @inheritParams brokenstick
#' @keywords internal
#' @export
brokenstick.formula <- function(x, y = NULL, npts = 6, start = NULL, 
                                na.action, ...) {
  # Format input data from "formula" ("x" arg) and "data" ("y" arg) syntax
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("x", "y"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  names(mf) <- c('', 'formula', 'data')[seq_along(mf)]
  mf <- eval(mf, parent.frame())
  
  # Apply na.action
  if (missing(na.action)) na.action <- options("na.action")[[1]]
  
  # Launch algorithm as usual
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
update.bsm <- function(object, npts, allow.dup = FALSE, ...) {
  # Check if an update is needed
  np <- length(object$pts.x)
  if (np == npts) return(object)
  
  # If Hi-Res data are not available update possible with smaller number of break points
  if (!"data" %in% names(object)) {
    if (npts > np) {
      stop('Cannot improve BSM resolution without data (set "eco.mem" to 0).')
    } else {
      # Just select the first npts breakpoints but cannot get fitted values and residuals
      return(brokenstick(object$pts.x, object$pts.y, npts, eco.mem = 4))
    }
  }
  
  # If Hi-Res data are available...
  if (np < npts) {
    # Resume process at last iteration and go on to requested number of break points
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
        # Stop if duplicated break points not allowed ("allow.dup" arg)
        if (!allow.dup) {
          pts <- pts[!dup] ; pts.no <- pts.no[!dup]
          break
        }
      }
      np <- length(pts)
    }
  } else {
    # Just select the first npts breakpoints
    cnd <- object$pts.no <= (npts - 1)
    pts <- object$pts[cnd]
    pts.no <- object$pts.no[cnd]
  }
  
  # Refit BSM to get the correct coefficients, fitted values and residuals
  x <- bsmfit(object$data, pts, ...)
  x$pts.no <- pts.no
  
  # Format output
  x$na.action <- "if"("na.action" %in% names(object), object$na.action, NULL)
  class(x) <- c("bsm", "list")
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
  # If "newdata" not provided, return fitted values else make new prediction
  if (missing(newdata)) {
    if (is.null(fitted(object))) {
      if (!"data" %in% names(object)) {
        stop("Impossible without x values. ", 'Consider add theses or set "eco.mem" less than 2.')
      } else {
        # Preventive programming, this case should not occur:
        # "data" slot is missing before "fitted" slot when "eco.mem" argument is used
        n <- diff(object$pts)
        n[1] <- n[1] + 1 
        return(rep(object$slope, n) * object$data$x + rep(object$intercept, n))
      }
    } else {
      return(fitted(object))
    }
  } else {
    # Find out which coefficients must be applied
    stks <- which.stick(object, newdata)
    if (is.POSIXct(newdata)){newdata <- as.numeric(newdata)}
    return(object$slope[stks] * newdata + object$intercept[stks])
  }
  
  # Preventive programming, for bug reports
  stop('Unexpected case occured in "predict.bsm".')
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
which.stick <- function(object, pts, type = c("x", "i")) {
  bsm.pts <- switch(match.arg(type), x = object$pts.x, i = object$pts)
  eql <- lapply(pts, function(x) x == bsm.pts)
  grt <- lapply(pts, function(x) max(which(x >= bsm.pts) %else% NA))
  lst <- lapply(pts, function(x) min(which(x <= bsm.pts) %else% NA))
  
  # Function returning scitck number given the location of data in comparison to 
  # break points
  .f <- function(eql, grt, lst) {
    if (any(eql)) {ifelse(which(eql) == length(eql), which(eql) - 1, which(eql))}
    else if (is.na(grt) || is.na(lst)) {NA}
    else {grt}
  }
  as.integer(mapply(.f, eql, grt, lst))
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
  # Just format corresponding slots into a data.frame output
  out <- data.frame(intercept = object$intercept, slope = object$slope)
  row.names(out) <- paste0("seg", seq(nrow(out)))
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
#' @param data Should the data used to fit the BSM be plotted as well ?
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
plot.bsm <- function(x, type = "b", lwd = 2, ylim = rev(range(y)), 
                     add = FALSE, col = (add || data) + 1, 
                     enumerate = FALSE, data = FALSE, ...) {
  # Generate BSM abstracted profile
  y <- predict(x, newdata = x$pts.x)
  xy <- xy.coords(x$pts.x, y, "x", "y")
  
  # Draw it
  if (add) lines(xy, type = type, lwd = lwd, col = col, ...)
  else plot(xy, type = type, lwd = lwd, ylim = ylim, col = col, ...)
  
  # If break points numbering is requested
  if (enumerate) text(xy, labels = x$pts.no, adj = c(1.5, 1.5), col = col, cex = .8)
  
  # Hi-Res data can be added if requested and available
  if (data) {
    "data" %in% names(x) || stop('"data" slot is missing in "x".')
    lines(x$data)
  }
  
  invisible(NULL)
}

#' Extract brokenstick model residuals
#' 
#' @param object \code{bsm} object, typically result from \code{\link{brokenstick}} 
#' or \code{\link{optBrokenstick}}.
#' @param type To choose in \code{c('normal', 'absolute')}. The second choice 
#' returning the absolute value of the first.
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
#' dv <- tdrply(identity, 1:2, no = 400, obj= exses)[[1]]
#' bsm <- brokenstick(dv)
#' plot(residuals(bsm)) ; abline(v = bsm$pts, h = 0)
residuals.bsm <- function(object, type = c("normal", "absolute"), newdata, ...) {
  # If "newdata" not provided, return residual slot else make new prediction
  if (!missing(newdata)) {
    if (length(newdata) < 2) {stop('Please provide "x" and "y" data in "newdata" argument.')}
    ypred <-  predict.bsm(object, newdata = newdata[ , 1])
    out <- newdata[ , 2] - ypred
  } else {
    if (all(!"residuals" %in% names(object))) {
      stop('Cant return residuals without "residuals" slot or "newdata" argument.')
    } else {
      out <- object$residuals
    }
  }
  
  # Format output to requested "type"
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
  # Check inputs consitency
  if (npmin < length(start)) 
    stop('"npmin" must be equal or greater than the length of "start".')
  if (npmax <= length(start)) 
    stop('"npmax" must be greater than the length of "start".')
  if (npmax <= npmin) 
    stop('"npmax" must be greater than "npmin".')
  
  # Set defaults for "na.action" and "threshold"
  if (missing(na.action)) na.action <- options("na.action")[[1]]
  if (missing(threshold)) threshold <- -Inf
  
  # Initiate algorithm
  bsm0 <- brokenstick(x, y, npmin, start, na.action, ...)
  S0 <- cost(bsm0)
  # Start iterations
  repeat {
    # Stop if max number of break point or if is the threshold cost is reached
    if (S0 <= threshold || length(bsm0$pts) == npmax) break
    bsm <- update(bsm0, npts = length(bsm0$pts) + 1)
    S <- cost(bsm)
    # Stop if threshold does not make sense or cost increases
    if (is.infinite(threshold) && S >= S0) break
    bsm0 <- bsm ; S0 <- S
  }
  bsm <- update(bsm0, npts = length(bsm0$pts), ...)
  
  # Format output
  bsm$cost <- cost(bsm)
  bsm
}

#' Get the maximun residual of a BSM at a given iteration 
#' 
#' @param x \code{bsm} object as returned by \code{\link{brokenstick}} or 
#' \code{\link{optBrokenstick}}.
#' @param iter the iteration number for which the maximum residual is to be 
#' returned (Ri). If \code{NULL} then \code{iter} is set to the last iteration 
#' in \code{x}
#' @inheritParams residuals.bsm
#' @export
#' @keywords internal brokenstick
#' @examples 
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' bsm <- brokenstick(dv) 
#' max_residual(bsm, 5)
#' 
#' \dontrun{
#' bsm <- brokenstick(dv, eco.mem = 4) 
#' max_residual(bsm, 5) # error
#' max_residual(bsm, 4)
#' }
max_residual <- function(x, iter = NULL, type = c("normal", "absolute")) {
  # Check inputs and set iter to default value when necessary
  is.bsm(x) || stop('x must be of class "bsm"')
  iter <- iter %else% (length(x$pts.x) - 1)
  
  # Update "x" to necessary number of break point given "iter"
  npts <- iter + 1
  bsm0 <- update(x, npts)
  
  # If data slot is available use it so that any iteration can be asked
  if ("data" %in% names(x)) {
    res <- resid(bsm0)
    out <- res[which.max(abs(res))]
  } else {
    bsm1 <- try(update(x, npts + 1), TRUE) %else% stop('Requested iteration ', 
                                                       'number "iter" is too high.')
    rk <- which.max(bsm1$pts.no)
    out <- predict(bsm0, newdata = bsm1$pts.x[rk]) - bsm1$pts.y[rk]
  }
  
  # Format output to requested "type"
  switch(match.arg(type), normal = out, absolute = abs(out))
}

#' Compute goodness of fit of brokenstick models: Dive Zone Index (DZI)
#' 
#' An index of the goodness of fit for brokenstick models. Value of the dive zone 
#' index ranges from 0 (perfect fit) to 1. See original article in references.
#' 
#' @param x a \code{bsm} object as returned by \code{\link{brokenstick}} or 
#' \code{\link{optBrokenstick}}.
#' @param iter the iteration number for which the DZI is to be computed. If 
#' \code{iter = NULL} then it is set to last BSM iteration.
#' @param n Optional. The number of points to use when calculating the dive zone limits. 
#' if \code{NULL} then \code{n} is set to the number of record in the dataset used to 
#' fit the BSM (when available) or 500 (when not available)
#' @export
#' @keywords brokenstick
#' @details If the original TDR data are not available to the function the DZI 
#' of the last BSM iteration will be computed using the maximum residual from the 
#' previous iteration.
#' @return A \code{dzi} object with:
#' \itemize{
#'   \item dzi The dive zone index obtained at each iteration.
#'   \item max_res The maximum residuals at each iteration.
#'   \item dz_Lbnd Dive zone lower bound.
#'   \item dz_Ubnd Dive zone upper bound.
#'   \item dz_width The difference between the two previous slots i.e. the vertical 
#'   width of the dive zone.
#'   \item dz_Xval A vector giving the x values corresponding to \code{dz_Lbnd}, 
#'   \code{dz_Ubnd} and \code{dz_width}.
#'   \item no_seg A vector giving the number of BSM segments to which the 
#'   \code{dz_Xval} belong.
#'   \item seg_width Vertical width covered by BSM segments.
#'   \item seq_length Duration of BSM segments.
#'   \item pts.x, pts.y, pts.no Breakpoints information inherited from \code{x}.
#'   \item data The raw data of input BSM when available.
#' }
#' @references
#' Photopoulou, T., Lovell, P., Fedak, M. A., Thomas, L. and Matthiopoulos, J. (2015). 
#' Efficient abstracting of dive profiles using a broken-stick model. 
#' Methods Ecol Evol 6, 278-288. Github repo: https://github.com/theoniphotopoulou/brokenstickmodel.git
#' @examples 
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' bsm <- brokenstick(dv) 
#' dzi <- dive_zone_index(bsm)
#' 
#' dzi
#' str(dzi)
#' plot(dzi)
dive_zone_index <- function(x, iter = NULL, n = NULL) {
  is.bsm(x) || stop('x must be of class "bsm"')
  
  # Compute constants
  npts_ini <- max(x$pts.no)
  npts_ini > 1 || stop('"x" must have at least 3 break points.')
  iter <- iter %else% npts_ini
  rng_y <- range(x$pts.y)
  if (has_data <- ("data" %in% names(x))) {
    n <- nrow(x$data)
    newx <- x$data[ , 1]
  } else {
    if (iter >= npts_ini) warning('Maximum residual (Ri) is not computable for the requested iteration number "iter".', 
                                  'The last Ri available will be used instead.')
    n <- n %else% 500
    newx <- seq(min(x$pts.x), max(x$pts.x), length.out = n)
  }
  
  # Initiate self-dependent variables
  R <- dz_index <- NULL
  Ubnd <- -(Lbnd <- rep(Inf, n))
  
  for (i in seq(1, iter)) {
    bsm_i <- update(x, npts = i + 1) # BSM at iteration i
    bsm_i_pred <- predict(bsm_i, newx)
    R <- c(R , try(max_residual(x, i, type = "absolute"), TRUE) %else% NULL)
    Ri <- R[length(R)]
    
    # Compute Dive Zone limits and the corresponding index
    Ubnd <- pmax(bsm_i_pred - Ri, rng_y[1], Ubnd)
    Lbnd <- pmin(bsm_i_pred + Ri, rng_y[2], Lbnd)
    dz_width <- Lbnd - Ubnd
    dz_index <- c(dz_index, sum(dz_width) / (diff(rng_y) * n))
  }
  
  # Format output
  out <- list(
    dzi = setNames(dz_index, paste0("DZI", seq(1, iter))), 
    max_res = setNames(R, paste0("R", seq(1, length(R)))), 
    dz_Lbnd = Lbnd, dz_Ubnd = Ubnd, dz_width = dz_width, 
    dz_Xval = newx, no_seg = which.stick(bsm_i, newx), 
    seg_width = abs(diff(bsm_i$pts.y)), seg_length = abs(diff(bsm_i$pts.x)), 
    pts.x = x$pts.x, pts.y = x$pts.y, pts.no = x$pts.no, 
    data = "if"(has_data, x$data, NULL)
  )
  class(out) <- c("dzi", "list")
  out
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
  max_residual(object, type = 'absolute')
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

#' @rdname max_dist_cost
#' @inheritParams max_dist_cost
#' @details \code{rss_cost} In this function the statistic is the sum of 
#' squared residuals.
#' @keywords internal brokenstick
rss_cost <- function(object) {
  sum(residuals.bsm(object)^2)
}

#' @rdname max_dist_cost
#' @inheritParams max_dist_cost
#' @details \code{dzi_cost} In this function the statistic is the Dive Zone Index. 
#' See details in \code{\link{dive_zone_index}}
#' @keywords internal brokenstick
dzi_cost <- function(object) {
  if (max(object$pts.no) == 1) stop('"object" must have at least 3 break points. ', 
                                    'Set "npmin" to 3 in "optBrokenstick"')
  dive_zone_index(object)
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
  n <- length(x$pts.x)
  df <- data.frame(st_tm = x$pts.x[-n], ed_tm = x$pts.x[-1], no_seg = seq(1, n-1),
                   bsm_slope = x$slope, intercept = x$intercept, duration = diff(x$pts.x))
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

#' General S3 utils for bsm objects
#' 
#' @param pts.x x coordinates of breakpoints. If \code{pts.x} is a list then 
#' it is interpreted as being a \code{"bsm"} object and is returned as is with 
#' updated class.
#' @param pts.y y coordinates of breakpoints
#' @param ... other BSM slots such as \code{"data"}, \code{"fitted"} or 
#' \code{"residuals"}. See details of slots in \code{\link{brokenstick}}.
#' @export
#' @keywords internal brokenstick
#' @examples
#' data(exses)
#' dv <- tdrply(identity, 1:2, no = 400, obj = exses)[[1]]
#' bsm <- brokenstick(dv) 
#' 
#' as.bsm(bsm$pts.x, bsm$pts.y)
#' as.bsm(bsm$pts.x, bsm$pts.y, residuals = "dummy")
#' as.bsm(bsm)
as.bsm <- function(pts.x, pts.y = NULL, ...) {
  if (is.bsm(pts.x) || is.list(pts.x)) {
    class(pts.x) <- c("bsm", "list")
    return(pts.x)
  } else {
    slts <- list(...)
    out <- brokenstick(pts.x, pts.y, length(pts.x), eco.mem = 4)
    out[names(slts)] <- slts
  }
  out
}

#' @rdname as.bsm
#' @param x an object to test.
#' @export
#' @keywords internal brokenstick
#' @examples 
#' is.bsm(bsm)
is.bsm <- function(x) is(x, "bsm")

#' Print method for bsm objects
#' @param x a bsm object
#' @param ... for generic compatibility
#' @export
#' @keywords internal brokenstick
print.bsm <- function(x, ...) print(as.data.frame(x))

#' Print method for dzi objects
#' @param x a dzi object
#' @param ... for generic compatibility
#' @export
#' @keywords internal brokenstick
print.dzi <- function(x, ...) print(x$dzi)

#' Plot method for "dzi" objects
#' @param x a \code{"dzi"} object
#' @param dz_col color of dive zone area
#' @param dz_border color of the border of the dive zone area
#' @param enumerate should the order of BSM break points be enumerated.
#' @param ... Arguments to be passed to methods, such as graphical 
#' parameters (see \code{\link{par}}).
#' @export
#' @keywords brokenstick
plot.dzi <- function(x, dz_col = "lightblue", dz_border = "blue", enumerate = TRUE, ...) {
  has_data <- ("data" %in% names(x))
  iter <- length(x$dzi)
  
  # Make plot drawing area
  if (has_data) 
    plot(x$data, type = "n", ylim = rev(range(x$data[ , 2])), ...)
  else 
    plot(x$pts.x, x$pts.y, type = "n", ylim = rev(range(x$pts.y)), ...)
  
  # Plot dive zone
  polygon(c(x$dz_Xval, rev(x$dz_Xval)), c(x$dz_Lbnd, rev(x$dz_Ubnd)), 
          col = dz_col, border = dz_border)
  
  # Add dive profiles (hi-res if available + abstracted profile)
  if (has_data) lines(x$data)
  plot(brokenstick(x$pts.x, x$pts.y, npts = iter + 1), enumerate = enumerate, add = TRUE)
  
  # Add title
  title(paste0("Iteration ", iter, " (", iter + 1, " brkpts): ", 
               "R", length(x$max_res), " = ", round(x$max_res[length(x$max_res)], digits = 2), ", ", 
               "DZI", iter, " = ", round(x$dzi[iter], digits = 2)))
  
  invisible(NULL)
}
