#' Count/Extract wiggles in 2D dataset.
#' 
#' According to Halsey et al. 2007 (see references): Wiggles are a particular 
#' pattern in the dive profile over time during a dive where an increase in depth 
#' over time changes to a decrease in depth and then back to an increase in 
#' depth. This creates a short period in the dive profile that is concave in 
#' shape. Wiggles are defined as elements of the dive profile during which at 
#' three points the vertical speed passes below 0 ms–1. (NB: If useful, 
#' certain wiggles could be ignored, e.g. using a threshold based on their 
#' depth range or duration.). This function implements this definition of the 
#' wiggle and is primarily intended to be used on dive profile but it can be 
#' useful to extract the inversions in any kind of 2D data having a monotonous 
#' x variable.
#' 
#' @param x The x data. A monotonous variable such as the time sequence of a TDR 
#' dataset. Can also be a list of \code{x} and \code{y} 
#' (processed by \code{\link{xy.coords}}). If x is a data frame with motre than 
#' two columns the fist two columns are used.
#' @param y The y data.
#' @param thres.y minimum y difference within a wiggle for it to be taken into 
#' account. The default values are usually appropriate if Y is a depth variable 
#' from a southern elephant seal dataset. If one value is provided wiggles are kept 
#' if the y differences are greater than this threshold. If two values are 
#' provided wiggles are kept when the y differences lie between these thresholds. 
#' This threshold is used by \code{\link{optBrokenstick}} in conjunction with 
#' \code{\link{max_dist_cost}} so it can not be set to \code{NULL} unless a value 
#' is passed to the \code{bsm} argument.
#' @param thres.x minimum x difference within a wiggle for it to be taken into 
#' account. The default value is usually appropriate if X is a time variable 
#' from a southern elephant seal dataset. If one value is provided wiggles are 
#' kept when the x differences are greater than this threshold. If two values are 
#' provided wiggles are kept if the x differences lie between these thresholds. 
#' \code{NULL} is equivalent to \code{c(0, Inf)}.
#' @param plt Should graphics about processing be plotted ?
#' @param output Should the function return the number of wiggles ("wig-count"), 
#' the number of steps ("stp-count") or a data frame with width, height and 
#'  height/width ratio for each transit/step/wiggle identified ("table"). 
#' @param bsm To speed up the process you can provide a brokenstick model to use 
#' directly instead of computing a new one from x and y data.
#' @param step A vertical speed threshold defining "steps" (0.35 m/s for king 
#' penguin). If \code{NULL} then steps are ignored.
#' @param step.thres.x similar to \code{thres.x} but applies to steps only.
#' @param step.thres.y similar to \code{thres.y} but applies to steps only.
#' @export
#' @keywords behavior
#' @references Halsey, L.G., Bost, C.-A. & Handrich, Y. (2007) A thorough and 
#' quantified method for classifying seabird diving behaviour. 
#' Polar Biology, 30, 991–1004.
#' @examples
#' data(exses)
#' 
#' # Number of wiggles can be used as a proxy of the foraging activity
#' sunflowerplot(tdrply(wiggles, ty = '_', obj = exses), exses$stat$pca)
#' sunflowerplot(tdrply(wiggles, ty = '_', obj = exses, step = 0.35), exses$stat$pca)
#' 
#' # Identifying steps as well
#' tdrply(wiggles, c("time", "depth"), ty = '_', no = 65, obj = exses, 
#'        step = 0.35, output = "table", plt = TRUE)
wiggles <- function(x, y = NULL, thres.y = 2.5, thres.x = c(10, Inf), 
                    step = NULL, step.thres.y = NULL, step.thres.x = c(10, Inf), 
                    output = c("wig-count", "stp-count", "table"), 
                    plt = FALSE, bsm = NULL) {
  output <- match.arg(output, output)
  if (is.null(step) && output == "stp-count") 
    stop('"step" canot be NULL when step count is required as output.')
  xy <- as.data.frame(xy.coords(x, y)[1:2])
  if (is.null(bsm)) 
    bsm <- optBrokenstick(xy, threshold = thres.y, cost = max_dist_cost)
  
  # Check if there is no wiggle at all
  bsm_slp <- coef(bsm)$slope
  if (!is.null(step)) 
    bsm_slp[bsm_slp %bw% c(0, step)] <- 0
  slp <- per(sign(bsm_slp))
  
  # Function to check if successive rows of a "slp" table can be a wiggle
  check_wiggles <-function(slp) {
    # Wiggles must start with negative slope (or zero) (is toward surface if y = depth)
    if (slp$value[1] == 1) slp <- slp[-1, ]
    if (nrow(slp) <= 1) return(NULL)
    # Wiggles must end with positive slope (or zero) (is toward benthos if y = depth)
    if (slp$value[nrow(slp)] == -1) slp <- slp[-nrow(slp), ]
    if (nrow(slp) <= 1) return(NULL)
    slp
  }
  
  # Delimitate steps
  slp$no_stp <- ifelse(slp$value == 0, cumsum(slp$value == 0), 0)
  if (output == "stp-count" && all(slp$no_stp == 0)) return(0)
  slp$no_wig  <- 0
  # Delimitate wiggles between steps
  tmp <- per(slp$value == 0)
  tmp <- tmp[!tmp$value & tmp$length >= 2, ]
  wigs <- Map(function(st, ed) check_wiggles(slp[st:ed, ]), tmp$st_idx, tmp$ed_idx)
  # Give wiggles a number
  cnd <- !sapply(wigs, is.null)
  if (any(cnd)) {
    wigs <- wigs[cnd]
    n_wig <- sapply(wigs, nrow) / 2
    no_wig <- cumsum(n_wig)
    wigs <- Map(function(x, n, no) {
      x$no_wig <- rep(seq(no-n+1, no), each = 2)
      x}, wigs, n_wig, no_wig)
    # Update slp table
    for (ii in seq_along(wigs)) {
      slp[row.names(wigs[[ii]]), ] <- wigs[[ii]]
    }
  } else {
    if (output == "wig-count") return(0)
  }
  
  # Function to compute stats given a colum with id numbers
  pts_stck <- which.stick(bsm, bsm$pts, type = 'i')
  compute_stats <- function(x, by) {
    idx <- by(x, by, function(x) seq(min(x$st_idx), max(x$ed_idx)+1), 
              simplify = FALSE)
    idx <- idx[names(idx) %w/o% "0"]
    st <- sapply(idx, function(x) min(bsm$pts[pts_stck %in% x]))
    ed <- sapply(idx, function(x) max(bsm$pts[pts_stck %in% x]))
    # Comptute stats
    out <- data.frame(start_x = xy$x[st], end_x = xy$x[ed], no = seq_along(st))
    out$width <- mapply(function(st, ed) diff(xy$x[c(st, ed)]), st, ed)
    out$height <- mapply(function(st, ed) diff(range(xy$y[st:ed], na.rm = TRUE)), st, ed)
    out$ratio <- out$height / out$width
    out
  }
  
  # Get wiggles stats & check if match thres.x and thres.y conditions 
  if (any(slp$no_wig != 0)) {
    out_wig <- compute_stats(slp, slp$no_wig)
    out_wig$type <- "wiggle"
    if (length(thres.x) < 2) {
      cond_x <- out_wig$width >= (thres.x %else% 0)
    } else {
      cond_x <- out_wig$width %bw% (thres.x %else% c(0, Inf))
    }
    if (length(thres.y) != 2) {
      cond_y <- out_wig$height >= (thres.y %else% 0)
    } else {
      cond_y <- out_wig$height %bw% thres.y
    }
    cnd <- cond_x & cond_y
    no_rejected <- out_wig$no[!cnd]
    out_wig <- out_wig[cnd, ]
    out_wig$no <- seq_along(out_wig$start_x)
    slp$no_wig[slp$no_wig %in% no_rejected] <- 0
  } else {
    out_wig <- data.frame(start_x = numeric(), end_x = numeric(), no = numeric(), 
                          width = numeric(), height = numeric(), type =character())
  }
  
  # Get steps stats & check if match thres.x and thres.y conditions
  if (any(slp$no_stp != 0)) {
    out_stp <- compute_stats(slp, slp$no_stp)
    out_stp$type <- "step"
    if (length(step.thres.x) < 2) {
      cond_x <- out_stp$width >= (step.thres.x %else% 0)
    } else {
      cond_x <- out_stp$width %bw% (step.thres.x %else% c(0, Inf))
    }
    if (length(step.thres.y) != 2) {
      cond_y <- out_stp$height >= (step.thres.y %else% 0)
    } else {
      cond_y <- out_stp$height %bw% step.thres.y
    }
    cnd <- cond_x & cond_y
    no_rejected <- out_stp$no[!cnd]
    out_stp <- out_stp[cnd, ]
    out_stp$no <- seq_along(out_stp$start_x)
    slp$no_stp[slp$no_stp %in% no_rejected] <- 0
  } else {
    out_stp <- do.call(data.frame, setNames(lapply(out_wig, function(x) vector(mode(x))), names(out_wig)))
  }
  
  # Get transit stats
  # Update transit id with rejected wiggles and steps
  slp$no_trn <- (slp$no_wig == 0 & slp$no_stp == 0)
  slp$no_trn <- ifelse(slp$no_trn, cumsum(slp$no_trn), 0)
  if (any(slp$no_trn != 0)) {
    out_trn <- compute_stats(slp, slp$no_trn)
    out_trn$type <- "transit"
  } else {
    out_trn <- do.call(data.frame, setNames(lapply(out_wig, function(x) vector(mode(x))), names(out_wig)))
  }
  
  # Merge results
  out <- rbind(out_wig, out_stp, out_trn)
  out <- out[order(out$start_x), ]
  out$type <- factor(out$type, levels = c("transit", "wiggle", "step"))
  out <- na.omit(out)
  cnd <- duplicated(out$end_x)
  if (sum(cnd) > 1) browser()
  if (any(cnd)) out$end_x[which(cnd) - 1] <- out$start_x[which(cnd)]
  
  # Plot results
  if (plt) {
    yl <- rev(range(xy$y))
    plot(xy, type = 'l', ylim = yl)
    points(xy, col = out$type[which.bw(xy$x, out)], pch = 19)
    lims <- unique(c(out$start_x, out$end_x))
    abline(v = lims, lty = 2, col = 'gray')
    legend("topleft", legend = levels(out$type), col = 1:3, lwd = 2)
    ys <-  rnorm(nrow(out), max(xy$y), abs(diff(yl))*0.01)
    segments(x0 = out$start_x, x1 = out$end_x, y0 = ys, lwd = 2, col =  out$type)
    lims <- unique(c(out$start_x, out$end_x))
    abline(v = lims, lty = 2, col = 'gray')
    legend("topleft", legend = levels(out$type), col = 1:3, lwd = 2)
  }
  
  # Return result according to output argument
  if (output == "wig-count") {
    cnd <- out$type == "wiggle"
    out <- "if"(any(cnd), max(out$no[cnd]), 0)
  } else if (output == "stp-count") {
    cnd <- out$type == "step"
    out <- "if"(any(cnd), max(out$no[cnd]), 0)
  }
  out
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

#' Count the number of Prey Catch Attemps (PCA)
#' 
#' @param x a logical vector indicating at each timestamp if it was associated 
#' to a PCA event.
#' @details A continuous succession of \code{TRUE} is considered as a single PCA.
#' @export
#' @examples 
#' data(exses)
#' btt_pca <- tdrply(pca_count, "is_pca", ty = "_", obj = exses)
pca_count <- function(x) try(sum(per(x)$value)) %else% NA
