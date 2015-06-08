#' Sampling rate of a TDR dataset
#' 
#' @param x The time sequence.
#' @param n the number of replicates to do.
#' @param type Should the result be returned in seconds ('period') or hertz ('frequence') ?
#' @keywords internal
#' @export
time_reso <- function(x, n = 10, type = c('period', 'frequence')) {
  iis <- sample(1:length(x), n, replace = FALSE)
  dT <- round(unlist(lapply(Map(seq, iis, iis + 3), function(I) diff(x[I]))))
  reso  <- median(dT)
  if (is.na(reso) | !all(dT == reso))
    warning('Some missing time stamps were found. Check the dataset !')
  switch(match.arg(type), period = reso, frequence = 1 / reso)
}

#' Find the dives start and end
#' 
#' @param obj An object of class \code{'ses'}, \code{'statdives'} or \code{'tdr'}.
#' @param nms The names of the variables to use. Depends on the method.
#' @param thres.dur The minimum time (is seconds) a dive should last before it is 
#' indeed considered a dive.
#' @param thres.depth The depth threshold (in sensor units): diving 
#' period implies greater values.
#' @param warn Should a Warning column be added ?
#' @param ... Other arguments
#' @return A data frame with the following variable: indice of period start, 
#' indice of period end, type of period, duration (s), the dive number and a optional 
#' Warning column. Attributes of the data frame include some information about 
#' the processing.
#' @seealso \code{\link{bottom_delim}}
#' @export
#' @examples
#' data(exses)
#' dvs <- dive_delim(exses$tdr)
#' 
#' # Processing info
#' unique(dvs$warning)
#' attr(dvs, 'thres.depth')
#' attr(dvs, 'thres.dur')
dive_delim <- function(obj, nms = c('time', 'depth'), thres.dur = 300, thres.depth = 15, 
                       warn = TRUE, ...) {
  UseMethod('dive_delim')
}

#' @inheritParams dive_delim
#' @export
#' @keywords internal
dive_delim.default <- function(obj, nms = c('time', 'depth'), thres.dur = 300, 
                               thres.depth = 15, warn = TRUE, ...) {
  tm <- nms[1] ; dp <- nms[2]
  obj[ , tm] <- as.numeric(obj[ , tm])
  
  # ID diving/surface periods given a depth threshold
  dvs <- per(abs(obj[ , dp]) < abs(thres.depth))
  
  # Correction according to duration threshold
  # A dive has to last more than 'thres.dur' second
  dt <- round(obj[dvs$ed_idx, tm] - obj[dvs$st_idx, tm]) + 1
  cond <- !dvs$value & dt < thres.dur
  dvs$value[cond] <- TRUE
  out <- per(rep(dvs$value, times = dvs$length))
  
  # Add a Warning column to store possible problems
  if (warn) out$warn <- list(NA)
  
  # Check for missing time stamps
  dr <- (out$length - 1) * time_reso(obj[ , tm])
  dt <- round(obj[out$ed_idx, tm] - obj[out$st_idx, tm])
  missingtime <- dt != dr
  if (any(missingtime)) {
    warning('The datset has missing time stamps')
    if (warn) out$warn[missingtime] <- list('missing time stamp')
  }
  
  # Add a "no_dive" column
  out$no_dive <- cumsum(!out$value)
  out$no_dive <- ifelse(!out$value, out$no_dive, -out$no_dive)
  
  # Warning if last line is not a surface
  if (out$value[nrow(out)] == 'dive') {
    warning('The dataset ends with a dive. This last dive may be truncated.')
    if (warn) out$warn[nrow(out)] <- ifelse(is.na(out$warn[nrow(out)]), list('may be truncated'), 
                                            list(c(out$warn[nrow(out)], 'may be truncated')))
  }
  
  # Format & Set attributes with processing information
  out <- setNames(out, c("st_idx", "ed_idx", "type", "dive_dur", 'warning', "no_dive"))
  out_nms <- c('st_idx', 'ed_idx', 'no_dive', `if`(warn, 'warning', NULL))
  out <- out[ , out_nms]
  attr(out, 'thres.depth') <- thres.depth
  attr(out, 'thres.dur') <- thres.dur
  attr(out, 'ignored_dives') <- list(st_idx = dvs$st_idx[cond], dur = dur <- dt[cond], 
                                     hist = hist(dur, plot = FALSE))
  out
}

#' @inheritParams  dive_delim
#' @export
#' @keywords internal
dive_delim.ses <- function(obj, nms = c('time', 'depth'), thres.dur = 300, 
                           thres.depth = 15, warn = TRUE, ...) {
  if ('delim' %in% names(obj)) {
    dvs <- obj$delim
    warnCol <- `if`('warning' %in% dvs, which(names(dvs) == 'warning'), NULL)
    return(dvs[ , c(1:3, warnCol)])
  }
  dive_delim.default(obj, nms = c('time', 'depth'), thres.dur = 300, 
                     thres.depth = 15, warn = TRUE, ...)
}

#' Find the dives' bottom start and end
#' 
#' @param obj An object of class \code{'ses'}, \code{'statdives'} or \code{'tdr'}.
#' @param nms The names of the variables to use. Depends on the method.
#' @param thres.spd The vertical speed threshold.
#' @param w The window width for rolling mean applied to vertical speed.
#' @param min.depth The minimum depth (reported to the dive maximun depth) 
#' allowed for bottom limits.
#' @param ... Other arguments to be passed to \code{\link{dive_delim}}.
#' @export
#' @seealso \code{\link{dive_delim}}
#' @examples
#' data(exses)
#' dvs <- bottom_delim(exses$tdr)
bottom_delim <- function(obj, nms = c('time', 'depth'), thres.spd = 0.9, w = 12, 
                         min.depth = 0.4, ...) {
  UseMethod('bottom_delim')
}

#' @inheritParams bottom_delim
#' @export
#' @keywords internal
bottom_delim.default <- function(obj, nms = c('time', 'depth'), thres.spd = 0.9, w = 12, 
                                 min.depth = 0.4, ...) {
  dvs <- dive_delim(obj, nms, ...)
  dv_attr_nms <- c('thres.depth', 'thres.dur', 'ignored_dives')
  dv_attr <- attributes(dvs)[dv_attr_nms]
  tm <- nms[1] ; dp <- nms[2]
  obj[ , tm] <- as.numeric(obj[ , tm])
  obj <- as.ses(tdr = obj, delim = dvs)
  
  # Compute smoothed vertical speed
  dd <- c(NA, diff(obj$tdr[ , dp]))
  dt <- c(NA, diff(obj$tdr[ , tm]))
  obj$tdr$spd <- rollapply(dd / dt, mean, w, na.rm = TRUE)
  # Correct for "time jumps" > 10 s
  obj$tdr$spd <- ifelse(dt > 10, NA, obj$tdr$spd)
  
  # Bottom delim using 4 deg. poly.
  tm <- ifelse(is.numeric(tm), names(obj)[tm], tm) # Make sure tm is a character
  tmp <- dvs[isdv <- (dvs$no_dive > 0), ] 
  bttNms <- c("btt_st_idx", "btt_ed_idx")
  dvs[, bttNms] <- NA
  bttPoly <- function(x, offset) {
    # Simplify signal using a model (4 deg poly)
    mod <- lm(x[ , 2] ~ poly(x[ , 1], degree = 4, raw = TRUE))
    # Use threshold on predictions
    btt <- per(abs(predict(mod)) < abs(thres.spd))
    if (sum(btt$value) == 0)
      c(NA, NA)
    else 
      unlist(btt[btt$value, -(3:4)]) + rep(offset, 2)
  }
  tmp <- tdrply(bttPoly, c(tm, 'spd'), ty = tmp[ , 1:2], obj = obj, 
                la = list(offset = tmp$st_idx - 1))
  dvs[isdv, bttNms] <- Reduce(rbind, tmp)
  
  # Check that bottom limits are deeper than 'min.depth'*100 % of dive max depth
  thr <- abs(tdrply(max, dp, ty = dvs[ , 1:2], obj = obj) * min.depth)
  chkok <- abs(obj$tdr[dvs$btt_st_idx, dp]) >= thr & abs(obj$tdr[dvs$btt_ed_idx, dp]) >= thr
  chkok <- ifelse(is.na(chkok), TRUE, chkok)
  
  # Old fashion delim: idx of first and last depth >= 80% of Max depth
  tmp <- dvs[!chkok, ] ; nodv_old <- dvs$no_dive[which(!chkok)]
  thr <- (thr[!chkok] / min.depth) * 0.8 
  bttOld <- function(x, thr, offset) {
    rk <- which(abs(x) >= thr)
    if (length(rk) == 0)
      rk <- which.max(x)
    range(rk) + offset
  }
  tmp <- tdrply(bttOld, dp, ty = tmp[ , 1:2], obj = obj, 
                la = list(thr = thr, offset = tmp$st_idx -1))
  dvs[!chkok, bttNms[1:2]] <- Reduce(rbind, tmp)
  
  # Format & Set attributes
  out_nms <- c('st_idx', 'ed_idx', 'no_dive', 'btt_st_idx', 'btt_ed_idx', 
               `if`('warning' %in% names(dvs), 'warning', NULL))
  out <- dvs[ , out_nms]
  for (att in dv_attr_nms) {attr(out, att) <- dv_attr[[att]]}
  attr(out, 'old_bottom_delim') <- list(Dive.id = nodv_old)
  attr(out, 'thres.spd') <- thres.spd
  attr(out, 'w') <- w
  attr(out, 'min.depth') <- min.depth
  out
}

#' @inheritParams bottom_delim
#' @export
#' @keywords internal
bottom_delim.ses <- function(obj, nms = c('time', 'depth'), thres.spd = 0.9, w = 12, 
                             min.depth = 0.4, ...) {
  if ('delim' %in% names(obj))
    return(obj$delim)
  bottom_delim.default(obj, nms = c('time', 'depth'), thres.spd = 0.9, w = 12, 
                       min.depth = 0.4, ...)
}

#' Use brokensticks to identify potential drifts
#' 
#' \code{drift_stat} takes an object and return a data frame of all the brokenstick 
#' segments along the trip with statistics such as their vertical speed (average 
#' and standard deviation), their duration and so on. This data frame is to be 
#' filtered to identify drift segments.
#' 
#' @param object a \code{SES} object including a TDR dataset, a "dives statistics" 
#' table and a "delim" table.
#' @param thres.bsm a threshold to use with \code{\link{max_dist_cost}} when 
#' running \code{\link{optBrokenstick}}.
#' @param bsm A list of brokenstick models to use directly instead of computing a 
#' new one with \code{thres.bsm} as parameter.
#' @export
#' @return A data frame of all brokenstick segments with start and end index of 
#' the segment, no_seg_tot and ID number of the segement, no_seg_dive the number 
#' of the segment within its dive, no_dive the ID of the dive conataining the segment, 
#' time the starting time of segement, drift_rate the vertical speed of the segment, 
#' dur the duration of the segment, min_depth and max_depth the minimum and maximum 
#' depths of the segment, acceleration statistics (if acceleration data available) 
#' average roll, pitch and swimming effort and PCA rate, the mean squared residuals 
#' of the segment.
#' @details To get average pitch and roll angles (degrees), the TDR table must include 
#' static acceleration with variable names "axG", "ayG" and "azG". Similarly, for 
#' swimming effort and PCA rate, the TDR table must have columns "swm_eff" (numeric) 
#' and "is_pca" (logical). Use functions \code{\link{static_acc}}, 
#' \code{\link{swimming_effort}} and \code{\link{prey_catch_attempts}} to compute 
#' these variables from the raw acceleration readings.
#' @keywords drift
#' @examples
#' data(exses)
#' tab <- drift_stat(exses)
#' 
#' \dontrun{
#' require(manipulate)
#' manipulate(
#' {
#' tab$day <- floorPOSIXct(tab$time)
#'  cnd <- tab$dur >= dur * 60 & 
#'         abs(tab$drift_rate) <= rate & 
#'         tab$swm_eff <= swm_eff &
#'         abs(tab$pitch) <= pitch
#'  cnd <- "if"(no_pca, cnd & tab$pca_rate == 0, cnd)
#'  belly_up <- abs(tab$roll) >= 90
#'  if (up_only) cnd <- cnd & belly_up
#'  na_cnd <- is.na(cnd)
#'  cnd <- ifelse(is.na(cnd), TRUE, cnd)
#'  pts_col <- ifelse(belly_up[cnd], "red", "black")
#'  pts_col <- ifelse(na_cnd, "blue", pts_col)
#'  pts_shape <- ifelse(na_cnd, 20, 1)
#'  plot(drift_rate ~ time, tab[cnd, ], col = pts_col, pch = pts_shape, 
#'       xlim = range(tab$time), ylim = c(-3, 3))
#'  abline(h = 0, col = 'gray')
#' },
#'  dur = slider(0, 30, initial = 2, step = 0.25, label = "Duration (minutes)"), 
#'  rate = slider(0, 3, initial = 3, step = 0.05, label = "Max. vertical speed (m/s)"),
#'  swm_eff = slider(0, 2.5, initial = 2.5, step = .01, label = "Max. swimming effort (m/s³)"),
#'  pitch = slider(0, 90, initial = 90 , step = 5, label = "Max pitch angle (deg)"), 
#'  up_only = checkbox(FALSE, label = "Show belly up only"), 
#'  no_pca = checkbox(FALSE, label = "Show no PCA only")
#' )
#' 
#' # Drift dive example
#' tdrply(plot, 1:2, no = 300, type = 'l', obj = exses)
#' }
drift_stat <- function(object, thres.bsm = 5, bsm = NULL) {
  if (is.null(bsm))
    bsm <- tdrply(optBrokenstick, 1:2, cost = max_dist_cost, 
                  threshold = thres.bsm, obj = object)
  
  # Compute basic stats of each segment
  bsm_df <- as.data.frame(rbindlist(lapply(bsm, as.data.frame)))
  bsm_df$time <- as.POSIXct(bsm_df$st_tm, origin = "1970-01-01", tz = attr(object$tdr[ , 1], "tz"))
  bsm_df$dur <- apply(bsm_df[ , 1:2], 1, diff.default)
  bsm_df[ , 1:2] <- lapply(bsm_df[ , 1:2], which.row, obj = object)
  bsm_df$no_seg_tot <- seq_along(bsm_df[ , 1])
  bsm_df$no_dive <- which.dive(bsm_df$time, object)
  bsm_df$min_depth <- tdrply(min, 2, ty = bsm_df[ , 1:2], obj = object)
  bsm_df$max_depth <- tdrply(max, 2, ty = bsm_df[ , 1:2], obj = object)
  
  # Add "bsm fit" stats 
  bsm_res <- lapply(bsm, residuals)
  bsm_pts <- df_search(lapply(bsm, function(x) x$pts))
  mean_squared_residuals <- function(res, pts) {
    mapply(function(st, ed) mean(res[st:ed]^2), st = pts[ , 1], ed = pts[ ,2])
  }
  bsm_df$msr <- unlist(Map(mean_squared_residuals, bsm_res, bsm_pts))
  
  # If static acceleration available: add pitch and roll info
  nms_acc <- c()
  acc_cols <- c("ayG", "azG")
  if (all(acc_cols %in% names(object$tdr))) {
    roll <- atan2(object$tdr$ayG, object$tdr$azG)
    roll <- roll - agl_mean(roll)
    object$tdr$roll <- agl_rescale(roll)
    bsm_df$roll <- tdrply(agl_mean, "roll", ty = bsm_df[ , 1:2], obj = object) * 180/pi
    nms_acc <- "roll"
  }
  acc_cols <- c("axG", "ayG", "azG")
  if (all(acc_cols %in% names(object$tdr))) {
    object$tdr$pitch <- -atan(object$tdr$axG / sqrt(object$tdr$ayG^2 + object$tdr$azG^2))
    bsm_df$pitch <- tdrply(agl_mean, "pitch", ty = bsm_df[ , 1:2], obj = object) * 180/pi
    nms_acc <- c(nms_acc, "pitch")
  }
  # If swimming effort available: add swimming effort info
  if ("swm_eff" %in% names(object$tdr)) {
    bsm_df$swm_eff <- tdrply(mean, "swm_eff", ty = bsm_df[ , 1:2], obj = object)
    nms_acc <- c(nms_acc, "swm_eff")
  }
  # If prey catch attempts available: pca rate info
  if ("is_pca" %in% names(object$tdr)) {
    bsm_df$pca_rate <- tdrply(function(x) nrow(subset(per(x), value == TRUE)), 
                              "is_pca", ty = bsm_df[ , 1:2], obj = object) / bsm_df$dur
    nms_acc <- c(nms_acc, "pca_rate")
  }
  
  
  # Set names
  names(bsm_df) <- c("st_idx", "ed_idx", "no_seg_dive", "drift_rate", 
                     "intercept", "time", "dur", "no_seg_tot", "no_dive", 
                     "min_depth", "max_depth", "msr", nms_acc)
  # Reorder columns
  bsm_df[ , c("st_idx", "ed_idx", "no_seg_tot", "no_seg_dive", "no_dive", 
              "time", "drift_rate", "dur", "min_depth", "max_depth", 
              nms_acc, "msr")]
}

