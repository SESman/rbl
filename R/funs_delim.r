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

#' Find drift dives using hight-resolution TDR data
#' 
#' \code{find_drift_tdr} takes an object and return a data frame with all the brokenstick 
#' segments along the trip with statistics such as their vertical speed (average 
#' and standard deviation), their duration and so on. This data frame is to be used to 
#' identify the drift segments manually (see examples for an example of method).
#' 
#' @param object a \code{SES} object including a TDR dataset, a "dives statistics" 
#' table and a "delim" table.
#' @param thres.bsm a threshold to use with \code{\link{max_dist_cost}} when 
#' running \code{\link{optBrokenstick}}.
#' @param bsm A list of brokenstick models to use directly instead of computing a 
#' new one with \code{thres.bsm} as parameter.
#' @export
#' @keywords drift
#' @examples
#' data(exses)
#' tab <- find_drift_tdr(exses)
#' 
#' \dontrun{
#' # The main criterion is the duration of segments. Then the R2 and the average vertical speed.
#' require(manipulate)
#' manipulate(
#' {
#'  cnd <- tab$bsm_r2 >= r2 & tab$dur >= dur * 60 & abs(tab$bsm_slope) <= slope
#'  cnd <- ifelse(is.na(cnd), FALSE, cnd)
#'  plot(bsm_slope ~ time, tab[cnd, ], xlim = range(tab$time), ylim = c(-0.8, 0.8))
#'  abline(h = 0, col = 'gray')
#' },
#'  r2 = slider(.95, 1, initial = .99, step = 0.0005, label = "bsm R^2"), 
#'  dur = slider(0, 30, initial = 5, step = 0.25, label = "Duration (minutes)"), 
#'  slope = slider(0, 4, initial = 0.6, step = 0.05, label = "Max. vertical speed (m/s)")
#' )
#' 
#' # Drift dive example
#' tdrply(plot, 1:2, no = 300, type = 'l', obj = exses)
#' }
find_drift_tdr <- function(object, thres.bsm = 5, bsm = NULL) {
  if (is.null(bsm))
    bsm <- tdrply(optBrokenstick, 1:2, cost = max_dist_cost, 
                  threshold = thres.bsm, obj = object)
  bsm_coef <- lapply(bsm, function(x) coef(x))
  n_seg    <- sapply(bsm_coef, nrow)
  coef_a   <- unlist(lapply(bsm_coef, function(x) x[ , 1]))
  coef_b   <- unlist(lapply(bsm_coef, function(x) x[ , 2]))
  
  bsm_rawidx <- df_search(lapply(bsm, function(x) x$pts))
  bsm_idx  <- df_search(Map(function(x, offset) x$pts + offset, bsm, 
                            offset = sapply(ty_delim(obj = object), function(x) x[1, 1]) - 1))
  seg_tab           <- as.data.frame(data.table::rbindlist(bsm_idx, use.names = FALSE))
  seg_tab$dur       <- apply(seg_tab, 1, diff) * time_reso(object$tdr[ , 1])
  seg_tab$bsm_slope <- coef_b
  seg_tab$no_dive   <- rep(numIn(bsm, TRUE), n_seg)
  seg_tab$time      <- object$tdr$time[seg_tab[ , 1]]
  seg_tab$min_depth <- mapply(function(st, ed) min(object$tdr[st:ed, 2]), 
                              seg_tab[ , 1], seg_tab[ , 2])
  seg_tab$max_depth <- mapply(function(st, ed) max(object$tdr[st:ed, 2]), 
                              seg_tab[ , 1], seg_tab[ , 2])
  seg_tab$bsm_msr   <- unlist(Map(function(x, df) Map(function(st, ed) mean((x$resid^2)[st:ed], na.rm = TRUE), 
                                                      st = df$V1, ed = df$V2), bsm, bsm_rawidx))
  
  bsm_var_res <- unlist(Map(function(x, df) Map(function(st, ed) var(x$resid[st:ed], na.rm = TRUE), 
                                            st = df$V1, ed = df$V2), bsm, bsm_rawidx))
  bsm_var_tot <- unlist(Map(function(x, df) Map(function(st, ed) var(x$data[st:ed, 2], na.rm = TRUE), 
                                                st = df$V1, ed = df$V2), bsm, bsm_rawidx))
  percent_res_var <- bsm_var_res / bsm_var_tot
  seg_tab$bsm_r2    <- 1 - ifelse(percent_res_var > 1, NA, percent_res_var)
  
  names(seg_tab) <- c('st_idx', 'ed_idx', names(seg_tab[ , -(1:2)]))
  seg_tab[ , c('st_idx','ed_idx', 'no_dive', 'time', 'bsm_slope', 
               'dur', 'bsm_msr', 'min_depth', 'max_depth', 'bsm_r2')]
}

#' Find drift dives using roll angle from static accelerometry
#' 
#' @param object a "ses" object
#' @param thres.roll threshold to apply to roll in order to discriminate 
#' "on belly" and "on back" attitude
#' @param thres.dur minimum duration for a continuous "on back" period to be 
#' considered as a drifting period.
#' @return A data frame one row per drift period with 
#' the start/end indices (refering to TDR rows) of the drift periods, 
#' dive number where the drifting period was observed, the starting time of the 
#' drift and the corresponding drift rate (m/s).
#' @export
#' @keywords drift
#' @examples 
#' data(exses)
#' tab <- find_drift_roll(exses)
#' plot(drift_rate ~ time, tab, ylim = c(-0.8, 0.8))
#' tdrply(plot, 1:2, no = sample(tab$no_dive, 1), obj = exses)
find_drift_roll <- function (object, thres.roll = pi/2, thres.dur = 400) { #...
  # Compute and center roll
  acc_cols <- c("axG", "ayG", "azG")
  roll <- atan2(object$tdr[ , acc_cols[2]], object$tdr[ , acc_cols[3]])
  roll <- roll - mean(roll, na.rm = TRUE)
  roll <-  atan2(sin(roll), cos(roll))
  
  # Apply roll threshold
  tmp <- per(abs(roll) >= thres.roll)
  
  # Correct test assuming animal is more often on belly than on the back
  ctrl <- diff(tapply(tmp$length, tmp$value, sum))
  if (ctrl > 0) tmp$value <- !tmp$value
  
  # Apply duration threshold
  tmp$value[tmp$length < thres.dur] <- FALSE
  tmp <- per(rep(tmp$value, tmp$length))
  tmp <- tmp[tmp$value, ]
  
  # Compute drift rate and pairing with dive numbers
  slps <- tdrply(function(x) coef(brokenstick(x, npts = 2))[1, 2], 1:2, 
                 ty = tmp[ , 1:2], obj = object)
  tm <- object$tdr[tmp$st_idx, 1]
  dvs <- which.dive(tm, object)
  data.frame(tmp[ , 1:2], no_dive = dvs, time = tm, drift_rate = slps)
}
