#' Sampling rate of a TDR dataset
#' 
#' @param x The time sequence.
#' @param n the number of replicates to do.
#' @param type Should the result be returned in seconds ("period") or hertz ("frequence") ?
#' @keywords internal
#' @export
time_reso <- function(x, n = 10, type = c("period", "frequence")) {
  iis <- sample(1:length(x), n, replace = FALSE)
  dT <- round(unlist(lapply(Map(seq, iis, iis + 3), function(I) diff(x[I]))))
  reso  <- median(dT)
  switch(match.arg(type), period = reso, frequence = 1 / reso)
}

#' Find the dives start and end
#' 
#' @param obj An object of class \code{'ses'}, \code{'statdives'} or \code{'tdr'} 
#' with time and depth variables named "time" and "depth".
#' @param thres.dur The minimum time (seconds) a dive should last to be 
#' considered as a dive. If \code{thres.dur} has two elements the second is the 
#' maximum time (seconds) a dive should last to be considered as a dive. If 
#' \code{thres.dur = NULL} then this threshod is ignored.
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
#' table(dvs$warning)
#' attr(dvs, "thres.depth")
#' attr(dvs, "thres.dur")
dive_delim <- function(obj, thres.dur = c(300, 3000), thres.depth = 15, 
                       warn = TRUE, ...) {
  UseMethod("dive_delim")
}

#' @inheritParams dive_delim
#' @export
#' @keywords internal
dive_delim.default <- function(obj, thres.dur = c(300, 3000), thres.depth = 15, ...) {
  obj[ , "time"] <- as.numeric(obj[ , "time"])
  
  # ID diving/surface periods given a depth threshold
  dvs <- per(abs(obj[ , "depth"]) < abs(thres.depth))
  
  # Correction according to duration threshold
  # A dive has to last more than "thres.dur" second
  if (!is.null(thres.dur)) {
    dt <- round(obj[dvs$ed_idx, "time"] - obj[dvs$st_idx, "time"]) + 1
    cond_short <- !dvs$value & dt <= thres.dur[1]
    cond_long <- if (length(thres.dur) == 2) !dvs$value & dt >= thres.dur[2]
    else rep(FALSE, nrow(dvs))
    cond <- cond_short | cond_long
    dvs$value[cond] <- TRUE
    idx_short <- pmean(dvs$st_idx[cond_short], dvs$ed_idx[cond_short])
    idx_long <- pmean(dvs$st_idx[cond_long], dvs$ed_idx[cond_long])
    out <- per(rep(dvs$value, times = dvs$length))
  } else {
    cond <- rep(FALSE, nrow(dvs))
  }
  
  # Add a Warning column to store possible problems
  out$warn <- NA
  if (!is.null(thres.dur)) {
    out$warn[which.bw(idx_short, out[ , 1:2])] <- 'include a "short" dive.'
    out$warn[which.bw(idx_long, out[ , 1:2])] <- 'include a "long" dive.'
  }
  
  # Check for missing time stamps
  dr <- (out$length - 1) * time_reso(obj[ , "time"])
  dt <- round(obj[out$ed_idx, "time"] - obj[out$st_idx, "time"])
  missingtime <- dt != dr
  if (any(missingtime)) {
    warning("The datset has missing time stamps")
    out$warn[missingtime] <- upd_warn(out$warn[missingtime], "missing time stamp")
  }
  
  # Add a "no_dive" column
  out$no_dive <- cumsum(!out$value)
  out$no_dive <- ifelse(!out$value, out$no_dive, -out$no_dive)
  
  # Warning if last line is not a surface
  if (out$value[nrow(out)] == "dive") {
    warning("The dataset ends with a dive. This last dive may be truncated.")
    out$warn[nrow(out)] <- upd_warn(out$warn[nrow(out)], "may be truncated")
  }
  
  # Format & Set attributes with processing information
  out <- setNames(out, c("st_idx", "ed_idx", "type", "dive_dur", "warning", "no_dive"))
  out_nms <- c("st_idx", "ed_idx", "no_dive", "warning")
  out <- out[ , out_nms]
  attr(out, "thres.depth") <- thres.depth
  attr(out, "thres.dur") <- thres.dur
  attr(out, "ignored_dives") <- list(st_idx = dvs$st_idx[cond], dur = dur <- dt[cond])
  
  # Preallocate bottom limits columns
  bttNms <- c("btt_st_idx", "btt_ed_idx")
  out[, bttNms] <- NA
  out
}

#' @inheritParams  dive_delim
#' @export
#' @keywords internal
dive_delim.ses <- function(obj, thres.dur = c(300, 3000), thres.depth = 15, 
                           warn = TRUE, ...) {
  if ("delim" %in% names(obj)) {
    dvs <- obj$delim
    warnCol <- "if"("warning" %in% dvs, which(names(dvs) == "warning"), NULL)
    return(dvs[ , c(1:3, warnCol)])
  }
  dive_delim.default(obj, thres.dur, thres.depth, warn, ...)
}

#' Find the dives' bottom start and end
#' 
#' @param obj An object of class \code{'ses'} or \code{'tdr'}. In the first case 
#' the function just returns the delim table of the object if it exists else it 
#' computes the delim table from the TDR data.
#' with time and depth variables named "time" and "depth".
#' @param method the method to used. See help for functions in see also section 
#' to get details about possible options. When \code{"vspd"} or \code{"halsey"} 
#' methods fails \code{"std"} is used instead.
#' @param dvs a delim table such as return by \code{\link{dive_delim}}. If NULL 
#' \code{\link{dive_delim}} will be used with default arguments to compute one.
#' @param ... arguments to be passed to \code{\link{bottom_delim_vspd}}, 
#' \code{\link{bottom_delim_halsey}}, \code{\link{bottom_delim_bsm}} or 
#' \code{\link{bottom_delim_std}} according to the specified method.
#' @export
#' @examples
#' data(exses)
#' dvs_vspd   <- bottom_delim(exses$tdr)
#' dvs_halsey <- bottom_delim(exses$tdr, method = "halsey")
#' dvs_bsm    <- bottom_delim(exses$tdr, method = "bsm")
#' dvs_std    <- bottom_delim(exses$tdr, method = "std")
#' 
#' # Example using the four methods
#' n_dv <- 100
#' opar <- par(no.readonly = TRUE) ; par(mfrow = c(2, 2))
#' exses$delim <- dvs_vspd
#' tdrply(plot, 1:2, no = n_dv, obj = exses, main = 'method = "vspd"')
#' tdrply(points, 1:2, ty = "_", no = n_dv, obj = exses, col = "red")
#' exses$delim <- dvs_halsey
#' tdrply(plot, 1:2, no = n_dv, obj = exses, main = 'method = "halsey"')
#' tdrply(points, 1:2, ty = "_", no = n_dv, obj = exses, col = "red")
#' exses$delim <- dvs_bsm
#' tdrply(plot, 1:2, no = n_dv, obj = exses, main = 'method = "bsm"')
#' tdrply(points, 1:2, ty = "_", no = n_dv, obj = exses, col = "red")
#' exses$delim <- dvs_std
#' tdrply(plot, 1:2, no = n_dv, obj = exses, main = 'method = "std"')
#' tdrply(points, 1:2, ty = "_", no = n_dv, obj = exses, col = "red")
#' par(opar)
bottom_delim <- function(obj, method = c("vspd", "haslsey", "bsm", "std"), dvs = NULL, ...) {
  UseMethod("bottom_delim")
}

#' @inheritParams bottom_delim
#' @export
#' @keywords internal
bottom_delim.default <- function(obj, method = c("vspd", "halsey", "bsm", "std"), 
                                 dvs = NULL, ...) {
  method <- match.arg(method, method) 
  delim_fun <- switch(method, bottom_delim_std, "bsm" = bottom_delim_bsm, 
                       "vspd" = bottom_delim_vspd, "halsey" = bottom_delim_halsey)
  
  # First delimitate dives and surfaces
  if (is.null(dvs)) dvs <- dive_delim(obj)
  dv_attr_nms <- c("thres.depth", "thres.dur", "ignored_dives")
  dv_attr <- attributes(dvs)[dv_attr_nms]
  obj[ , "time"] <- as.numeric(obj[ , "time"])
  obj <- as.ses(tdr = obj, delim = dvs)
  bttNms <- c("btt_st_idx", "btt_ed_idx")
  isdv <- (dvs$no_dive > 0)
  
  if (method == "halsey") {
    tmp <- dvs[isdv, ] 
    offset <- tmp$st_idx - 1
    tmp <- tdrply(delim_fun, c("time", "depth"), ty = tmp, obj = obj, ...)
    tmp <- as.data.frame(rbindlist(tmp))
    dvs[isdv, bttNms] <- sweep(tmp[ , 1:2], 1, offset, "+")
    cnd <- ifelse(isdv, !tmp$success, FALSE)
    dvs$warning[cnd] <- upd_warn(dvs$warning[cnd], "used std method for bottom delim")
  }
  
  if (method == "bsm") {
    tmp <- dvs[isdv, ] 
    offset <- tmp$st_idx - 1
    tmp <- tdrply(delim_fun, c("time", "depth"), ty = tmp, obj = obj, ...)
    tmp <- as.data.frame(rbindlist(tmp))
    dvs[isdv, bttNms] <- sweep(tmp[ , 1:2], 1, offset, "+")
    cnd <- ifelse(isdv, !tmp$success, FALSE)
  }
  
  if (method == "vspd") {
    # Bottom delim 
    tmp <- dvs[isdv, ] 
    offset <- tmp$st_idx - 1
    tmp <- tdrply(delim_fun, c("time", "depth"), ty = tmp, obj = obj, ...)
    dvs[isdv, bttNms] <- sweep(as.data.frame(rbindlist(tmp)), 1, offset, "+")
    
    # Check if delim return NAs
    chkok_vspd <- ifelse(isdv, !is.na(dvs$btt_st_idx), TRUE)
    dvs$warning[!chkok_vspd] <- upd_warn(dvs$warning[!chkok_vspd], "vspd delim failed.")
    # Check if bottom limits deeper than 40% of max depth
    ledge_depth <- abs(tdrply(max, "depth", ty = dvs[ , 1:2], obj = obj) * 0.40)
    st_depth <- abs(obj$tdr[dvs$btt_st_idx, "depth"])
    ed_depth <- abs(obj$tdr[dvs$btt_ed_idx, "depth"])
    chkok_depth <-  st_depth >= ledge_depth & ed_depth >= ledge_depth
    chkok_depth <- ifelse(isdv, chkok_depth, TRUE)
    chkok_depth <- ifelse(is.na(chkok_depth), FALSE, chkok_depth)
    cnd <- !chkok_depth & chkok_vspd
    dvs$warning[cnd] <- upd_warn(dvs$warning[cnd], "Shallow bottom limits.")
  } else {
    chkok_depth <- dvs$no_dive <= 0
  }
  
  if (method %in% c("vspd", "std")) {
    # Bottom delim
    tmp <- dvs[!chkok_depth, ]
    nodv_std <- dvs$no_dive[which(!chkok_depth)]
    offset <- tmp$st_idx - 1
    if (method == "vspd") tmp <- tdrply(bottom_delim_std, "depth", ty = tmp, obj = obj)
    else tmp <- tdrply(delim_fun, "depth", ty = tmp, obj = obj, ...)
    dvs[!chkok_depth, bttNms[1:2]] <- sweep(as.data.frame(rbindlist(tmp)), 1, offset, "+")
  }
  
  # Format & Set attributes
  out_nms <- c("st_idx", "ed_idx", "no_dive", "btt_st_idx", "btt_ed_idx", "warning")
  out <- dvs[ , out_nms]
  for (att in dv_attr_nms) {attr(out, att) <- dv_attr[[att]]}
  attr(out, "method") <- method
  attr(out, "method.args") <- list(...) %else% formals(delim_fun)
  out
}

#' @inheritParams bottom_delim
#' @export
#' @keywords internal
bottom_delim.ses <- function(obj, method = c("vspd", "haslsey", "bsm", "std"), dvs = NULL, ...) {
  if ("delim" %in% names(obj))
    return(obj$delim)
  bottom_delim.default(obj, method, dvs, ...)
}

#' Delimitate bottom phase of dive using brokenstick models
#' 
#' The bottom is defined as the dive period between the first brokenstick 
#' segment (descent) and the last brokenstick segment (ascent).
#' 
#' @param time time readings, sorted in chronological order.
#' @param depth depth readings, sorted in chronological order.
#' @param pts The number of points to use in the brokenstick model. The default 
#' (\code{npts = 6}) is used by CTD-SRDL tags (Sea Mammal Research unit, St Andrews 
#' University) to summarize dive profiles.
#' @keywords internal
#' @export
#' @examples 
#' data(exses)
#' ind(exses)
#' 
#' n <- 65
#' idx <- tdrply(bottom_delim_bsm, 1:2, "!_/", no = n)[[1]]
#' tdrply(plot, 1:2, "!_/", no = n, main = n)
#' tdrply(function(x, st, ed, ...) points(x[st:ed, ]), 1:2, no = n, la = idx)
bottom_delim_bsm <- function(time, depth = NULL, npts = 6) {
  bsm <- try(brokenstick(time, depth, npts), TRUE)
  if (is.error(bsm)) list(st = NA, ed = NA, success = FALSE)
  else list(st = bsm$pts[2], ed = bsm$pts[npts - 1], success = TRUE)
}

#' Delimitate bottom phase of dive using wiggles and steps
#' 
#' Improves \code{\link{bottom_delim_std}} by defining a lower ledge (which implies 
#' longer bottoms) but reducing the bottom phase to a more relevant period 
#' (where specific behavior seems to indicate intensive foraging search). This 
#' is achieved by selecting the period between the first and last step/wiggle 
#' (see function \code{\link{wiggles}}) deeper than the ledge.
#' 
#' @param time time readings, sorted in chronological order.
#' @param depth depth readings, sorted in chronological order.
#' @param ledge depth threshold, specified as a percentage which is compared 
#' to the maximum depth. 0.75 has been used for king pengins. For 
#' elephant seals: about 95\% of prey catch attemps occur at depth greater 
#' than 50\% of the maximum depth (default value 0.50).
#' @param vert_vel A threshold used to define the steps. 0.35 m/s has been used for 
#' king pengins.
#' @references Halsey, L.G., Bost, C.-A. & Handrich, Y. (2007) A thorough and 
#' quantified method for classifying seabird diving behaviour. 
#' Polar Biology, 30, 991–1004.
#' @keywords internal
#' @seealso \code{\link{wiggles}}
#' @export
#' @examples 
#' data(exses)
#' ind(exses)
#' 
#' n <- 65
#' idx <- tdrply(bottom_delim_halsey, 1:2, "!_/", no = n)[[1]]
#' tdrply(plot, 1:2, "!_/", no = n, main = n)
#' tdrply(function(x, st, ed) points(x[st:ed, ]), 1:2, no = n, la = idx)
bottom_delim_halsey <- function(time, depth = NULL, ledge = 0.50, vert_vel = 0.35) {
  xy <- as.data.frame(xy.coords(time, depth)[1:2])
  rks0 <- bottom_delim_std(xy$y, ledge)
  df  <- xy[do.call(seq, unname(rks0)), ]
  tbl <- wiggles(df, step = vert_vel, output = "table")
  tbl <- tbl[tbl$type %in% c("step", "wiggle"), ]
  if (nrow(tbl) == 0) return(c(bottom_delim_std(xy$y, 0.80), list(success = FALSE)))
  list(st = which(time == min(tbl$start_x)), ed = which(time == max(tbl$end)), 
       success = TRUE)
}

#' Delimitate bottom phase of dive using vertical speed threshold
#' 
#' Delimitate bottom which is defined as the period between the first and last 
#' moment in the dive were vertical speed (simplified signal using \code{\link{poly}}) 
#' exceeds a threshold.
#' 
#' @param time time readings, sorted in chronological order.
#' @param depth depth readings, sorted in chronological order.
#' @param vspd smoothed verical speed, sorted in chronological order. Optional. 
#' If provided, function is faster as it does not need to compute it from 
#' depth and time.
#' @param vert_vel the verical speed threshold.
#' @param w a window width to use for vertical speed smoothing (moving average) 
#' when \code{vspd} is not provided.
#' @keywords internal
#' @export
#' @examples
#' data(exses)
#' ind(exses)
#' 
#' n <- 65
#' idx <- tdrply(bottom_delim_vspd, 1:2, "!_/", no = n)[[1]]
#' tdrply(plot, 1:2, "!_/", no = n, main = n)
#' tdrply(function(x, st, ed) points(x[st:ed, ]), 1:2, no = n, la = idx)
bottom_delim_vspd <- function(time = NULL, y = NULL, y.type = c("depth", "vspd"), 
                              vert_vel = 0.75, w = 21) {
  y.type <- match.arg(y.type, y.type)
  if (is.null(y) && length(time) == 1) 
    stop("time-depth or vertical speed must be provided")
  xy <- as.data.frame(xy.coords(time, y)[1:2])
  if (y.type == "depth") {
    xy$y <- rollapply(c(NA, diff(xy$y))/c(NA, diff(xy$x)), mean, w, na.rm = TRUE)
  }
  mod <- lm(y ~ stats::poly(x, degree = 4, raw = TRUE), xy)
  btt <- per(abs(predict(mod)) < abs(vert_vel))
  if (sum(btt$value) == 0) c(NA, NA)
  else setNames(as.list(btt[btt$value, 1:2]), c("st", "ed"))
}

#' Delimitate bottom phase of dive using depth threshold
#' 
#' Classic bottom delimitation method that defines bottom as the period between 
#' the first and last passage at a depth specified as a percentage of the 
#' maximum depth in the dive.
#' 
#' @param depth depth readings, sorted in chronological order.
#' @param ledge depth threshold, specified as a percentage which is compared 
#' to the maximum depth.
#' @keywords internal
#' @export
#' @examples 
#' data(exses)
#' ind(exses)
#' 
#' n <- 65
#' idx <- tdrply(bottom_delim_std, 2, "!_/", no = n)[[1]]
#' tdrply(plot, 1:2, "!_/", no = n, main = n)
#' tdrply(function(x, st, ed) points(x[st:ed, ]), 1:2, no = n, la = idx)
bottom_delim_std <- function(depth, ledge = 0.80) {
  x <- abs(depth)
  ledge_depth <- max(x, na.rm = TRUE) * ledge
  rk <- which(x >= ledge_depth)
  if (length(rk) == 0) rk <- which.max(x) + c(-1, 1)
  setNames(as.list(range(rk)), c("st", "ed"))
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

