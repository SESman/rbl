#' Check if a dive profile is complete or truncated
#' 
#' @param x The depth sequence of a dive profile.
#' @param max_diff The maximum depth difference tolerated between the start and 
#' the end of the dive before the dive is considered as truncated.
#' @param error Should the function give an error if a truncated dive is 
#' encountered ?
#' @export
#' @examples 
#' \dontrun{
#' data(exses)
#' tdrply(is_dive_truncated, "depth", error = TRUE, obj = exses)
#' }
is_dive_truncated <- function(x, max.diff = 5, error = FALSE) {
  out <- abs(first(x) - last(x)) > max.diff
  if (error && any(out)) stop("Truncated dive detected")
  invisible(out)
}

#' Resize an individual
#' 
#' \code{resize} cuts a \code{ses} object: given a set of dives to discard, it 
#' remove the according rows from the TDR dataset and the statdive table. The dive
#' and bottom delimitation ("delim" table) is updated as well.
#' 
#' @param x An object of class \code{ses}.
#' @param no_dv_min The number of the first dive to keep.
#' @param no_dv_max The number of the last dive to keep.
#' @param rm Optional. A logical vector of length = \code{nrow(x$delim)}, where 
#' \code{TRUE} indicate that dives are to be removed.
#' @export
#' @keywords correction
#' @examples
#' data(exses)
#' range(exses$tdr$time)
#' range(abs(exses$delim$no_dive))
#' range(exses$stat$no_dive)
#' 
#' new.ses <- resize(exses, 50, 55)
#' range(new.ses$tdr$time)
#' range(abs(new.ses$delim$no_dive))
#' range(new.ses$stat$no_dive)
resize <- function(x, no_dv_min = 1, no_dv_max = max(x$delim$no_dive), rm = NULL) {
  dvs <- x$delim
  nr <- nrow(dvs)
  keep <- if (!is.null(rm)) {!rm} else {rep(TRUE, nr)}
  
  # Delete surface at start of the dataset
  if (dvs$no_dive[1] <= 0) keep[1] <- FALSE
  # Delete surface at end of the dataset
  if (dvs$no_dive[nr] <= 0) keep[nr] <- FALSE
  # Delete last dive if no surface after
  if (dvs$no_dive[nr] > 0) keep[nr] <- FALSE
  # Delete dives & surfaces > no_dv_max (but keep post no_dv_max surface)
  if (no_dv_max < max(dvs$no_dive)) {
    rk <- which(dvs$no_dive == no_dv_max)
    keep[(rk+2):nr] <- FALSE
  }
  # Delete dives & surfaces < no_dv_min (but keep surface before no_dv_min)
  if (no_dv_min > min(dvs$no_dive[dvs$no_dive <= 0])) {
    rk <- which(dvs$no_dive == no_dv_min)
    keep[1:(rk-2)] <- FALSE
  }
  
  # Indice correction due to deletions
  nokp <- per(!keep)
  nokp <- nokp[nokp$value, ]                   # Deleted periods
  nokp$st_idx <- dvs$st_idx[nokp$st_idx]       # their start ..
  nokp$ed_idx <- dvs$ed_idx[nokp$ed_idx]       # and end indice
  nokp$length <- nokp$ed_idx - nokp$st_idx + 1 # number of rows deleted
  
  coridx <- function(x) {
    # Test if x's values belong to a deleted period
    COND <- mapply(function(st, ed) x %bw% c(st, ed), nokp$st_idx, nokp$ed_idx)
    isDlt <- apply(COND, 1, any) 
    # Number of deleted rows before x's values
    COND <- sapply(nokp$ed_idx, function(ed) ed < x)
    nDlt <- apply(COND, 1, function(x) sum(nokp$length[which(x)])) 
    ifelse(isDlt, NA, x - nDlt)
  }

  dvs_att <- attributes(dvs)
  dvs_att$ignored_dives$st_idx <- coridx(dvs_att$ignored_dives$st_idx)
  
  # Subset tdr data
  dlt_rows <- unlist(Map(seq, nokp$st_idx, nokp$ed_idx))
  x$tdr <- x$tdr[-dlt_rows, ]
  
  # Subset and format dvs
  dvs[ , grep('idx', names(dvs))] <- lapply(dvs[ , grep('idx', names(dvs))], coridx)
  dvs <- dvs[keep, ]
  for (att in (names(dvs_att) %w/o% 'row.names')) {attr(dvs, att) <- dvs_att[[att]]}
  x$delim <- dvs
  
  # Subset statdives
  x$stat <- x$stat[x$stat$no_dive %in% dvs$no_dive, ]
  
  x
}

#' Correct depth sequence from drift of values over time
#' 
#' @param x The depth sequence.
#' @param dur A window width for rolling functions. Make sure to choose that 
#' argument so that the animal went back to the surface at least once in this 
#' window.
#' @param plt Should a plot about processing should be drawn ?
#' @param span Argument passed to \code{\link{loess}}.
#' @param ... Other arguments to be passed to \code{\link{loess}}.
#' @export
#' @keywords correction
#' @examples
#' # correct_depth(exses$tdr$depth, plt = TRUE)
correct_depth <- function (x, dur = 5000, plt = FALSE, span = 0.15, ...) 
{
  wds <- seq(1, length(x), by = dur)
  rks <- rollapply(wds, function(st, ed) (which.min(x[st:ed]) + st - 1) %else% NA, 
                   2, aty = "m")
  vals <- rollapply(wds, function(st, ed) min(x[st:ed], na.rm = TRUE), 
                    2, aty = "m")
  gam <- loess(vals ~ rks, span = span, family = "symmetric", ...)
  if (plt) {
    plot(vals ~ rks)
    lines(na.omit(rks), fitted(gam), col = "red", lwd = 2)
  }
  pred <- predict(gam, seq_along(x))
  out <- x - ifelse(is.na(pred), 0, pred)
  ifelse(out < 0, 0, out)
}

#' Fix jumps and duplicates in TDR time stamps
#' 
#' @param obj The TDR dataset
#' @param time_seq A variable to use in order to subset \code{obj} and extract 
#' the time column. The time variable has to be in \code{POSIXct} format.
#' @param verbose A logical indicating if the function should be verbose.
#' @export
#' @details Can be use for sampling frequency <= 1 Hz.
#' @keywords correction
#' @import data.table
#' @examples
#' data(exses)
#' exses$tdr <- correct_time(exses$tdr, verbose = TRUE)
correct_time <- function(obj, time_seq = 1, verbose = FALSE) {
  x <- obj[ , time_seq]
  is.POSIXct(x) || stop("'x' ust be of class POSIXct.")
  seqs <- per(diff(round(x)))
  
  # Duplicates
  if (any(seqs$value == 0)) {
    nms_tot <- names(obj)
    tm_nm <- names(obj[ , time_seq, drop = FALSE])
    nms <- nms_tot %w/o% tm_nm
    
    dup <- seqs[seqs$value == 0, ]
    dup_obj <- Map(function(st, ed) obj[st:ed, ], dup$st_idx, dup$ed_idx + 1)
    
    fmla <- as.formula(paste0('cbind(', paste(nms, collapse = ','), ') ~ ', tm_nm))
    dup_obj <- lapply(dup_obj, function(x) aggregate(fmla, x, mean)[ , nms_tot])
    for (ii in seq_along(dup_obj)) {
      obj[dup$st_idx[ii], ] <- dup_obj[[ii]]
    }
    obj <- obj[-(dup$ed_idx[ii] + 1), ]
    x <- obj[ , time_seq]
    seqs <- per(diff(x))
    if (verbose)
      message('Duplicated time stamp(s) were detected. The implied rows are averaged.')
  }
  
  # Find sampling frequency
  tmp <- aggregate(length ~ value, seqs, sum)
  reso <- tmp$value[which.max(tmp$length)]
  
  # Find jumps
  cond <- seqs$value == reso
  ok <- seqs[cond, ]
  ano <- seqs[!cond, ]
  message('Sampling frequency: ', 1 / as.numeric(reso), ' Hz. ')
  if (verbose) {
    msg_tab <- setNames(ano[ , 1:3], c('from_row', 'to_row', 'time_difference'))
    msg_tab[ , 2] <- msg_tab[ , 2] + 1
    message(nrow(msg_tab), ' temporal jump(s) detected: ')
    if (nrow(msg_tab) != 0) print(msg_tab)
  }
  
  if (nrow(ano) == 0) {
    if (verbose)
      message('No missing time stamp, object returned as is.')
    return(invisible(obj))
  } else {
    miss_area <- Map(seq, from = x[ano$st_idx], 
                     to = x[ano$st_idx + 1] - reso, 
                     length.out = floor(as.numeric(ano$value) / as.numeric(reso)))
    
    warning('Missing time stamp(s). ', sum(sapply(miss_area, length) - 1), ' rows added.')
    miss_area <- lapply(miss_area, function(miss_area) {
      out <- setNames(as.data.frame(matrix(nrow = length(miss_area) - 1, ncol = ncol(obj))), 
                      names(obj))
      out[ , time_seq] <- miss_area[-1]
      out
    })
    ok_area <- Map(function(st, ed) obj[st:ed, ], ok$st_idx, ok$ed_idx + 1)
    
    min_time <- function(x) min(x[ , time_seq])
    seqs_order <- order(c(sapply(ok_area, min_time), sapply(miss_area, min_time)))
    tot <- c(ok_area, miss_area)[seqs_order]
    out <- data.table::rbindlist(tot, use.names = FALSE)
    out <- as.data.frame(out)
    row.names(out) <- seq(1, nrow(out))
    return(invisible(out))
  }
}
