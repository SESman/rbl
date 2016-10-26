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
is_dive_truncated <- function(x, max_diff = 5, error = FALSE) {
  out <- abs(first(x) - last(x)) > max_diff
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

#' Fix duplicates in TDR time stamps
#' 
#' @param obj The TDR dataset
#' @param time_seq A variable to use in order to subset \code{obj} and extract 
#' the time column. The time variable has to be in \code{POSIXct} format.
#' @param verbose A logical indicating if the function should be verbose.
#' @author Simon Wotherspoon, Yves Le Bras
#' @details Table with time stamp and number of replicate of duplicates is 
#' returned in attributes.
#' @export
#' @keywords internal correction
#' @examples
#' data(exses)
#' tmp <- exses$tdr[c(1,1,1,2:30), ]
#' tmp <- correct_duplicated_time(tmp)
#' attr(tmp, "correct_time")
correct_duplicated_time <- function(obj, time_seq = 1, verbose = FALSE) {
  # Check time specification
  if (!(length(time_seq) == 1) && is.POSIXct(obj[ , time_seq])) {
    stop("time_seq must specify a column of POSIXct times")
  } else {
    tzone <- attr(obj[ , time_seq], "tzone")
  }
  # Round times to nearest second
  obj[ , time_seq] <- .POSIXct(round(as.numeric(obj[ , time_seq])), tz = tzone)
  
  # Fix time duplicates
  dup <- duplicated(obj[ , time_seq])
  if (any(dup)) {
    # Extract indices of any duplicated row
    rs <- which(obj[ , time_seq] %in% obj[dup, time_seq])
    obj.dup <- obj[rs, ]
    # Aggregate duplicates according to data type of non time columns
    for (jj in setdiff(seq_along(obj.dup), time_seq)) {
      if (is.numeric(obj.dup[, jj])) # mean
        obj.dup[, jj] <- ave(obj.dup[ , jj], obj.dup[, time_seq])
      else if (is.logical(obj.dup[ , jj])) # or
        obj.dup[, jj] <- ave(obj.dup[ , jj], obj.dup[ , time_seq], FUN = any)
      else { # count for factor/character
        warning('"obj" contains non numeric/logical data. 
                Time duplicates aggregated with table()')
        obj.dup[, jj] <- ave(obj.dup[ , jj], obj.dup[ , time_seq], FUN = list %c% table)
      }
    }
    # Replace in original then remove duplicated
    obj[rs, ] <- obj.dup
    obj <- obj[!dup, ]
    # Get duplicates info
    msg_tab <- table(obj.dup[ , time_seq])
    if (verbose) {
      message('Duplicated time stamp(s) were detected. These rows are averaged.
              Details:')
      print(msg_tab)
    }
  } else {
    msg_tab <- "No duplicate found"
  }
  # Save info in output attributes
  attr(obj, "correct_time")$duplicates <- msg_tab
  
  invisible(obj)
}

#' Fix jumps and duplicates in TDR time stamps
#' 
#' @param obj The TDR dataset
#' @param time_seq A variable to use in order to subset \code{obj} and extract 
#' the time column. The time variable has to be in \code{POSIXct} format.
#' @param verbose A logical indicating if the function should be verbose.
#' @export
#' @author Simon Wotherspoon, Yves Le Bras
#' @details Output rows will be chronologically ordered. Table with time stamp 
#' and number of replicate of duplicates is returned in attributes, as well as 
#' a data frame listing "time jumps". 
#' @keywords correction
#' @examples
#' data(exses)
#' tmp <- rbind(exses$tdr[c(1,1,2), ], exses$tdr[-c(10, 20:30), ])
#' tmp <- correct_time(tmp)
#' attr(tmp, "correct_time")
correct_time <- function(obj, time_seq = 1, verbose = FALSE) {
  # Check time specification
  if (!(length(time_seq) == 1) && is.POSIXct(obj[ , time_seq])) {
    stop("time_seq must specify a column of POSIXct times")
  } else {
    tzone <- attr(obj[ , time_seq], "tzone")
  }
  # Round times to nearest second
  obj[ , time_seq] <- .POSIXct(round(as.numeric(obj[ , time_seq])), tz = tzone)
  
  # Fix duplicated times 
  obj <- correct_duplicated_time(obj, time_seq,verbose)
  
  # Get sampling interval
  dt <- time_reso(obj[, time_seq])
  
  # Generate complete and ordered sequence of times
  full_tms <- seq(min(obj[ , time_seq]), max(obj[ , time_seq]), dt)
  # Expand data and fill in missing times.
  idx <- match(full_tms, obj[ , time_seq])
  # Get "jumps" locations and save info in output attributes
  msg_tab <- subset(per(is.na(idx)), value, select = 1:2)
  msg_tab <- data.frame(from = full_tms[msg_tab$st_idx - 1], 
                        to = full_tms[msg_tab$ed_idx + 1])
  msg_tab$time_diff <- msg_tab[,2] - msg_tab[,1] - 1
  attr(obj, "correct_time")$jumps <- msg_tab
  if (verbose) {
    message(nrow(msg_tab), ' temporal jump(s) detected: ')
    if (nrow(msg_tab) != 0) print(msg_tab)
  }
  # Fix jumps & sort chronologically
  obj <- obj[idx, ]
  obj[ , time_seq] <- full_tms
  # Fix row names
  row.names(obj) <- seq_along(full_tms)
  
  invisible(obj)
}
