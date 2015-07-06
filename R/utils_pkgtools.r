#' pmean
#' 
#' Returns the (parallel) average of the input values.
#' @param ... numeric or logical arguments
#' @param na.rm	a logical indicating whether missing values should be removed.
#' @keywords internal
#' @export
#' @examples
#' pmean(1:3, 3:1)
#' pmean(1:3, 0)
pmean <- function (..., na.rm = FALSE) {
  tmp <- Map(c, ...)
  sapply(tmp, mean, na.rm = na.rm)
}

#' Update warning column in delim table
#' 
#' @param x Subset of warn column
#' @param  msg Message to append
#' @keywords internal
#' @export
upd_warn <- function(x, msg) {
  paste0(ifelse(is.na(x), "", x), ifelse(is.na(x), "", "; "), msg)
}

#' Angle average
#' 
#' @param x angle in radians.
#' @keywords internal angle
#' @export
#' @examples 
#' agl_mean(c(-pi, pi))
agl_mean <- function(x) {
  sinr <- sum(sin(x), na.rm = TRUE)
  cosr <- sum(cos(x), na.rm = TRUE)
  atan2(sinr, cosr)
}

#' Rescale angle to [-pi; pi]
#' 
#' @param x angle in radians.
#' @keywords internal angle
#' @export
#' @examples 
#' agl_rescale(5*pi / c(-4, 4))
agl_rescale <- function(x) {
  atan2(sin(x), cos(x))
}

#' from time stamp to row number
#' 
#' @param x a POSIXct vector
#' @export
#' @keywords internal
#' @details assumes that TDR data is sorted by time.
#' @examples 
#' data(exses)
#' ind(exses)
#' x <- sample(1:nrow(exses$tdr), 100)
#' identical(which.row(exses$tdr$time[x]), x)
which.row <- function(x, obj = ind()) {
  obj$tdr <- data.table(obj$tdr, key = "time")
  out <- obj$tdr[ , .I[time %in% x]]
  out[order(order(x))]
}

#' Match values against a data.frame with start and end values 
#' 
#' @param x the values to be matched against \code{ref}
#' @param ref a data.frame with start values in the first column end values in 
#' the second column and an optional id number in the third column.
#' @return for each \code{x} value, the row number of \code{ref} where \code{x} 
#' lies between start and end values. If \code{ref} has a third column (an id) 
#' its value is returned instead of the row number. When x value matches a start 
#' and a end value the priority is given to the start.
#' @keywords internal
#' @export
which.bw <- function(x, ref) {
  first_ed_greater  <- sapply(x, function(x) {
    tmp <- which(x < ref[ , 2])
    vals <- ref[tmp, 2]
    if (length(tmp) == 0) 0 else tmp[which.min(vals)]
  })
  last_st_less_eq <- sapply(x, function(x) {
    tmp <- which(x >= ref[ , 1])
    vals <- ref[tmp, 1]
    if (length(tmp) == 0) NA else tmp[which.max(vals)]
  })
  first_ed_eq  <- sapply(x, function(x) {
    tmp <- which(x == ref[ , 2])
    vals <- ref[tmp, 2]
    if (length(tmp) == 0) NA else tmp[which.min(vals)]
  })
  rks <- ifelse(first_ed_greater == last_st_less_eq, last_st_less_eq, NA)
  rks <- ifelse(is.na(rks), first_ed_eq, rks)
  if (ncol(ref) == 3) ref[rks, 3] else rks
}

#' Find to which specific dive/surface a instant belongs to
#' 
#' @param x The time (format \code{POSIXct}) or a integer giving the row number.
#' @param object A \code{ses} object such as returned by \code{\link{as.ses}}.
#' @export
which.dive <- function(x, object = ind()) {
  if (is.POSIXct(x)) {
    ref <- data.frame(
      st = object$tdr[object$delim[ , 1], 1], 
      ed = object$tdr[object$delim[ , 2], 1], 
      id = object$delim[ , 3])
  } else {
    ref <- object$delim[ , 1:3]
  }
  which.bw(x, ref)
}

#' x with(in/out) y
#' 
#' @param x Vector or NULL: the values to be matched.
#' @param y Vector or NULL: the values to be matched against. 
#' @export
#' @keywords internal
#' @examples
#' (1:10) %w/i% c(3,7,12) # 3 7
'%w/i%' <- function(x, y) x[x %in% y]

#' @rdname grapes-w-slash-i-grapes
#' @inheritParams grapes-w-slash-i-grapes
#' @export
#' @keywords internal
#' @examples
#' (1:10) %w/o% c(3,7,12) # 1  2  4  5  6  8  9 10
'%w/o%' <- function(x, y) x[!x %in% y]

#' Scale a series between two values
#' 
#' \code{rescale} is a utility to resize the range of values while keeping 
#' the original spacing between values.
#' 
#' @param x Numeric vector.
#' @param to Output range.
#' @param from Input range to be rescaled to \code{to}. Default is the range of \code{x}.
#' @keywords internal
#' @export
#' @examples
#' x <- -10:10
#' rescale(x)
#' rescale(x, to = c(-1, 3))
#' rescale(x, from = c(5, max(x)), to = c(0, 10))
rescale <- function (x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  if (length(to)   > 2) to   <- range(to)
  if (length(from) > 2) from <- range(from)
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

#' Extract numbers in character strings
#' 
#' @param x Atomic vector or list.
#' @param simplify Logical or character string. Should the result be simplified 
#' to a vector, matrix or higher dimensional array if possible? 
#' @keywords internal
#' @export
#' @examples
#' # Atomic character
#' x <- levels(cut(1:100, 3))
#' (out <- numIn(x))
#' 
#' # Atomic factor is coerced to character
#' x <- unique(cut(1:100, 3))
#' identical(numIn(x), out) 	# TRUE
#' 
#' # Works on list as well
#' x <- do.call(list, as.list(x))
#' identical(numIn(x), out) 	# TRUE
#' 
#' # When type is not character or factor the names are used
#' x <- do.call(list, as.list(1:3))
#' names(x) <- unique(cut(1:100, 3))
#' identical(numIn(x), out) 	# TRUE
#' 
#' # If names is NULL or empty the row.names are used instead
#' x <- matrix(1:6, 3)
#' row.names(x) <- levels(cut(1:100, 3))
#' is.null(names(x)) 		# TRUE
#' identical(numIn(x), out) 	# TRUE
numIn <- function(x, simplify = FALSE) {
  if (is.recursive(x)) {
    if (any(sapply(x, function(x) !is.character(x)))) {
      x <- if (all(sapply(x, is.factor))) lapply(x, as.character)
      else names(x) %else% row.names(x)
    }
  } else {
    if (is.numeric(x)) x <- names(x) %else% row.names(x)
  }
  m <- gregexpr('-?[0-9]+\\.?([0-9]*e(\\+|-))?[0-9]*', x)
  mtch <- if (is.list(x)) mapply(function(x, m) regmatches(x, list(m)), x, m)
  else regmatches(x, m)
  sapply(mtch, as.numeric, simplify = simplify)
}

#' Special operator to test if numeric values belong to a given range
#' 
#' %bw% for "between". Values are evaluated against the upper and lower 
#' bounds with \code{<=} and \code{>=} operators.
#' 
#' @param x numeric values
#' @param int range. Can have more than two elements. Atomic vectors are interpreted 
#' as a single condition while lists as a list of conditions (recycled if needed).
#' @export
#' @keywords internal
#' @examples
#' 1:10 %bw% c(2, 9)
#' 1:10 %bw% 2:10
#' 1:10 %bw% list(1:4, 1:2)
'%bw%' <- function (x, int) {
  if (is.atomic(int)) int <- list(int)
  .f <- function(x, int) { 
    if (all(is.na(int)) || all(is.na(int))) NA
    else x >= min(int, na.rm = TRUE) & x <= max(int, na.rm = TRUE)
  }
  mapply(.f, x, int)
}

#' Replace values in an atomic vector.
#' @param x The atomic vector
#' @param na.0 The value to be replaced. Default is NaN.
#' @param na.1 The replacement. Default is NA.
#' @keywords internal
#' @export
#' @examples
#' x <- sample(c(1:3,NaN), 20, replace=TRUE)
#' x 
#' replaceMissing(x)
replaceMissing <- function(x, na.0 = NaN, na.1 = NA) {
  if (is.nan(na.0)) x[is.nan(x)] <- na.1
  else if (is.na(na.0)) x[is.na(x)] <- na.1
  else x[is.na(x)] <- na.1
  x
}

#' Count the number of NAs in a vector
#' 
#' Shortcut for \code{compose(sum, is.na, unlist)}
#' 
#' @param x a vector whose elements are to be tested.
#' @return Return the number of \code{NA} in \code{x}.
#' @details As any number different from 0 return a \code{TRUE} when coerced to
#' logical, this function can be used in \code{if} statements.
#' @export
#' @keywords internal
#' @examples
#' x <- c(rep(NA, 3), 1:3)
#' nNA(x)
#' if (nNA(x)) {TRUE} else {FALSE}
#' if (nNA(1:3)) {TRUE} else {FALSE}
nNA <- function(x) sum(is.na(unlist(x)))

#' Count the number of unique values in a vector
#' 
#' Shortcut for \code{compose(length, unique)}. Count the number of distinct 
#' values in an atomic vector.
#' 
#' @param x a vector whose unique elements are to be counted.
#' @export
#' @keywords internal
#' @examples
#' nUN(rep(1:5, 5:1)) # 5
nUN <- function(x) length(unique(x))

#' Else special operator
#' 
#' Discard first value if \code{FALSE}, \code{NULL}, empty or \code{"try-error"}
#' 
#' @param val Normal output.
#' @param def Default output when \code{val} is \code{FALSE}, \code{NULL} or empty.
#' @export
#' @keywords internal
#' @examples
#' "abc" %else% "Another value is returned"
#' NULL %else% "Another value is returned"
#' try(log("abc"), silent = TRUE) %else% "Another value is returned"
'%else%' <- function (val, def = NA){
  if (identical(val, FALSE) || is.null(val) || length(val) == 0 || is.error(val)) def else val
}

#' Depth of an R object
#' @param x The object to analyse.
#' @export
#' @details function \code{plotrix::maxDepth}
#' @keywords internal
list_depth <- function (x) {
  if (is.list(x)) {
    if (identical(x, list())) return(0)
    maxdepth <- 1
    for (lindex in 1:length(x)) {
      newdepth <- list_depth(x[[lindex]]) + 1
      if (newdepth > maxdepth) 
        maxdepth <- newdepth
    }
  }
  else maxdepth <- 0
  return(maxdepth)
}

#' nstr
#' 
#' Recursive extraction of names (such as \code{names(c(x, recursive=TRUE))}) but
#' stops when a subelement is atomic (avoid long long run when launched on 
#' large object such as a TDR dataset).
#' 
#' @param x The object to analyse.
#' @export
#' @keywords internal
#' @examples
#' x <- data.frame(X=1:10, Y=10:1)
#' names(c(x, recursive=TRUE))
#' nstr(x)
nstr <- function(x) {
  n <- list_depth(x)
  name.vec <- c()
  if (n == 1){
    return(names(x))
  } else if (n > 1){
    for (i in seq_along(x)){
      name.vec <- c(name.vec,
                    names(x), 
                    paste(names(x)[i], nstr(x[[i]]), sep='.'))
    }
  }
  return(unique(name.vec[!grepl('\\.$', name.vec)]))
}

#' Search recurssively to a data.frame
#' 
#' This function is a helper designed to be used in \code{tdrply}. It searches 
#' recurssively in an object for a list of data.frames at a given level of depth. 
#' If the search ends to atomic vectors, the function builds a list of 
#' data.frames by taking their elements by two successively.
#' 
#' @param .idx An object.
#' @export
#' @keywords internal
#' @seealso \code{\link{tdrply}}
#' @examples
#' .idx <- list(1:10, list(1:10))
#' # df_search(.idx) # error
#' .idx <- list(1:10, data.frame(1:10, 1:10))
#' # df_search(.idx) # error
#' 
#' .idx <- data.frame(1:10, 10:1)
#' df_search(.idx)
#' .idx <- list(a = .idx, b = .idx)
#' df_search(.idx)
#' .idx <- list(a = 1:3, b = 1:10)
#' df_search(.idx)
df_search <- function(.idx) {
  if (is(.idx, 'data.frame')) return(.idx)
  # Test the type of the elements
  tfuns <- list(df = function(x) is.data.frame(x), 
                lst = function(x) inherits(x, 'list'), 
                atm = function(x) is.atomic(x))
  tres <- lapply(tfuns, function(f) sapply(.idx, f))
  
  # Check if the results are even for each type
  tresHomo <- sapply(tres, function(x) Reduce(identical, x == x[1]))
  if (any(!tresHomo)) stop('.idx must have an evenly nested structure')
  
  # Get matching type and check it is unique
  type <- names(tres)[sapply(tres, unique)]
  if (length(type) != 1) stop('Unexpected type(s) found in .idx', str(.idx))
  
  # Apply function to elements
  .roll <- function(x) rollapply(x, c, 2, aty = 'm', nas = FALSE, simplify = TRUE)
  switch(type, df = .idx, lst = .roll(lapply(.idx, df_search)), 
         atm = lapply(.idx, function(x) as.data.frame(t(.roll(x)))))
}

#' Set and get the current individual
#' 
#' @param value If provided this value becomes te current individual. If omited 
#' the function return the last declared individual.
#' @param cache Should the object be copied in a cache rather than a link to 
#' the object ?
#' @export
#' @seealso \code{ind} is convenient to use with \code{\link{tdrply}}.
#' @examples
#' data(exses)
#' ind(exses)
#' exses$test <- "test!"
#' identical(ind(), exses)
ind <- function(value, cache = FALSE) {
  if (missing(value)) {
    if (cache == TRUE) {
      cache <- get("cache", envir = .GlobalEnv)
      return(cache$ind)
    } else {
      cache <- get("cache", envir = .GlobalEnv)
      return(eval(cache$link$val, cache$link$env))
    }
  } else {  
    if (!exists("cache", .GlobalEnv))
      assign("cache", list(), envir = .GlobalEnv) 
    if (cache == TRUE) {
      cache$ind <<- value
    } else {
      cache$link <<- list(env = parent.frame() , val = substitute(value))
    }
  }
  invisible(NULL)
}

#' is.error
#' 
#' @param x The objet to proceed
#' @export
#' @keywords internal
is.error <- function(x) inherits(x, "try-error")

#' floorPOSIXct
#' 
#' @param x The POSIXct vector
#' @param units How to cut the values. The units are partially matched in 
#' \code{c('secs', 'mins', 'hours', 'days')}. A number can precede the unit.
#' @param offset To use in the case where a cut occurs at a inconvenient
#' date (see examples).
#' @export
#' @keywords internal
#' @examples
#' data(exses)
#' x <- exses$stat$time - 304*(24*3600)
#' plot(x, x, type = 'l')
#' lines(x, floorPOSIXct(x, '2days'), col = 'blue', type = 's')
#' # To force the cut to occur on the 1st January
#' lines(x, floorPOSIXct(x, '2days', '1d'), col = 'lightblue', type = 's')
#' lines(x, floorPOSIXct(x, 'days'), col = 'red', type = 's')
#' lines(x, floorPOSIXct(x, '0.5d'), col = 'green', type = 's')
floorPOSIXct <- function(x, units = "days", offset = '0 days', ...) {
  if (is.numeric(units)) stop("'units' must be a character string.")
  opt <- list(units, offset)
  n <- sapply(opt, function(x) unlist(numIn(x)) %else% 1)
  u <- mapply(function(x, n) gsub(paste0(as.character(n), '|\\ '), '', x), opt, n)
  chc <- c('secs' = 1, 'mins' = 60, 'hours' = 3600, 'days' = 86400)
  o <- chc[pmatch(u, names(chc), NA, duplicates.ok = TRUE)] * n 
  if (nNA(o)) stop('Unknown unit found: ', paste(u, collapse=', '))
  as.POSIXct(floor((as.numeric(x) - o[2]) / o[1]) * o[1]  + o[2], 
             tz = attr(x, 'tzone'), origin = '1970-01-01')
}

#' Decompose an atomic vector to its successive values and their length.
#' 
#' The reverse of 'base::rep()' function: decompose an atomic vector to its successive 
#' values and their length.
#' 
#' @param x The atomic vector to examine.
#' @param idx Should the indexes (start and end) of homogeneous sequences be 
#' returned as well ?
#' @return A data frame with values and lengths of the homogeneous sequences 
#' of x. The class of the column 'value' is copied from the input.
#' @keywords internal
#' @export
#' @examples
#' (x <- rep(LETTERS[1:10], 10:1))
#' (y <- per(x))
#' identical(rep(y$value, y$length), x)   # TRUE
#' inherits(y$value, class(x))            # TRUE
per <- function(x, idx = TRUE) {
  x.org <- x
  if (is.logical(x) || is.factor(x)) {x <- as.numeric(x)}
  else if (is.character(x)) {x <- as.numeric(as.factor(x))}
  
  chg <- diff(x)
  end <- c(which(chg != 0), length(x))
  start <- c(1, end[-length(end)] + 1)
  
  out <- if (idx) 
    data.frame(st_idx = start, ed_idx = end, 
               value = x.org[start], length = end - start + 1, 
               stringsAsFactors = FALSE)
  else 
    data.frame(value = x.org[start], length = end - start + 1, 
               stringsAsFactors = FALSE)
  
  class(out) <- c("per", "data.frame")
  out
}

#' is.POSIXct
#' 
#' @param x The objet to proceed
#' @export
#' @keywords internal
is.POSIXct <- function (x) is(x, "POSIXct")

#' Linear interpolation
#' 
#' @param x A vector with missing values to interpolate.
#' @param n_max The maximun number of successive missing values to interpolate.
#' @export
#' @keywords internal
li <- function(x, n_max = NULL) {
  to_interpolate <- is.na(x) | is.nan(x)
  seqs <- per(to_interpolate)
  if (!is.null(n_max)) {
    seqs$value[seqs$value & seqs$length > n_max] <- FALSE
    seqs <- per(rep(seqs$value, seqs$length))
  }
  
  st_idx <- seqs$ed_idx[!seqs$value][-sum(!seqs$value)]
  ed_idx <- seqs$st_idx[!seqs$value][-1]
  n_vals  <- ed_idx - st_idx + 1
  
  li_out <- Map(seq, from = x[st_idx], to = x[ed_idx], length.out = n_vals)
  for (ii in seq_along(li_out)) {
    x[st_idx[ii]:ed_idx[ii]] <- li_out[[ii]]
  }
  x
}
