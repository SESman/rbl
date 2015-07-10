#' Subset dive numbers with values as using indices
#' 
#' @param x a numeric or logical.
#' @param obj a "ses" object.
#' @return vector of selected dive numbers.
#' @export
#' @keywords internal tdr_interpreter
#' @examples
#' data(exses)
#' ind(exses)
#' no_interpreter(NULL)    # all dive numbers returned
#' 
#' # Define default values using "no_dive" slot of the ses object
#' exses$no_dive <- 1:15
#' no_interpreter(NULL)
#' 
#' no_interpreter(1:10)    # returns dive numbers 1:10
#' no_interpreter(-(1:10)) # all dive number returned but 1:10
no_interpreter <- function(x, obj = ind()) {
  no <- obj$no_dive
  if (is.null(no)) no <- obj$delim$no_dive[obj$delim$no_dive > 0]
  if (is.null(x)) {
    return(no)
  } else if (is.logical(x)) {
    any(x) || stop('At least one element of "no" should be TRUE.')
    return(no[x])
  } else {
    if (is.numeric(x)) {
      if (any(x > 0) && any(x < 0)) 
        stop("only 0's may be mixed with negative subscripts.")
      if (all(x < 0)) return(no %w/o% -x)
      else return(x)
    } else {
      stop('"x" must be logical or numeric')
    }
  }
}

#' Get detailed parse information from "ty" regular expressions
#' 
#' @param x a kind of regular expression to subet rows of TDR datasets.
#' @param all should both leaf and branches be returned.
#' @return a data frame with parsed text data.
#' @details Take advantage of \code{\link{getParseData}} to decompose 
#' "ty" expression.
#' @export
#' @keywords internal tdr_interpreter
#' @examples 
#' View(ty_interpreter('-|(!(_/))|#########(_)'))
ty_interpreter <- function(x, all = TRUE) {
  # Remove invalid characters
  x <- gsub("[^\\(|\\)|&|\\||!|/|~|_|-]", "", x)
  
  # Insert operators
  # "CONCATENATION" operator = "*": "CONCAT" priority > "OR" priority
  # "OR" replaced by "+": "OR" priority > "AND" priority
  while (grepl("[-!_/~\\)][-!_/~\\(]", x)) {
    x <- sub("([-!_/~\\)])([-!_/~\\(])", "\\1*\\2", x)
  }
  x <- gsub("\\|", "+", x)
  
  # Translate symbols to strings
  x <- gsub("-", "psf", x)
  x <- gsub("!", "dsc", x)
  x <- gsub("_", "btt", x)
  x <- gsub("/", "asc", x)
  x <- gsub("~", "nsf", x)
  
  # Parse expression
  df <- getParseData(parse(text = x), includeText = TRUE)
  
  # Translate strings to symbols
  df$text <- gsub("psf", "-", df$text)
  df$text <- gsub("dsc", "!", df$text)
  df$text <- gsub("btt", "_", df$text)
  df$text <- gsub("asc", "/", df$text)
  df$text <- gsub("nsf", "~", df$text)
  
  if (!all) df[df$terminal, ] else df
}

#' Extract indices corresponding to a symbol from a "delim" table
#' 
#' @param elt A single character to match against delim table.
#' @param no dives numbers to extract.
#' @param delim a delim table (delim slot of a "ses" object).
#' @return a data frame with row number of sart/end of the desired periods.
#' @export
#' @keywords internal tdr_extract
#' @examples
#' data(exses)
#' elt_delim('-', 60:62, exses$delim)
#' elt_delim("/", no = c(62, 62:60), delim = exses$delim)
elt_delim <- function(elt, no, delim) {
  is.character(elt) || stop('"elt" must be character')
  nchar(elt) == 1 || stop('"elt" must be ONE character')
  
  # Raw subset
  del_no <- with(delim, switch(elt, no_dive, "-" = 1 - no_dive, "~" = -no_dive))
  del_rw <- del_no %in% no
  del_cl <- list("-" = 1:2, "!" = c(1, 4), "_" = 4:5, "/" = c(5, 2), "~" = 1:2)
  out <- data.frame(delim[del_rw, c(del_cl[[elt]], 3)], no = del_no[del_rw])
  
  # Avoid redondancy between bottom and descent/ascent
  if (elt == "_"){
    out[ , 1] <-  as.numeric(out[ , 1]) + 1
    out[ , 2] <-  as.numeric(out[ , 2]) - 1
  }
  
  # Insert rows where there is no match for "no"
  no_abs <- no %w/o% out$no
  if (length(no_abs) != 0) {
    fac <- apply(outer(out$no, no_abs, ">"), 1, sum) * 2
    no_pres <- split(out, fac)
    no_abs <- lapply(no_abs, function(x) as.data.frame(as.list(c(rep(NA, 3), x))))
    subset_order <- order(c(unique(fac), seq_along(no_abs)))
    out <- setnames(rbindlist(c(no_pres, no_abs)[subset_order]), c('st', 'ed', 'no_dive', 'no'))
  }
  out <- as.data.frame(out)
  
  # Match info: "elt" and "elt_mtch"
  out$elt <- rep(elt, nrow(out))
  out$elt_mtch <- apply(out[ , 1:2], 1, function(x) !any(is.na(x)))
  
  out
}

#' Extract indices from "delim" table using a regexp-like syntax
#' 
#' @param ty a kind of regular expression to subet diving periods of a TDR dataset.
#' @param no the dive numbers to extract.
#' @param obj a "ses" object.
#' @param ... arguments to be passed to \code{\link{check_nomtch}}
#' @export
#' @import data.table
#' @keywords internal tdr_extract
#' @examples
#' data(exses)
#' ind(exses)
#' 
#' str(tmp <- ty_delim('_', no = 60))
#' str(tmp <- ty_delim('_&!', no = 60))
#' 
#' str(tmp <- ty_delim('_', no = 60:62))
#' str(tmp <- ty_delim('!_/', no = 60:62))
#' 
#' str(tmp <- ty_delim('-&!&_&/&~', no = 60:62))
#' str(tmp <- ty_delim('-!&/~', no = 60:62))
#' 
#' str(tmp <- ty_delim('-!|/~', no = 60:62))
#' str(tmp <- ty_delim('-(!|/)~', no = 60:62))
#' 
#' str(tmp <- ty_delim('-(!|/)&~', no = 60:62))
#' 
#' str(tmp <- ty_delim("-&(__)|(/~)&//", no = 60:62))
ty_delim <- function(ty = "!_/", no = NULL, obj = ind(), no_match = "ignore") {
  all(sapply(ty, is.character)) || stop('"ty" must be a character')
  if (length(ty) > 1) return(lapply(ty, tydelim, no = no, obj = obj, no_match = no_match))
  no <- no_interpreter(no, obj)
  ty_bak <- ty
  ty <- ty_interpreter(ty, all = TRUE)
  
  # Find command related to an id
  child <- function(id) {
    ty$id[which(ty$parent == id)]
  }
  brother <- function(id) {
    rk <- which(ty$id == id)
    child(ty$parent[rk]) %w/o% id
  }
  parent <- function(id) {
    ty$parent[which(ty$id == id)]
  }
  
  # Functions to check input data
  check_grp <- function(x) {
    # Make sure all data frames in have a "grp" column
    if (is.null(x)) return(x)
    if (!is.data.frame(x)) {
      if (length(x) > 1) return(lapply(x, check_grp))
      else return(check_grp(x[[1]]))
    } 
    if (is.null(x$grp)) data.frame(x, grp = unique(x$elt), stringsAsFactors = FALSE)
    else x
  }
  check_mtch <- function(x) {
    # Check grp match 
    if (is.null(x)) return(x)
    if (!is.data.frame(x)) return(lapply(x, check_mtch))
    if (!is.null(x$grp_mtch)) return(x)
    fac <- paste0(x$grp, "#", x$no)
    mtch <- tapply(x$elt_mtch, fac, mean)
    data.frame(x, grp_mtch = unsplit(mtch, fac))
  }
  check_no <- function(x) {
    # Check that is a list of element including a unique "no" identifier
    if (is.null(x)) return(x)
    if (!is.data.frame(x)) return(lapply(x, check_no))
    else "if"(nUN(x$no) == 1, x, split(x, paste0(x$grp, "#", x$no)))
  }
  check <- function(x) check_no(check_mtch(check_grp(x)))
  
  # Operators
  cat_op <- function(x) {
    # Merge two groups into a single
    x <- lapply(x, check)
    if (any(sapply(x, is.null))) return(NULL)
    while (list_depth(x) > 2) x <- unlist(x, recursive = FALSE)
    x <- as.data.frame(rbindlist(x))
    x$grp <- unsplit(tapply(x$elt, x$no, function(x) do.call(paste0, as.list(x))), x$no)
    split(x, paste0(x$grp, "#", x$no))
  }
  or_op <- function(...) {
    # Choose best match for each no between input groups
    x <- lapply(as.list(...), check)
    if (any(sapply(x, is.null))) return(NULL)
    nUN(sapply(x, length)) == 1 || stop("Inputs must have the same length. ", 
                                        'Mixed "|" and "&" operators ?')
    out <- x[[1]]
    for (rk in seq_along(out)) {
      xs <- lapply(x, function(x) x[[rk]])
      mtch <- sapply(xs, function(x) unique(x$grp_mtch))
      out[[rk]] <- xs[[which.max(mtch)]]
    }
    out
  }
  and_op <- function(x) {
    # Concatenate inputs in a list
    x <- lapply(x, check)
    if (any(sapply(x, is.null))) return(NULL)
    if (list_depth(x) == 2) return(x)
    unlist(x, recursive = FALSE)
  }
  
  # Order "ty" by token type (from basic to master)
  tokens <- c(xtract = "SYMBOL", or = "'+'", opar = "'('", cpar = "')'", 
              cat = "'*'", and = "AND", other = "expr")
  token_priority <- factor(ty$token, tokens)
  ty$parent[ty$parent == 0] <- Inf
  ty <- ty[order(token_priority, ty$parent), ]
  
  # Evaluate ty, command after command until master is evaluated
  master_id <- ty$id[which(is.infinite(ty$parent))]
  tmp <- replicate(max(ty$id), NULL, simplify = FALSE)
  while (is.null(tmp[[master_id]])) {
    for (id in ty$id) {
      if (!is.null(tmp[[id]])) next
      rk <- which(ty$id == id)
      tok <- ty$token[rk]
      txt <- ty$text[rk]
      
      if (tok == tokens["xtract"]) {
        tmp[[id]] <- elt_delim(txt, no, obj$delim)
      } else if (tok == tokens["cat"]) {
        tmp[[id]] <- cat_op
      } else if (tok == tokens["or"]) {
        tmp[[id]] <- or_op
      } else if (tok == tokens["and"]) {
        tmp[[id]] <- and_op
      } else if (tok %in% tokens[grep("par", names(tokens))]) {
        tmp[[id]] <- identity
      } else if (tok == tokens["other"]) {
        if (is.infinite(parent(id)) && any(sapply(tmp[child(id)], is.null))) next
        if (txt %in% c("-", "!", "_", "/", "~")) {
          tmp[[id]] <- check(tmp[[child(id)]])
        } else {
          chld <- tmp[child(id)]
          fun <- Filter(is.function, chld) %else% next
          tmp[[id]] <- fun[[1]](Filter(is.list, chld))
        }
      }
    }
  }
  out <- check_nomtch(tmp[[max(ty$id)]], no_match, obj)
  if (is.data.frame(out)) 
      out <- setNames(list(out), paste0(ty_bak, "#", no))
  if (is.null(names(out))) 
      out <- setNames(out, paste0(sapply(out, function(x) paste0(x$elt)), "#", no))
  out <- out[order(sapply(strsplit(names(out), "#"), function(x) as.numeric(x[[2]])))]
  class(out) <- c(class(out), "ty")
  out
}

#' Deal with missing periods
#' 
#' \code{check_nomatch} search for periods without entry in "delim" table and 
#' check if a warning/error is to be thrown.
#' 
#' @param x a data.frame with "grp" match information.
#' @param no_match a string indicating how to handle periods with no match.
#' @param obj a "ses" object.
#' @return updated x, warning /error.
#' @export
#' @keywords internal tdr_extract
check_nomtch <- function(x, no_match, obj = ind()) {
  if (!is.data.frame(x)) {
    return(lapply(x, check_nomtch, no_match = no_match))
  }
  ref_order <- c('-', '!', '_', '/', '~')
  # List periods with no match
  mtch <- x[!x$elt_mtch, c("elt", "no")]
  
  if (nrow(mtch) != 0) {
    msg <- paste('Some periods have no matching entry in "delim" table of "obj".', 
                 'Found for "no" in ', paste(unique(mtch$no), collapse = ", "))
    # If "no match" found check if is included in a dive cycle
    order_divecycle <- function(x) {
      tmp <- factor(x, levels = ref_order)
      order(tmp)
    }
    is.in_divecycle <- function(x) {
      tmp <- unique(x$elt)[order_divecycle(unique(x$elt))]
      grepl(do.call(paste0, as.list(tmp)), do.call(paste0, as.list(ref_order)))
    }
    guess1 <- function(x) {
      # Try to guess looking to indices in x
      # Check if min & max indices are first & last part of this dive cycle
      cnd_min <- which.min(x[ , 1]) == order_divecycle(unique(x$elt))[1]
      cnd_max <- which.max(x[ , 2]) == rev(order_divecycle(unique(x$elt)))[1]
      if (!cnd_min || !cnd_max) {
        return(NULL)
      } else {
        return(data.frame(
          st_idx = min(x[ , 1], na.rm = TRUE), ed_idx = max(x[ , 2], na.rm = TRUE), 
          no_dive = "?", no = unique(x$no), elt = unique(x$grp), 
          elt_mtch = TRUE, grp = unique(x$grp), 
          grp_mtch = 1, stringsAsFactors = FALSE
        ))
      }
    }
    guess2 <- function(x) {
      # Try to guess looking in delim table
      first_elt <- x$elt[order_divecycle(unique(x$elt))[1]]
      last_elt <-  x$elt[rev(order_divecycle(unique(x$elt)))[1]]
      previous_elt <- which(ref_order == first_elt) - 1
      next_elt <- which(ref_order == last_elt) + 1
      previous_elt > 0 && next_elt <= length(ref_order) || return(NULL)
      previous_elt <- ref_order[previous_elt]
      next_elt <- ref_order[next_elt]
      no <- unique(x$no)
      out <- data.frame(
        st_idx = elt_delim(previous_elt, no, delim = obj$delim)[ , 2] + 1, 
        ed_idx = elt_delim(next_elt, no, delim = obj$delim)[ , 1] - 1, 
        no_dive = "?", no = no, elt = unique(x$grp), elt_mtch = TRUE, 
        grp = unique(x$grp), grp_mtch = 1, stringsAsFactors = FALSE
      )
      out
    }
    
    # Try to see if indices can be found (e.g. "!_/" not requires "_" indices)
    if (is.in_divecycle(x)) {
      tmp <- try(guess1(x), silent = TRUE) %else% try(guess2(x), silent = TRUE)
      if (!is.null(tmp) && !is.error(tmp)) x <- tmp
    }
  }
  out <- na.omit(x)
  if (nrow(out) == 0) {
    if (no_match == 'error') stop(msg) else warning(msg)
    return(NULL)
  }
  class(out) <- c(class(out), "ty")
  out
}

#' Expand a summary variable to the length of TDR data
#' 
#' @param x variable to expand.
#' @param ty a kind of regular expression that explains how to 
#' duplicate x values along TDR rows.
#' @param na_value the NA value to be used (0 is can be useful).
#' @param obj a "ses" object
#' @return a vector of the same length as the number of rows in the tdr table of 
#' "obj".
#' @export
#' @import data.table
#' @seealso \code{\link{tdrply}}
#' @examples
#' data(exses)
#' ind(exses)
#' exses$tdr$no_btt <- tdrexpand(exses$stat$no_dive, "_")
#' exses$tdr$no_dive  <- tdrexpand(exses$stat$no_dive, "!_/")
tdrexpand <- function(x, ty = "!_/", na_value = NA, obj = ind()) {
  if (is.character(ty)) {
    length(ty) == 1 || stop('"ty" must be of length = 1')
    !all(c(grepl("-", ty), grepl("~", ty))) || stop('"ty" can not contain "-" and "~"')
    ty <- ty_delim(ty, no = NULL, obj, no_match = "ignore")
    n <- length(ty) 
    ty <- rbindlist(ty)
    ty$mrg <- paste0(ty$grp, ty$no)
  } else if (is.data.frame(ty)) {
    ty$grp <- "cst"
    ty$no <- try(ty[ , 3], silent = TRUE) %else% seq_along(ty$grp)
    ty$mrg <- paste0(ty$grp, ty$no)
    ty <- as.data.table(ty)
  } else {
    stop('"ty" must be a character or a data frame.')
  }
  n == length(x) || stop('"x" has an incorrect length')
  ty <- merge(ty, data.table(mrg = unique(ty$mrg), val = x), by = "mrg")
  ty <- as.data.frame(ty)
  
  # Pre-allocate memory
  out <- data.table(as(rep(na_value, nrow(obj$tdr)), typeof(x)))
  
  # Replace content
  for (ii in seq(1, nrow(ty))) {
    idx <- seq(ty[ii, 2], ty[ii, 3])
    val <- ty$val[ii]
    out[idx, V1 := val]
  }
  
  as.vector(out$V1)
}

#' Apply function to subsets of a TDR dataset
#' 
#' \code{tdrply} is a functional programing utility to apply functions to specific 
#' parts and variables of TDR datasets. It is based on a call to \code{\link{mapply}}.
#' 
#' @param f function(s) to apply. The function has
#' to be written considering that it's first argument will be a subset of the TDR data
#' (columns according to \code{cl}, rows according to \code{ty}). 
#' @param cl character of numeric, the columns to select in the \code{tdr} table. 
#' Enter \code{cl = NULL} or \code{cl = .} (default) to keep all columns.
#' @param ty a kind of regular expression to subset diving periods of \code{tdr} table. 
#' The pattern can be any part of this dive cycle representation \code{'-!_/~'}:
#' \itemize{
#'  \item \code{'-'} surface preceding the dive
#'  \item \code{'!'} descent of the dive
#'  \item \code{'_'} bottom of the dive
#'  \item \code{'/'} ascent of the dive
#'  \item \code{'~'} surface following the dive
#' }
#' Several symbols can be juxtaposed to bluid more complex groups e.g. 
#' \code{ty = "!_/"} (default which represents a dive). \code{tdrply} return the 
#' result of \code{f} for each group. \code{'&'} operator can be used to provide 
#' several groups in a same \code{ty} expression. \code{'|'} operator can be used 
#' within a group in order to match a group in priority but provide rescue cases 
#' (e.g \code{'-|~'} matches \code{'-'} first but will also match \code{'~'} 
#' if \code{'-'} is not found). \code{'()'} are implemented so \code{'(-|~_)'} is 
#' different from \code{'(-|~)_'}). 
#' @param no the dive numbers to process. Keep default \code{no = NULL} for all 
#' dives available in the \code{delim} table of \code{obj}. Negative values can 
#' be used to exlude dives; \code{-0} syntax is accepeted to remove the 
#' dive number 0.
#' @param la list of additional arguments of \code{f} whose values depend 
#' on the period involved.
#' @param no_match how to handle no match or partial match. See last example for 
#' details.
#' @param obj a "ses" object.Optional if a default individual has been declared 
#' with \code{\link{ind}}.
#' @param ... additional arguments to be passed to \code{f}. Arguments passed 
#' through \code{...} are recycled for each group defined by \code{ty} and \code{no}.
#' @export
#' @import data.table
#' @seealso \code{\link{tdrexpand}}
#' @examples
#' data(exses)
#' ind(exses)
#' 
#' ## Apply function to each dive
#' tdrply(function(x) max(x$depth, na.rm = TRUE))
#' # or, using "cl" to subset columns and "..." to set "na.rm"
#' tdrply(max, "depth", na.rm = TRUE)
#' # other examples with "cl" and "..."
#' bsm <- tdrply(brokenstick, 1:2, npts = 10)
#' tdrply(plot, 1:2, no = 33) ; plot(bsm[["!_/#33"]], add = TRUE, enumerate = TRUE)
#' 
#' # use "no" to specify what dive numbers should be processed
#' (tmp <- tdrply(max, "depth", no = 111, na.rm = FALSE)) # dive no 111 only
#' tdrply(max, "depth", no = -111, na.rm = FALSE) # all dives but no 111
#' tdrply(max, "depth", no = exses$stat$max_depth == tmp, na.rm = FALSE) # logicals are accepeted
#' 
#' ## How to use "ty": few examples
#' # apply "f" to each group delimited by &
#' tdrply(max, "light", ty = "!&_&/", no = 50:51)
#' # apply "f" previous surface or, if not found, next surface 
#' tdrply(max, "light", ty = "-|~", no = 50:51)
#' # mix operators and use parentheses as desired
#' tdrply(max, "light", ty = "_!&_/&_(-|~)", no = 50:51)
#' # When "ty" syntax is not enought just give a data frame instead
#' df <- data.frame(start = seq(1, 5000, 1000), 
#'                  end = seq(1001, 6000, 1000))
#' tdrply(max, "depth", df) # notice that names start with "cst" as "custom"
#' # An id for output names can be provided in third column
#' (df$id <- sample(10:20, 5, replace = FALSE))
#' tdrply(max, "depth", df)
#' 
#' ## Provide different arguments to each group with "la" (List of Arguments)
#' opar <- par(no.readonly = TRUE); par(mfrow = c(2, 2))
#' tdrply(plot, 1:2, no = 50:53, la = list(col = 1:4))
#' par(opar)
#' 
#' ## tdrply is vectorized over "f"
#' funs <- c(min_depth = min, max_depth = max)
#' tmp <- tdrply(funs, "depth", "_", no = 50:53, na.rm = TRUE)
#' as.data.frame(tmp)
#' 
#' ## tdrply is vectorized over "ty" as well
#' tmp <- tdrply(funs, "depth", c("!", "/", surf = "-~"), no = 50:53, na.rm = TRUE)
#' as.data.frame(unlist(tmp, recursive = FALSE))
#' 
#' \dontrun{
#' ## Assuming that some dive bottom could not be properly defined
#' exses$delim[exses$delim$no_dive == 111, c("btt_st_idx", "btt_ed_idx")] <- c(NA, NA)
#' # choose how to handle missing periods using "no_match".
#' # a warning is printed anyway
#' tdrply(max, "depth", ty = "_", no = 110:112, no_match = "na") # default
#' tdrply(max, "depth", ty = "_", no = 110:112, no_match = "ignore")
#' tdrply(max, "depth", ty = "_", no = 110:112, no_match = "error")
#' # but not if the missing period can be ignored e.g:
#' tdrply(max, "depth", ty = "!_/", no = 110:112, no_match = "error")
#' }
tdrply <- function(f, cl = ., ty = "!_/", no = NULL, la = NULL, 
                   no_match = "na", obj = ind(), ...) {
  # Vectorize function over "f" and "ty" arguments
  if (is.character(ty) && length(ty) > 1) {
    args <- as.list(match.call(expand.dots = TRUE))
    args <- eval(args[-c(1, 4)], parent.frame())
    out <- lapply(ty, function(x, ...) do.call(tdrply, c(list(ty = x), args)))
    # Set names to output according to ty content of name
    nms <- sapply(ty, function(x) switch(x, "-"="psf","!"="dsc","_"="btt","/"="asc","~"="nsf","!_/"="dv","!/"="trs", x))
    nms <- "if"(is.null(names(ty)), nms, mapply(function(x, y) "if"(x == "", y, x), names(ty), nms))
    return(setNames(out, nms))
  }
  if (is.list(f) || length(f) > 1) {
    args <- as.list(match.call(expand.dots = TRUE))
    args <- eval(args[-(1:2)], parent.frame())
    return(lapply(f, function(x, ...) do.call(tdrply, c(list(f = x), args))))
  }
  
  no <- no_interpreter(no, obj)
  cl <- "if"(deparse(substitute(cl)) == "." || is.null(cl), names(obj$tdr), cl)
  ty <- "if"(is.character(ty), ty_delim(ty, no, obj, no_match = no_match), ty)
  
  if (is.null(ty)) {
    switch(no_match, ignore = return(NULL), na = return(NA), error = stop('Empty "ty" found (no match).'))
  } else if (is.data.frame(ty)) {
    data <- Map(function(st, ed) identity(obj$tdr[st:ed, cl]), ty[ , 1], ty[ , 2])
    if (is(ty, "ty")) {
      # Concatenate subset of TDR
      data <- "if"(list_depth(data) == 1, unlist(data), as.tdr(rbindlist(data)))
      return(f(data, ...))
    } else {
      # Loop on subset of TDR
      args <- c(list(FUN = f, data, SIMPLIFY = FALSE), 
                if (is.null(la)) c() else la, 
                if (is.null(list(...))) c() else list(MoreArgs = list(...)))
      out <- do.call(mapply, args)
    }
  } else if (is.list(ty) && !is.data.frame(ty)) {
    # Apply "tdrply" to each element of list
    loop_f <- function(ty, ...) tdrply(f, cl, ty, no, ..., no_match = no_match, obj = obj)
    args <- c(list(FUN = loop_f, ty = ty, SIMPLIFY = FALSE), 
              if (is.null(la)) c() else la, 
              if (is.null(list(...))) c() else list(MoreArgs = list(...)))
    out <- do.call(mapply, args)
  } else {
    stop('"ty" must be a data frame or a list (of any depth) containg data frames.')
  }
  
  # Set names and order in output
  if (is.null(names(out))) {
    # "ty" argument was a data.frame
    nms <- paste0("cst#", try(ty[ , 3], silent = TRUE) %else% seq_along(ty[ , 1]))
    out <- setNames(out, nms)
  } else {
    # ty_delim was used: set order by "no" and "ty" groups
    df <- as.data.frame(do.call(rbind, strsplit(names(out), "#")))
    df <- df[sapply(no, function(x) which(df$V2 == x)), ]
    out <- out[paste0(df$V1, "#", df$V2)]
  }
  
  # Simplify output
  "if"(all(sapply(out, is.null)), return(invisible(out)))
  out <- "if"(any(cnd <- sapply(out, is.null)), out[!cnd], out)
  "if"(all(!sapply(out, is.recursive)) && nUN(sapply(out, length)) == 1, unlist(out, recursive = FALSE), out)
}
