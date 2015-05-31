#' Creates or tests for 'ses'.
#' 
#' @param ... Elements to group into an 'ses' object.
#' @details \code{as.ses} coerces an object to class \code{ses}. An \code{ses} 
#' object should contain:
#' \itemize{
#'   \item name The name of the individual.
#'   \item hash The md5sum hash of the original dataset.
#'   \item tdr The TDR data.
#'   \item stat The dives statistics.
#'   \item delim A \code{data.frame} that stores dive and bottom delimitation results.
#' }
#' @export
#' @examples
#' new.ses <- list(name = "Example", tdr = data.frame(), stat = data.frame())
#' is.ses(new.ses)
#' is.stat(new.ses$stat)
#' 
#' new.ses <- as.ses(new.ses)
#' is.ses(new.ses)
#' is.stat(new.ses$stat)
as.ses <- function(...) {
  args <- list(...)
  if (length(args) == 1 & any(sapply(args[[1]], is.data.frame)))
    args <- args[[1]]
  var <- c('name', 'hash', 'tdr', 'stat', 'delim')
  cls <- list(Ind.id = "character", hash = "character", 
              tdr = c("tdr", "data.frame"), 
              stat = c("statdives", "data.frame"), delim = c("per", "data.frame"))
  out <- setNames(replicate(length(var), list(NULL)), var)
  for (iobj in seq_along(args)) {
    nobj <- names(args)[iobj]
    if (nobj == "") {
      if (is(args[[iobj]], 'tdr')) out$tdr <- args[[iobj]]
      if (is(args[[iobj]], 'statdives')) out$stat <- args[[iobj]]
    } else {
      out[[nobj]] <- args[[iobj]]
      if (any(grepl(nobj, var)))
        class(out[[nobj]]) <- cls[[nobj]]
    }
  }
  class(out) <- c('ses', 'list')
  out
}

#' @rdname as.ses
#' @param x A object to test.
#' @details \code{as.ses} tests if an object belong to the \code{'ses'} class
#' @keywords internal
#' @export
is.ses <- function(x) inherits(x, 'ses')

#' @rdname as.ses
#' @details \code{is.tdr} Test if an object belong to the \code{'tdr'} class
#' @keywords internal
#' @export
is.tdr <- function(x) inherits(x, 'tdr')

#' @rdname as.ses
#' @details \code{is.stat} tests if an object belong to the \code{'statdives'} class
#' @keywords internal
#' @export
is.stat <- function(x) inherits(x, 'statdives')

#' @rdname as.ses
#' @details \code{as.tdr} coerces an object to class \code{tdr}. A TDR table should 
#' contain at least columns time and depth in that order.
#' @keywords internal
#' @export
as.tdr <- function(x) {
  x <- as.data.frame(x)
  class(x) <- c("tdr", "data.frame")
  x
}
