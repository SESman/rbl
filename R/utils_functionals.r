#' Apply a function along a vector or a list
#' 
#' @param obj A list or an atomic vector.
#' @param fun A function to be applied.
#' @param w Width of the window.
#' @param ... Optional arguments to \code{fun}.
#' @param wty Window type: to choose between \code{c('moving', 'fixed'))}.
#' @param simplify Logical or character string. Should the result be simplified 
#' to a vector, matrix or higher dimensional array if possible? 
#' @param aty Argument type: How the elements of \code{obj} should be given to \code{fun}. 
#' \code{'single'} if they should be regouped in one argument. \code{'multiple'} 
#' if they should be considered as separated arguments. See the last example.
#' @param nas Should NAs be added to the result to match the length of the input.
#' @export
#' @keywords internal
#' @examples
#' x <- seq(1, 3, length = 150) + rt(150, df = 2) / 3
#' plot(x)
#' av  <- rollapply(x, mean, 15)
#' lines(av, lwd = 2)
#' err <- rollapply(x, sd, 15)
#' lines(av + err, lwd = 2, lty = 2) ; lines(av - err, lwd = 2, lty = 2)
#' 
#' # 'wty' argument
#' lines(rollapply(x, mean, wty = 'f'), type = 's', col = 'red')
#' 
#' # 'aty' argument
#' x <- 2^(0:5)
#' out <- rollapply(x, seq, 2, aty = 'm')
#' names(out) <- rollapply(x, paste, 2, sep = ':', aty = 'm')
#' out
rollapply <- function(obj, fun, w = 5, ..., wty = c('moving', 'fixed'), simplify = TRUE, 
                      aty = c('single', 'multiple'), nas = pmatch(match.arg(aty), 'single', 0)) {
  # Make intervals
  int <- lapply(seq_along(obj), function(n) seq(w) + n - 1)
  int <- int[sapply(int, function(x) all(x <= length(obj)))]
  if (match.arg(wty) == 'fixed')
    int <- int[seq(from = 1, to = length(int), by = w)]
  
  # Make the work
  .f <- switch(match.arg(aty), 
               single = function(int, ...) fun(unlist(obj[int]),...), 
               multiple = function(int, ...) do.call(fun,c(as.list(obj[int]),list(...))))
  out <- sapply(int, .f, simplify = simplify, ...)
  
  # Final formatting
  if (match.arg(wty) == 'fixed')
    out <- rep(out, each = w)
  if (nas)
    out <- c(rep(NA, (w - 1) %/% 2), out, rep(NA, w %/% 2))
  out
}

#' Curry
#' 
#' @param FUN a function
#' @param ... arguments to be set
#' @export
#' @keywords functional
#' @examples 
#' f <- curry(mean, na.rm = TRUE)
#' f(c(1, NA))
curry <- function(FUN, ...) {
  .orig = list(...);
  function(...) do.call(FUN,c(.orig,list(...)))
}

#' Function composition
#' 
#' @param ... functions to be composed
#' @return \code{compose} returns the composed of functions listed in \code{...}.
#' Pay attention to the order, compose(f, g) returns \code{g(f())}.
#' @export
#' @details see the infix version \code{\link{\%c\%}}
#' @keywords functional
#' @examples 
#' f <- compose(is.na, sum)
#' f(c(1, NA, NA))
compose <- function(...) {
  fs <- list(...)
  function(...) Reduce(function(x, f) f(x), fs, ...)
}

#' Operator for function composition
#' 
#' @param g a function
#' @param f a function
#' @return \code{g \%c\% f} returns the composed \code{g(f())}.
#' @export
#' @keywords internal functional
#' @examples 
#' f <- sum %c% is.na
#' f(c(1, NA, NA))
"%c%" <- function(g, f) compose(f, g)

#' Create a function that return a result if number of valid obs >= n
#' @param f a function having a vector of data as first argument
#' @param n the minimum number of observation required to return a non-NA 
#' result. A number in \code{(0;1[} will be interpreted as a proportion.
#' @export
#' @keywords internal
#' @examples 
#' mean5 <- min_n(mean, 5)
#' mean5(c(1:4, NA))
#' mean5(c(1:5, NA))
#' 
#' mean90percent <- min_n(mean, 0.90)
#' mean90percent(c(1:8, NA, NA))
#' mean90percent(c(1:9, NA))
min_n <- function(f, n) {
  if (n >= 0 & n < 1)
    function(x, ...) "if"(mean(cnd <- is.finite(x)) < n, NA, f(x[cnd], ...))
  else
    function(x, ...) "if"(sum(cnd <- is.finite(x)) < n, NA, f(x[cnd], ...))
}
