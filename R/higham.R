#' Implementation of the Higham function
#'
#' A form of f(t) = t^2 corrupted by a user-specified level of computational
#' noise
#'
#' Reference: More' and Wild, "Estimating computational noise"
#' (2011) SIAM J. Sci. Comput.
#'
#' Implemented by Zach del Rosario, Feb. 2019
#'
#' @param t evaluation point
#' #param L number of iterations; more leads to greater noise
#' @return function value
#' @export
#' @examples
#'

higham <- function(t, L = 30) {
  f <- t
  for (i in 1:L) {
    f <- sqrt(f)
  }
  for (i in 1:L) {
    f <- f^2
  }
  f^2
}
