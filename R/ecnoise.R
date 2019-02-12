#' Run the ECnoise algorithm on a provided difference stencil
#'
#' Reference: More' and Wild, "Estimating computational noise"
#' (2011) SIAM J. Sci. Comput.
#'
#' Implemented by Zach del Rosario, Feb. 2019
#'
#' @param fval vector of function evaluations
#' @return list() bundling outputs:
#'         $fnoise = scalar noise level
#'         $level  = vector of differences
#'         $inform = string stating algorithm results
#' @export
#' @examples
#'

ecnoise <- function(fval) {
  ## Setup
  nf     = length(fval)
  dsgn   = rep(0, nf - 1)
  gamma  = 1

  result = list(
    "fnoise" = 0,
    "level"  = rep(0, nf - 1),
    "inform" = "noise detected"
  )

  ## Compute range
  fmin = min(fval)
  fmax = max(fval)

  ## Check for h too large
  if ((fmax - fmin) / max(abs(fmax), abs(fmin)) > 0.1) {
    result$inform = "h too large; try h / 100 next (early)"
    return(result)
  }

  ## Construct difference table
  for (j in 0:(nf - 3)) {
    for (i in 0:(nf - 2 - j)) {
      fval[i + 1] = fval[i + 2] - fval[i + 1]
    }

    ## h is too small only when half the function values are equal
    if ((j == 0) & (sum(fval == 0) >= (nf / 2))) {
      result$inform = "h too small; try h * 100 next"
      return(result)
    }
    gamma = 0.5 * ( (j + 1) / (2 * (j + 1) - 1) ) * gamma

    ## Compute the estimates for the noise level
    result$level[j + 1] = sqrt(gamma * mean(fval[1:(nf - j)]^2))

    ## Determine differences in sign
    emin = min(fval[1:(nf - j)])
    emax = max(fval[1:(nf - j)])
    if (emin * emax < 0.0) {
      dsgn[j + 1] = 1
    }
  }

  ## Determine the noise level
  for (k in 1:(nf - 5)) {
    emin = min(result$level[(k + 1):(k + 3)])
    emax = max(result$level[(k + 1):(k + 3)])

    if ((emax <= 4 * emin) & (dsgn[k + 1] == 1)) {
      result$fnoise = result$level[k + 1]
      result$inform = "noise detected"
      return(result)
    }
  }

  ## If noise not detected, then h is too large
  result$inform = "h too large; try h / 100 next (end)"
  return(result)
}
