#' @export
#' @aliases  betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @rdname betagpd

# random number generation for beta bulk with GPD for upper tail
rbetagpd <- function(n = 1, bshape1 = 1, bshape2 = 1, u = qbeta(0.9, bshape1, bshape2),
  sigmau = sqrt(bshape1 * bshape2 / (bshape1 + bshape2)^2 / (bshape1 + bshape2 + 1)),
  xi = 0, phiu = TRUE) {

  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(bshape1), length(bshape2), 
    length(u), length(sigmau), length(xi), length(phiu))

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")

  if (mode(bshape1) != "numeric" | mode(bshape2) != "numeric" | mode(u) != "numeric" |
    mode(sigmau) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(bshape1, bshape2, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (min(bshape1) < 0)
    stop("beta shape 1 must be non-negative")

  if (min(bshape2) < 0)
    stop("beta shape 2 must be non-negative")
  
  if (min(u) <= 0 | max(u) >= 1)
    stop("threshold must be between 0 and 1 (exclusive)")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")

  if (any(xi == 1))
    stop("shape cannot be 1")

  if (is.logical(phiu) & any(!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiu < 0) | any(phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }

  qbetagpd(runif(n), bshape1, bshape2, u, sigmau, xi, phiu)
}


