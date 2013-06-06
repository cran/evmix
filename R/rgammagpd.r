#' @export
#' @aliases  gammagpd dgammagpd pgammagpd qgammagpd rgammagpd
#' @rdname gammagpd

# random number generation for gamma bulk with GPD for upper tail
rgammagpd <- function(n = 1, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 1/gscale),
  sigmau = sqrt(gshape) * gscale, xi = 0, phiu = TRUE) {

  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(gshape), length(gscale), 
    length(u), length(sigmau), length(xi), length(phiu))

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")

  if (mode(gshape) != "numeric" | mode(gscale) != "numeric" | mode(u) != "numeric" |
    mode(sigmau) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(gshape, gscale, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (min(gshape) < 0)
    stop("gamma shape must be non-negative")

  if (min(gscale) < 0)
    stop("gamma scale must be non-negative")
  
  if (min(u) <= 0)
    stop("threshold must be non-negative")

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

  qgammagpd(runif(n), gshape, gscale, u, sigmau, xi, phiu)
}


