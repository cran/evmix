#' @export
#' @aliases  gammagpdcon dgammagpdcon pgammagpdcon qgammagpdcon rgammagpdcon
#' @rdname gammagpdcon

# inverse cumulative distribution function for gamma bulk with GPD for upper tail
qgammagpdcon <- function(p, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 1/gscale),
  xi = 0, phiu = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
  
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(p), length(gshape), length(gscale), 
              length(u), length(xi), length(phiu))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
  if (mode(gshape) != "numeric" | mode(gscale) != "numeric" | mode(u) != "numeric" |
    mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(gshape, gscale, u, xi, phiu))))
    stop("parameters must be numeric")
  
  if (min(gshape) < 0)
    stop("gamma shape must be non-negative")
  
  if (min(gscale) < 0)
    stop("gamma scale must be non-negative")
  
  if (min(u) <= 0)
    stop("threshold must be non-negative")
  
  if (is.logical(phiu) & any(!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiu < 0) | any(phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }
  
  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  if (!lower.tail) p = 1 - p
  
  p = rep(p, length.out = n)
  gshape = rep(gshape, length.out = n)
  gscale = rep(gscale, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  if (is.logical(phiu)) {
    phiu = 1 - pgamma(u, shape = gshape, scale = gscale)
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pgamma(u, shape = gshape, scale = gscale)
  
  sigmau = phiu / (phib * dgamma(u, shape = gshape, scale = gscale))
  
  if (any(!is.finite(sigmau)))
    stop("sigmau is not numeric")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  qgammagpd(p, gshape, gscale, u, sigmau, xi, phiu, lower.tail)
}
