#' @export
#' @aliases  betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @rdname betagpd

# inverse cumulative distribution function for beta bulk with GPD for upper tail
qbetagpd <- function(p, bshape1 = 1, bshape2 = 1, u = qbeta(0.9, bshape1, bshape2),
  sigmau = sqrt(bshape1 * bshape2 / (bshape1 + bshape2)^2 / (bshape1 + bshape2 + 1)),
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
  linputs = c(length(p), length(bshape1), length(bshape2), 
    length(u), length(sigmau), length(xi), length(phiu))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

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
  bshape1 = rep(bshape1, length.out = n)
  bshape2 = rep(bshape2, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  if (is.logical(phiu)) {
    phiu = 1 - pbeta(u, shape1 = bshape1, shape2 = bshape2)
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pbeta(u, shape1 = bshape1, shape2 = bshape2)
    
  q = p # this will pass through NA/NaN in p just as they are entered
  
  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)

  if (nb > 0) {q[whichb] = qbeta(p[whichb] / phib[whichb], shape1 = bshape1[whichb], shape2 = bshape2[whichb])}
  if (nu > 0) {q[whichu] = qgpd(p[whichu], u[whichu], sigmau[whichu], xi[whichu], phiu[whichu])}

  q
}

