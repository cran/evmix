#' @export
#' @aliases  gpd dgpd pgpd qgpd rgpd
#' @rdname gpd

# inverse cumulative distribution function for GPD
qgpd <- function(p, u = 0, sigmau = 1, xi = 0, phiu = 1, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
    
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(p), length(u), length(sigmau), length(xi), length(phiu))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

  if (mode(u) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric" | 
    mode(phiu) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")

  if (any(phiu < 0) | any(phiu > 1))
    stop("phiu must between 0 and 1 (inclusive)")

  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")

  if (lower.tail) p = 1 - p
  
  p = rep(p, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)
  
  # check for x values in range
  yind = (p < phiu) 

  q = p # this will pass through NA/NaN in p just as they are entered
  q[which(!yind)] = ifelse(phiu[which(!yind)] == 1, u[which(!yind)], NA) # NA is default, but if only GPD then gives threshold

  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)
    
  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    q[whichexp] = u[whichexp] - sigmau[whichexp] * log(p[whichexp] / phiu[whichexp])
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    q[whichxi] = u[whichxi] + sigmau[whichxi] * ((p[whichxi] / phiu[whichxi]) ^ (-xi[whichxi]) - 1) / xi[whichxi]
  }
    
  q
}
