#' @export
#' @aliases  gpd dgpd pgpd qgpd rgpd
#' @rdname gpd

# random number generation for GPD
rgpd <- function(n = 1, u = 0, sigmau = 1, xi = 0, phiu = 1) {

  # Check properties of inputs
  if (length(n) != 1 || mode(n) != "numeric") 
    stop("sample size must be positive integer")

  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 0)
    stop("sample size must be positive integer")

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(u), length(sigmau), length(xi), length(phiu))

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")

  if (mode(u) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric" | 
    mode(phiu) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")

  if (any(xi == 1))
    stop("shape cannot be 1")

  if (any(phiu < 0) | any(phiu > 1))
    stop("phiu must between 0 and 1 (inclusive)")
    
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)

  yind = runif(n) < phiu
  
  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)
  
  r = rep(NA, n)
  
  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    r[whichexp] = u[whichexp] + sigmau[whichexp] * rexp(u[whichexp])
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind  & yind)
    r[whichxi] = u[whichxi] + sigmau[whichxi] * (runif(length(whichxi)) ^ (-xi[whichxi]) - 1) / xi[whichxi]
  }
    
  r
}
