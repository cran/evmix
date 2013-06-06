#' @export
#' @aliases  gpd dgpd pgpd qgpd rgpd
#' @rdname gpd

# cumulative distribution function for GPD
pgpd <- function(q, u = 0, sigmau = 1, xi = 0, phiu = 1, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
    
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")

  q[is.infinite(q)]=NA # user will have to deal with infinite cases

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(q), length(u), length(sigmau), length(xi), length(phiu))
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

  q = rep(q, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)
  
  yu = pmax(q - u, 0) / sigmau # used when shape is zero
  syu = pmax(1 + xi*yu, 0)     # used when shape non-zero
  
  # check for x values in range
  yind = (yu > 0) 

  p = q # this will pass through NA/NaN in q just as they are entered
  p[which(!yind)] = ifelse(phiu[which(!yind)] == 1, 0, NA) # NA is default, but if only GPD then gives 0 below threshold
  
  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)
    
  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    p[whichexp] = 1 - exp(-yu[whichexp])
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    p[whichxi] = 1 - syu[whichxi] ^ (-1 / xi[whichxi])
  }
  
  p = 1 - (1 - p) * phiu
  
  if (!lower.tail) p = 1 - p
  
  p
}

