#' @export
#' @aliases  hpd dhpd phpd qhpd rhpd
#' @rdname hpd

# inverse cumulative distribution function for hybrid Pareto model
qhpd <- function(p, nmean = 0, nsd = 1, xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
  
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(p), length(nmean), length(nsd), length(xi))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, xi))))
    stop("parameters must be numeric")
  
  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")  
  
  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  if (!lower.tail) p = 1 - p
  
  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  xi = rep(xi, length.out = n)
  
  z = (1+xi)^2/(2*pi)
  wz = lambertW(z)
  
  u = nmean + nsd * sqrt(wz)
  sigmau = nsd * abs(1+xi) / sqrt(wz)
  
  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  r = 1 + pnorm(u, mean = nmean, sd = nsd)
  
  q=p # this will pass through NA/NaN in q just as they are entered
  phi=pnorm((u-nmean)/nsd)
  phiu=(phi)/(1+phi)
  
  whichb = which(q <= phiu)
  nb = length(whichb)
  whichu = which(q > phiu)
  nu = length(whichu)
  
  if (nb>0) {q[whichb] = qnorm((1 + phi[whichb]) * q[whichb], nmean[whichb], nsd[whichb]) }
  if (nu>0) {q[whichu] = qgpd((1 + phi[whichu]) * q[whichu] - phi[whichu], 0, sigmau[whichu], xi[whichu])+ u[whichu]}
    
  q
  
}



