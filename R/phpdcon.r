#' @export
#' @aliases  hpdcon dhpdcon phpdcon qhpdcon rhpdcon
#' @rdname hpdcon

# cumulative distribution function for hybrid Pareto model with a single continuity constraint
phpdcon <- function(q, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), xi = 0, lower.tail = TRUE) {
  
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
  linputs = c(length(q), length(nmean), length(nsd), length(u), length(xi))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" | mode(u) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, u, xi))))
    stop("parameters must be numeric")
  
  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")

  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  sigmau = 1/dnorm(u, mean = nmean, sd = nsd)
  
  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  r = 1 + pnorm(u, mean = nmean, sd = nsd)
  
  p = q  # this will pass through NA/NaN in q just as they are entered
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0) {p[whichb] = pnorm(q[whichb], nmean[whichb], nsd[whichb]) / r[whichb] }
  if (nu > 0) {p[whichu] = (pnorm(u[whichu], nmean[whichu], nsd[whichu]) + pgpd(q[whichu], u[whichu], sigmau[whichu], xi[whichu])) / r[whichu]}
  
  if (!lower.tail) p = 1 - p
  
  p
}
