#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname  bckdengpd

# cumulative distribution function for boundary corrected kernel density estimators for the bulk
# distribution upto the threshold and conditional GPD above threshold.
pbckdengpd <- function(q, kerncentres, lambda = NULL,
  u = as.vector(quantile(kerncentres, 0.9)), sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, 
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
  
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")
  
  q[is.infinite(q)]=NA # user will have to deal with infinite cases
  
  if (missing(kerncentres))
    stop("kerncentres must be a non-empty numeric vector")
  
  if (length(kerncentres) == 0 | mode(kerncentres) != "numeric") 
    stop("kerncentres must be a non-empty numeric vector")
  
  if (any(!is.finite(kerncentres)))
    warning("non-finite kernel centres are dropped")
  
  kerncentres = kerncentres[is.finite(kerncentres)]
  nk = length(kerncentres)
  
  if (any(kerncentres < 0))
    stop("kernel centres cannot be non-positive")
  
  if (is.null(lambda))
    stop("bandwidth (lambda) must be specified")
  
  linputs = c(length(lambda), length(u), length(sigmau), length(xi), length(phiu))
  
  if (sum(linputs != 1) > 0)
    stop("parameters must be scalar")
  
  if (mode(lambda) != "numeric" | mode(u) != "numeric" | mode(sigmau) != "numeric" |
      mode(xi) != "numeric")
    stop("parameters must be numeric")

  if (any(!is.finite(c(lambda, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")  
  
  if (lambda <= 0)
    stop("bandwidth must be non-negative")  
  
  if (sigmau <= 0)
    stop("scale must be non-negative")
  
  if (is.logical(phiu) & (!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
         or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if ((phiu < 0) | (phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }
  
  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  if (is.logical(phiu)) {
    phiu = 1 - pbckden(u, kerncentres, lambda = lambda, bcmethod = bcmethod,
      proper = proper, nn = nn, offset = offset, xmax = xmax)
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pbckden(u, kerncentres, lambda = lambda, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax)
 
  p = q # this will pass through NA/NaN in q just as they are entered
  
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0) {
    p[whichb] = phib*sapply(q[whichb], FUN = pbckden, kerncentres = kerncentres, lambda = lambda, 
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
  }
  if (nu > 0) p[whichu] = 1 - phiu + phiu*pgpd(q[whichu], u, sigmau, xi)
  
  if (!lower.tail) p = 1 - p
  
  p
}
