#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname bckdengpd

# inverse cumulative distribution function for boundary corrected kernel density estimators
# for the bulk distribution upto the threshold and conditional GPD above threshold.
qbckdengpd <- function(p, kerncentres, lambda = NULL,
  u = as.vector(quantile(kerncentres, 0.9)), sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, 
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
  
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")
  
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
  
  if (!lower.tail) p = 1 - p
  
  if (is.logical(phiu)) {
    phiu = 1 - pbckden(u, kerncentres, lambda = lambda, bcmethod = bcmethod,
      proper = proper, nn = nn, offset = offset, xmax = xmax)
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pbckden(u, kerncentres, lambda = lambda, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax)
  
  q = p # this will pass through NA/NaN in p just as they are entered
  
  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)
  
  if (nb > 0) {
    q[whichb] = qbckden(p[whichb] / phib, kerncentres, lambda = lambda, bcmethod = bcmethod,
      proper = proper, nn = nn, offset = offset, xmax = xmax)
  }
  if (nu > 0) q[whichu] = qgpd(p[whichu], u, sigmau, xi, phiu)
  
  q
  
}

