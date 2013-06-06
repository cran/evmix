#' @export
#' @aliases kdengpdcon dkdengpdcon pkdengpdcon qkdengpdcon rkdengpdcon
#' @rdname kdengpdcon

# inverse cumulative distribution function for kernel density estimation using normal kernel
# for the bulk distribution upto the threshold and conditional GPD above
# threshold with a single continuity constraint.
qkdengpdcon <- function(p, kerncentres, lambda = NULL, 
  u = as.vector(quantile(kerncentres, 0.9)), xi = 0, phiu = TRUE, lower.tail = TRUE) {
  
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
  
  if (is.null(lambda)){
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10){
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    lambda = bw.nrd0(kerncentres)
  }
  
  linputs = c(length(lambda), length(u), length(xi), length(phiu))

  if (sum(linputs != 1) > 0)
    stop("parameters must be scalar")
  
  if (mode(lambda) != "numeric" | mode(u) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")

  if (any(!is.finite(c(lambda, u, xi, phiu))))
    stop("parameters must be numeric")  
  
  if (lambda <= 0)
    stop("bandwidth must be non-negative")  
    
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
    phiu = 1 - pkdenx(u, kerncentres, lambda)
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pkdenx(u, kerncentres, lambda)
  
  du = kdenx(u, kerncentres, lambda)
  sigmau = phiu / (phib * du)
  
  if (!is.finite(sigmau))
    stop("sigmau is not numeric")

  if (sigmau <= 0)
    stop("scale must be non-negative")

  q = p # this will pass through NA/NaN in p just as they are entered

  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)
  
  if (nb > 0) q[whichb] = qkden(p[whichb] / phib, kerncentres, lambda)
  if (nu > 0) q[whichu] = qgpd(p[whichu], u, sigmau, xi, phiu)
   
  q
  
}
