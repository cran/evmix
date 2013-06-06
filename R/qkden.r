#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# inverse cumulative distribution function for kernel density estimator using normal kernel
qkden <- function(p, kerncentres, lambda = NULL, lower.tail = TRUE) {

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

  if (is.unsorted(kerncentres)) kerncentres = sort(kerncentres)
  
  if (is.null(lambda)) {
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10){
        stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    lambda = bw.nrd0(kerncentres)
  }

  if (length(lambda) != 1)
    stop("bandwidth must be scalar")

  if ((mode(lambda) != "numeric") & !is.finite(lambda))
    stop("bandwidth must be numeric")
  
  if (lambda <= 0)
      stop("bandwidth must be non-negative")  

  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")

  if (!lower.tail) p = 1 - p

  q = p
  
  # obtain quantile function but interpolation between kernel based CDF estimates
  # - CDF is quick to estimate, interpolation algorithms also quite fast
  # - an alternative solution is to numerical solve to find quantile, but is much slower
  qrange = range(kerncentres) + 5 * lambda * c(-1, 1)
  qk = seq(qrange[1], qrange[2], length.out = min(50 * nk, 1000))

  pk = sapply(qk, FUN = pkdenx, lambda = lambda, kerncentres = kerncentres) 

  qfun = splinefun(x = pk, y = qk)

  whichinterp = which((p >= min(pk)) & (p <= max(pk)))
  ninterp = length(whichinterp)
  which0 = which(p < min(pk))
  n0 = length(which0)
  which1 = which(p > max(pk))
  n1 = length(which1)
  
  if (ninterp > 0) q[whichinterp] = qfun(p[whichinterp])

  # if further out than 5 standard deviations from kernel centres then set:
  # lower tail to -Inf and upper tail set to Inf
  if (n0 > 0) q[which0] = -Inf
  if (n1 > 0) q[which1] = Inf
  
  q
}
