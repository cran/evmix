#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# random number generation for kernel density estimator using normal kernel
rkden <- function(n = 1, kerncentres, lambda = NULL) {

  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")

  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
  
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

  rnorm(rep(1, n), sample(kerncentres, n, replace = TRUE), sd = lambda)
}
