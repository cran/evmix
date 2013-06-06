#' @export
#' @aliases kdengpdcon dkdengpdcon pkdengpdcon qkdengpdcon rkdengpdcon
#' @rdname kdengpdcon

# random number generation for kernel density estimation using normal kernel
# for the bulk distribution upto the threshold and conditional GPD above
# threshold with a single continuity constraint.
rkdengpdcon <- function(n = 1, kerncentres, lambda = NULL,
  u = as.vector(quantile(kerncentres, 0.9)), xi = 0, phiu = TRUE) {
  
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

  qkdengpd(runif(n), kerncentres, lambda, u, sigmau, xi, phiu)
}



