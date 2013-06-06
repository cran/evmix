#' @export
#' @aliases lkdengpdcon nlkdengpdcon
#' @rdname lkdengpdcon

# negative cross-validation log-likelihood function for kernel density estimation Using normal kernel
# for the bulk distribution upto the threshold and conditional GPD above
# threshold with a single continuity constraint.
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlkdengpdcon <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {
  
  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
  
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")
  
  if (length(pvector) != 3)
    stop("Three parameters must be specified")
  
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")
  
  lambda = pvector[1]
  u = pvector[2]
  xi = pvector[3]
  
  nllh = -lkdengpdcon(x, lambda, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}


