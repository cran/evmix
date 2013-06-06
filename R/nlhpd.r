#' @export
#' @aliases  lhpd nlhpd
#' @rdname lhpd

# negative log-likelihood function for hybrid Pareto extreme value mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlhpd <- function(pvector, x, finitelik = FALSE) {
  
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
  
  nmean = pvector[1]
  nsd = pvector[2]
  xi = pvector[3]
    
  nllh = -lhpd(x, nmean, nsd, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}
