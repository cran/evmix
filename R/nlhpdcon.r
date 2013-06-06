#' @export
#' @aliases  lhpdcon nlhpdcon
#' @rdname lhpdcon

# negative log-likelihood function for hybrid Pareto extreme value mixture model
# with a single continuity constraint
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlhpdcon <- function(pvector, x, finitelik = FALSE) {
  
  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
  
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")
  
  if (length(pvector) != 4)
    stop("Four parameters must be specified")
  
  if (missing(x))
    stop("x must be a non-empty numeric vector")                              
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")
  
  nmean = pvector[1]
  nsd = pvector[2]
  u = pvector[3]
  xi = pvector[4]
    
  nllh = -lhpdcon(x, nmean, nsd, u, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}
