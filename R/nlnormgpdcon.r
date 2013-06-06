#' @export
#' @aliases lnormgpdcon nlnormgpdcon
#' @rdname lnormgpdcon

# negative log-likelihood function for normal bulk with GPD for upper tail with a continuity constraint
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlnormgpdcon <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

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

  nllh = -lnormgpdcon(x, nmean, nsd, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
