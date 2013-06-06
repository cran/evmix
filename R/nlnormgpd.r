#' @export
#' @aliases  lnormgpd nlnormgpd
#' @rdname lnormgpd

# negative log-likelihood function for normal bulk with GPD for upper tail
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlnormgpd <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
    
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")

  if (length(pvector) != 5)
    stop("Five parameters must be specified")

  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")

  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  nmean = pvector[1]
  nsd = pvector[2]
  u = pvector[3]
  sigmau = pvector[4]
  xi = pvector[5]

  nllh = -lnormgpd(x, nmean, nsd, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
