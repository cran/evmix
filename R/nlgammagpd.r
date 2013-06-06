#' @export
#' @aliases lgammagpd nlgammagpd
#' @rdname lgammagpd

# negative log-likelihood function for gamma bulk with GPD for upper tail
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgammagpd <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

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

  gshape = pvector[1]
  gscale = pvector[2]
  u = pvector[3]
  sigmau = pvector[4]
  xi = pvector[5]

  nllh = -lgammagpd(x, gshape, gscale, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
