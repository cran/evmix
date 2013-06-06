#' @export
#' @aliases  lgpd nlgpd
#' @rdname lgpd

# negative log-likelihood function for GPD
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgpd <- function(pvector, x, u = 0, phiu = 1, finitelik = FALSE) {

  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
    
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")

  if (length(pvector) != 2)
    stop("Two parameters must be specified")
  
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")

  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  sigmau = pvector[1]
  xi = pvector[2]

  nllh = -lgpd(x, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}

