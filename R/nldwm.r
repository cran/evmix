#' @export
#' @aliases ldwm nldwm
#' @rdname  ldwm

# negative log-likelihood function for dynamically weighted mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nldwm = function(pvector, x, finitelik = FALSE) {
  
  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
  
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")
   
  if (length(pvector) != 6)
    stop("Six parameters must be specified")
  
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")
  
  wshape = pvector[1]
  wscale = pvector[2]
  cmu = pvector[3]
  ctau = pvector[4]
  sigmau = pvector[5]
  xi = pvector[6]
  
  nllh = -ldwm(x, wshape, wscale, cmu, ctau, sigmau, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}
