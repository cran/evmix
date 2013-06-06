#' @export
#' @aliases lbckdengpd nlbckdengpd
#' @family  lbckdengpd

# negative cross-validation log-likelihood function for boundary corrected kernel density estimators GPD 
# tail extreme value mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlbckdengpd <- function(pvector, x, phiu = TRUE, finitelik = FALSE, 
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = 0, xmax = Inf) {
  
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
  
  lambda = pvector[1]
  u = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]
  
  nllh = -lbckdengpd(x, lambda, u, sigmau, xi, phiu, bcmethod, proper, nn, offset, xmax, log = TRUE) 

  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}

