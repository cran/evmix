#' @export
#' @aliases  lgkg nlgkg
#' @rdname lgkg

# negative cross-validation log-likelihood function for KDE for bulk with GPD's for both upper and lower tails
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgkg <- function(pvector, x, phiul = TRUE, phiur = TRUE, finitelik = FALSE) {
  
  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
  
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")
  
  if (length(pvector) != 7)
    stop("Seven parameters must be specified")
  
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")
  
  lambda = pvector[1]
  ul = pvector[2]
  sigmaul = pvector[3]
  xil = pvector[4]
  ur = pvector[5]
  sigmaur = pvector[6]
  xir = pvector[7]
  
  nllh = -lgkg(x, lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}

