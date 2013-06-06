#' @export
#' @aliases lgng nlgng
#' @rdname lgng

# negative log-likelihood function for normal bulk with GPD's for upper and lower tails
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgng <- function(pvector, x, phiul = TRUE, phiur = TRUE, finitelik = FALSE) {

  # Check properties of inputs
  if (missing(pvector))
    stop("pvector must be a non-empty numeric vector")
    
  if (mode(pvector) != "numeric") 
    stop("pvector must be a non-empty numeric vector")

  if (length(pvector) != 8)
    stop("Eight parameters must be specified")

  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")

  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  nmean = pvector[1]
  nsd = pvector[2]
  ul = pvector[3]
  sigmaul = pvector[4]
  xil = pvector[5]
  ur = pvector[6]
  sigmaur = pvector[7]
  xir = pvector[8]

  nllh = -lgng(x, nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
