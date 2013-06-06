#' @export
#' @aliases  lgngcon nlgngcon
#' @rdname lgngcon

# negative log-likelihood function for normal bulk with GPD's for upper and lower tails with Continuity Constraints
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgngcon <- function(pvector, x, phiul = TRUE, phiur = TRUE, finitelik = FALSE) {

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

  nmean = pvector[1]
  nsd = pvector[2]
  ul = pvector[3]
  xil = pvector[4]
  ur = pvector[5]
  xir = pvector[6]

  nllh = -lgngcon(x, nmean, nsd, ul, xil, phiul, ur, xir, phiur) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
                                               