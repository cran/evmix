#' @export
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @rdname gngcon

# random number generation for normal bulk with GPD's for upper and lower tails
rgngcon <- function(n = 1, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE) {
  
  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")                     
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(nmean), length(nsd),
              length(ul), length(xil), length(phiul),
              length(ur), length(xir), length(phiur))                                             
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" |
    mode(ul) != "numeric" | mode(xil) != "numeric" |
    mode(ur) != "numeric" | mode(xir) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, ul, xil, phiul, ur, xir, phiur))))
    stop("parameters must be numeric")
  
  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")
  
  if (any(ul >= ur))
    stop("lower threshold must be below upper threshold")
  
  if (any(xil == 1) | any(xir == 1))
    stop("shape cannot be 1")
  
  if ((is.logical(phiul) & any(!phiul)) | (is.logical(phiur) & any(!phiur))) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiul < 0) | any(phiul > 1) | any(phiur < 0) | any(phiur > 1))
      stop("phiu must between 0 and 1 (inclusive)")
    if (!is.logical(phiul) & !is.logical(phiur)) {
      if (any((phiul + phiur) > 1))
        stop("phiu + phiur must be less than 1")
    }
  }
  
  qgngcon(runif(n), nmean, nsd, ul, xil, phiul, ur, xir, phiur)
}


