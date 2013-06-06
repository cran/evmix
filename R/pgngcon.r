#' @export
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @rdname gngcon

# cumulative distribution function for normal bulk with GPD's for upper and lower tails
pgngcon <- function(q, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
    
  if (length(q) == 0 | mode(q) != "numeric")                                  
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")

  q[is.infinite(q)]=NA # user will have to deal with infinite cases

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(q), length(nmean), length(nsd),
    length(ul), length(xil), length(phiul),
    length(ur), length(xir), length(phiur))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

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

  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")

  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  ul = rep(ul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  xir = rep(xir, length.out = n)
  
  if (is.logical(phiul)) {
    phiul = pnorm(ul, mean = nmean, sd = nsd)
  } else {
    phiul = rep(phiul, length.out = n)
  }
  if (is.logical(phiur)) {
    phiur = 1 - pnorm(ur, mean = nmean, sd = nsd)
  } else {
    phiur = rep(phiur, length.out = n)
  }
  phib = (1 - phiul - phiur) / (pnorm(ur, mean = nmean, sd = nsd) - pnorm(ul, mean = nmean, sd = nsd))
  
  sigmaul = phiul / (phib * dnorm(ul, mean = nmean, sd = nsd))
  sigmaur = phiur / (phib * dnorm(ur, mean = nmean, sd = nsd)) 
  
  if (any(!is.finite(sigmaul)) | any(!is.finite(sigmaur)))
    stop("one the GPD scale parameters is not numeric")

  if (min(sigmaul) <= 0| min(sigmaur) <= 0)
    stop("scale must be non-negative")
    
  pgng(q, nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur, lower.tail)
  
}


