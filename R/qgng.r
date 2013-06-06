#' @export
#' @aliases  gng dgng pgng qgng rgng
#' @rdname gng

# inverse cumulative distribution function for normal bulk with GPD's for upper and lower tails
qgng <- function(p, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), sigmaul = nsd, xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), sigmaur = nsd, xir = 0, phiur = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
    
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(p), length(nmean), length(nsd),
    length(ul), length(sigmaul), length(xil), length(phiul),
    length(ur), length(sigmaur), length(xir), length(phiur))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" |
    mode(ul) != "numeric" | mode(sigmaul) != "numeric" | mode(xil) != "numeric" |
    mode(ur) != "numeric" | mode(sigmaur) != "numeric" | mode(xir) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur))))
    stop("parameters must be numeric")

  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")

  if (any(ul >= ur))
    stop("lower threshold must be below upper threshold")

  if (min(sigmaul) <= 0 | min(sigmaur) <= 0)
    stop("scale must be non-negative")
  
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

  if (!lower.tail) p = 1 - p

  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  ul = rep(ul, length.out = n)
  sigmaul = rep(sigmaul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  sigmaur = rep(sigmaur, length.out = n)
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
      
  q = p # this will pass through NA/NaN in p just as they are entered
  
  whichul = which(p < phiul)
  nul = length(whichul)
  whichb = which((p <= (1 - phiur)) & (p >= phiul)) 
  nb = length(whichb)
  whichur = which(p > (1 - phiur))
  nur = length(whichur)

  if (nul > 0) q[whichul] = -qgpd(1 - p[whichul], -ul[whichul], sigmaul[whichul], xil[whichul], phiul[whichul])
  if (nb > 0) q[whichb] = qnorm((p[whichb] - phiul[whichb]) / phib[whichb] + pnorm(ul[whichb], nmean[whichb], nsd[whichb]), nmean[whichb], nsd[whichb])
  if (nur > 0) q[whichur] = qgpd(p[whichur], ur[whichur], sigmaur[whichur], xir[whichur], phiur[whichur])
                  
  q
}


