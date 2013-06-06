#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname gkg

# inverse cumulative distribution function for KDE for bulk with GPD's for both upper and lower tails
qgkg <- function(p, kerncentres, lambda = NULL, 
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE, 
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
  
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")
  
  if (missing(kerncentres))
    stop("kerncentres must be a non-empty numeric vector")
  
  if (length(kerncentres) == 0 | mode(kerncentres) != "numeric") 
    stop("kerncentres must be a non-empty numeric vector")
  
  if (any(!is.finite(kerncentres)))
    warning("non-finite kernel centres are dropped")
  
  kerncentres = kerncentres[is.finite(kerncentres)]
  nk = length(kerncentres)
  
  if (is.null(lambda)){
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10){
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    lambda = bw.nrd0(kerncentres)
  }
  
  linputs = c(length(lambda), length(ul), length(sigmaul), length(xil), length(phiul),
                              length(ur), length(sigmaur), length(xir), length(phiur)) 
  
  if (sum(linputs != 1) > 0)
    stop("parameters must be scalar")
  
  if (mode(lambda) != "numeric"| 
      mode(ul) != "numeric" | mode(sigmaul) != "numeric" | mode(xil) != "numeric" |
      mode(ur) != "numeric" | mode(sigmaur) != "numeric" | mode(xir) != "numeric" |
      !(mode(phiul) %in% c("logical","numeric")) | !(mode(phiur) %in% c("logical","numeric")))
    stop("parameters must be numeric, phiu can be numeric or logical")

  if (any(!is.finite(c(lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur))))
    stop("parameters must be numeric")  
  
  if (lambda <= 0)
    stop("bandwidth must be non-negative")  
  
  if ((sigmaul <= 0) | (sigmaur <= 0))
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
  
  if (is.logical(phiul)) {
    phiul = pkdenx(ul, kerncentres, lambda)
  } else {
    phiul = phiul
  }
  if (is.logical(phiur)) {
    phiur = 1 - pkdenx(ur, kerncentres, lambda)
  } else {
    phiur = phiur
  }
  phib = (1 - phiul - phiur) / (pkdenx(ur, kerncentres, lambda) - pkdenx(ul, kerncentres, lambda))
  
  q = p # this will pass through NA/NaN in p just as they are entered
  
  whichul = which(p < phiul)
  nul = length(whichul)
  whichb = which((p <= (1 - phiur)) & (p >= phiul)) 
  nb = length(whichb)
  whichur = which(p > (1 - phiur))
  nur = length(whichur)
  
  if (nul > 0) q[whichul] = -qgpd(1 - p[whichul], -ul, sigmaul, xil, phiul)
  if (nb > 0) q[whichb] = qkden((p[whichb] - phiul) / phib + mean(pnorm(ul, kerncentres, lambda)), kerncentres, lambda)
  if (nur > 0) q[whichur] = qgpd(p[whichur], ur, sigmaur, xir, phiur)
    
  q
}

