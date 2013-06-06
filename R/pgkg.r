#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname gkg

# cumulative distribution function for KDE for bulk with GPD's for both upper and lower tails
pgkg <- function(q, kerncentres, lambda = NULL, 
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE, 
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
  
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")
  
  q[is.infinite(q)]=NA # user will have to deal with infinite cases
  
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
  
  p = q # this will pass through NA/NaN in q just as they are entered
    
  whichul = which(q < ul)
  nul = length(whichul)
  whichb = which((q <= ur) & (q >= ul)) 
  nb = length(whichb)
  whichur = which(q > ur)
  nur = length(whichur)
  
  if (nul > 0) p[whichul] = 1 - pgpd(-q[whichul], -ul, sigmaul, xil, phiul)
  if (nb > 0) p[whichb] = phiul + phib*(sapply(q[whichb], FUN = pkdenx, kerncentres = kerncentres, lambda = lambda) - pkdenx(ul, kerncentres, lambda))
  if (nur > 0) p[whichur] = pgpd(q[whichur], ur, sigmaur, xir, phiur)
  
  if (!lower.tail) p = 1 - p
  
  p
}
