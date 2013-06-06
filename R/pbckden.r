
#' @export
#' @aliases bckden dbckden pbckden qbckden
#' @rdname  bckden

# cumulative distribution function for boundary corrected KDE
pbckden <- function(q, kerncentres, lambda = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
    
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite cases are set to NaN")

  q[is.infinite(q)] = NaN # user will have to deal with infinite cases

  if (missing(kerncentres))
    stop("kerncentres must be a non-empty numeric vector")
    
  if (length(kerncentres) == 0 | mode(kerncentres) != "numeric") 
    stop("kerncentres must be a non-empty numeric vector")
  
  if (any(!is.finite(kerncentres)))
    warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]
  nk = length(kerncentres)

  if (any(kerncentres < 0))
    stop("kernel centres cannot be non-positive")

  if (is.unsorted(kerncentres)) kerncentres = sort(kerncentres)
  
  if (is.null(lambda))
    stop("bandwidth (lambda) must be specified")

  if (length(lambda) != 1)
    stop("bandwidth must be scalar")

  if ((mode(lambda) != "numeric") & !is.finite(lambda))
    stop("bandwidth must be numeric")
  
  if (lambda <= 0)
      stop("bandwidth must be non-negative")  

  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  allmethods = c("simple", "renorm", "reflect", "logtrans", 
    "beta1", "beta2", "gamma1", "gamma2", "copula")
  if (length(bcmethod) != 1)
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if ((mode(bcmethod) != "character"))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))
  
  if (!(bcmethod %in% allmethods))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if ((bcmethod == "copula") & (lambda >= 1))
      stop("bandwidth must between (0, 1) for copula method")  

  if (!is.logical(proper))
    stop("proper must be logical")
  
  if (length(proper) != 1)
    stop("proper must be of length 1")

  allnn = c("none", "zero", "jf96")
  if (length(nn) != 1)
    stop(paste("nn must be one of", allnn, collapse = " "))

  if ((mode(nn) != "character"))
    stop(paste("nn must be one of", allnn, collapse = " "))
  
  if (!(nn %in% allnn))
    stop(paste("nn must be one of", allnn, collapse = " "))
  
  if (length(offset) != 1)
    stop("offset must be scalar")

  if ((mode(offset) != "numeric") & !is.finite(offset))
    stop("offset must be numeric")
  
  if (offset < 0)
      stop("offset must be non-negative")  

  if ((offset != 0) & (bcmethod != "logtrans"))
    warning("offset only relevant for logtrans method")
    
  if (length(xmax) != 1)
    stop("xmax must be scalar")

  if ((mode(xmax) != "numeric"))
    stop("xmax must be numeric")
  
  if (xmax <= 0)
      stop("xmax must be positive")  

  upboundmethods = c("beta1", "beta2", "copula")
  if ((!is.infinite(xmax)) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for methods", upboundmethods, collapse = " "))
    
  if ((is.infinite(xmax)) & (bcmethod %in% upboundmethods))
    warning(paste("xmax cannot be infinite for methods", upboundmethods, collapse = " "))

  if ((bcmethod %in% upboundmethods) & (max(kerncentres) > xmax))
    stop("xmax must be higher than largest kernel centre as it is defined as upper bound")
   
  # problems caused in evaluation if no data near boundary
  # so don't evaluate cdf if no kernel centres near threshold
  if (bcmethod %in% upboundmethods) {
    minaccept = 0
    maxaccept = xmax
  } else if (bcmethod == "logtrans") {
    maxaccept = exp(log(max(kerncentres) + offset) + 5*lambda)
    minaccept = max(0, exp(log(min(kerncentres) + offset) - 5*lambda))
  } else {
    maxaccept = max(kerncentres) + 5*lambda
    minaccept = max(0, min(kerncentres) - 5*lambda)
  }
  
  p = q # this will pass through NA/NaN in x just as they are entered

  p[which(!is.na(q))] = ifelse(q[which(!is.na(q))] <= minaccept, 0, p)
  p[which(!is.na(q))] = ifelse(q[which(!is.na(q))] > maxaccept, 1, p) 
  # don't include upper limit so user can get integral upto maxaccept if not normalised
  
  # only select non-negative x-values for evaluation
  qok = q[!is.na(q)]
  qok = qok[(qok >= minaccept) & (qok <= maxaccept)]

  if (length(qok) > 0) {
    if (bcmethod == "simple") {
      # simple linear boundary correction method of Jones (1993), eq 3.4
      # which adapts normal kernels so they are sensibly behaved near boundary
      # use similar notation to original paper
      # (truncpoint is their p and lambda is their S_k)
      
      pok = sapply(qok, FUN = simplepbckdenx, kerncentres = kerncentres, lambda = lambda)

    } else if (bcmethod == "renorm") {
      # really crude renormalisation, based on how much of kernel is above boundary
      
      pok = sapply(qok, FUN = renormpbckdenx, kerncentres = kerncentres, lambda = lambda)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "reflect") {
      # really crude reflection, equivalent to using (kerncentres, -kerncentres)
      # only good if f'(0)=0
      
      pok = sapply(qok, FUN = reflectpbckdenx, kerncentres = kerncentres, lambda = lambda)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "logtrans") {
      # simple log transformation 
      # (Wand, Marron and Ruppert JASA, 1991, Marron and Ruppert JRSS B, 1994)
      # requires offset in min(x)=0
      if ((min(kerncentres) == 0) & (offset = 0)){
        stop("log transformation requires offset for zero kernel centre")
      }
      
      # transformation KDE is f(x) = g(t(x)) abs(t'(x))
      pok = sapply(qok, FUN = logpbckdenx, kerncentres = kerncentres, lambda = lambda, offset = offset)

      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta1") {
      # beta kernels by Chen (1999) CSDA
      
      pok = sapply(qok, FUN = beta1pbckdenx, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta2") {
      # modified beta kernels by Chen (1999) CSDA
      
      pok = sapply(qok, FUN = beta2pbckdenx, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "gamma1") {
      # gamma kernels by Chen (2000) AISM
      
      dok = sapply(qok, FUN = gamma1pbckdenx, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    } else if (bcmethod == "gamma2") {
      # modified gamma kernels by Chen (2000) AISM
      
      dok = sapply(qok, FUN = gamma2pbckdenx, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    } else if (bcmethod == "copula") {
      # bivariate normal copula based kernels by Jones and Henderson (2007) Biometrika
      
      pok = sapply(qok, FUN = copulapbckdenx, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      proper = FALSE # not needed
      nn = "none"    # not needed
    }

    # Apply non-negative correction first followed by renormalisation (if either requested)
    if (nn == "jf96") {
      pbar = sapply(qok, FUN = nnpbckdenx, kerncentres = kerncentres, lambda = lambda)
      pok = ifelse(pbar == 0, pok, pbar*exp(pok/pbar - 1))
    } else if (nn =="zero") {
      pok[pok < 0] = 0
    }
    
    if (proper) {
      if (bcmethod == "simple") {
        pxmax = simplepbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda)
      } else if (bcmethod == "beta1") {
        pxmax = beta1pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "beta2") {
        pxmax = beta2pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "gamma1") {
        pxmax = gamma1pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda)
      } else if (bcmethod == "gamma2") {
        pxmax = gamma2pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda)
      }
      pok = pok/pxmax
    }

    pok = ifelse(pok > 1, 1, pok) # numerical errors occassionally give cdf > 1
    
    p[ifelse(!is.na(q), (q >= minaccept) & (q <= maxaccept), FALSE)] = pok
  }
  
  if (!lower.tail) p = 1 - p
  p
}
