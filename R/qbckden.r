#' @export
#' @aliases bckden dbckden pbckden qbckden
#' @rdname  bckden

# inverse cumulative distribution function for boundary corrected KDE
qbckden <- function(p, kerncentres, lambda = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
    
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (any(is.infinite(p)))
    warning("infinite cases are set to NaN")

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
  
  if (!lower.tail) p = 1 - p

  q = p
  
  pok = p[which(!is.na(p))]
  
  qok = rep(NA, length(pok))

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
  
  qk = seq(minaccept, maxaccept, length.out = min(50 * nk, 1000))

  pk = sapply(qk, FUN = pbckden, kerncentres = kerncentres, lambda = lambda,
    bcmethod = bcmethod, proper = proper, nn = nn, xmax = xmax)

  qfun = splinefun(x = pk, y = qk, method = "monoH.FC") # good method for cdf as monotone

  whichinterp = which((pok > min(pk)) & (pok < max(pk)))
  ninterp = length(whichinterp)
  which0 = which(pok <= min(pk))
  n0 = length(which0)
  which1 = which(pok >= max(pk))
  n1 = length(which1)

  if (ninterp > 0) qok[whichinterp] = qfun(pok[whichinterp])

  # if further out then set:
  # lower tail to -Inf and upper tail set to Inf
  if (n0 > 0) qok[which0] = 0
  if (n1 > 0) qok[which1] = rep(ifelse(bcmethod %in% upboundmethods, xmax, Inf), length(which1))
  
  q[which(!is.na(p))] = qok

  q
}
