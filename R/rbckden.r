#' @export
#' @aliases bckden dbckden pbckden qbckden
#' @rdname  bckden

# random number generation for boundary corrected KDE
rbckden <- function(n, kerncentres, lambda = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf) {

  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")

  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
  
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

  qbckden(runif(n), kerncentres = kerncentres, lambda = lambda,
    bcmethod = bcmethod, proper = proper, nn = nn, xmax = xmax, offset = offset)
}
