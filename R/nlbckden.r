#' @export
#' @aliases  lbckden nlbckden
#' @rdname lbckden

# negative cross-validation log-likelihood function for boundary corrected KDE
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlbckden <- function(lambda, x, extracentres = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, finitelik = FALSE) {

  # Check properties of inputs
  if (missing(lambda))
    stop("lambda must be initial guess of bandwith (scalar)")
    
  if (mode(lambda) != "numeric") 
    stop("lambda must be initial guess of bandwith (scalar)")

  if (length(lambda) != 1)
    stop("lambda must be initial guess of bandwith (scalar)")

  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")

  if (!is.null(extracentres)) {
    if (length(extracentres) == 0 | mode(extracentres) != "numeric") 
      stop("extracentres must be a non-empty numeric vector")
    
    if (any(is.infinite(extracentres)))
      stop("infinite cases in extracentres must be removed")
  }

  allmethods = c("simple", "renorm", "reflect", "logtrans", 
    "beta1", "beta2", "gamma1", "gamma2", "copula")
  if (length(bcmethod) != 1)
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if ((mode(bcmethod) != "character"))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))
  
  if (!(bcmethod %in% allmethods))
    stop(paste("bcmethod must be one of", allmethods, collapse = "TRUE" ))

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

  if ((bcmethod %in% upboundmethods) & (max(x) > xmax))
    stop("xmax must be higher than largest kernel centre as it is defined as upper bound")
  
  if (!is.null(extracentres)){
    if ((bcmethod %in% upboundmethods) & (max(extracentres) > xmax))
      stop("xmax must be higher than largest kernel centre as it is defined as upper bound")
  }
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  nllh = -lbckden(x, lambda = lambda, extracentres = extracentres, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax) 

  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
