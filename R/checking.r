#' @name checking
#' 
#' @param x            scalar or vector of quantiles
#' @param prob         scalar or vector of probability  
#' @param n            scalar sample size
#' @param param        scalar or vector of parameters
#' @param phiu         scalar or vector of phiu (logical, NULL or 0-1 exclusive)
#' @param nparam       vector of lengths of parameter vectors
#' @param logicarg     logical input argument
#' @param textarg      character input argument
#' @param inputn       vector of input lengths
#' @param allowvec     logical, where TRUE permits vector
#' @param allownull    logical, where TRUE permits NULL values
#' @param allowmiss    logical, where TRUE permits NA and NaN values
#' @param allowna      logical, where TRUE permits NA values
#' @param allowinf     logical, where TRUE permits Inf values
#' @param allowfalse   logical, where TRUE permits FALSE (and TRUE) values
#' @param method       optimisation method (see \code{\link[stats:optim]{optim}})
#' @param control      optimisation control list (see \code{\link[stats:optim]{optim}})
#' @param bcmethod     boundary correction method
#' @param nn           non-negativity correction method (simple boundary correction only)
#' @param offset       offset added to kernel centres (logtrans only) or \code{NULL}
#' 
#' @title Functions for checking function input argument
#'
#' @description Functions for checking the input arguments to functions, so that main functions
#' are more concise. They will stop when an inappropriate input is found.
#' 
#' For likelihood functions you will often not want to stop on finding a non-positive values for
#' postive parameters, in such cases use \code{\link[evmix:checking]{check.param}} rather than 
#' \code{\link[evmix:checking]{check.posparam}}.
#' 
#' @return The checking functions will stop on errors and return no value. The only exception is
#' the \code{\link[evmix:checking]{check.inputn}} which outputs the maximum vector length.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @aliases checking check.quant check.prob check.n check.param check.posparam
#'  check.nparam check.logic check.text check.inputn
#' @family checking
#' 
NULL

#' @export
#' @rdname checking
check.quant <- function(x, allowvec = TRUE, allownull = FALSE, allowmiss = FALSE, allowinf = FALSE) {
  # vector of quantile accepted by default
  args = match.call()
  
  arg.name = args[[2]]

  if (missing(x)) stop(paste("quantile", arg.name, "is missing"))

  if (!allownull) {
    if (is.null(x)) stop(paste("quantile", arg.name, "is NULL"))  
  }  
  
  if (!is.null(x)) {
    if ((mode(x) != "numeric") & !((length(x) == 1) & is.na(x[1])))
      stop(paste("quantile", arg.name, "must be", ifelse(allowvec, "numeric vector", "scalar")))
        
    if (!allowvec & (length(x) > 1))
      stop(paste("quantile", arg.name, "must be scalar"))

    if (any(!is.finite(x))) {
      if (!allowinf & any(is.infinite(x)))
        stop(paste("infinite values for quantile", arg.name, "not permitted"))
      
      if (!allowmiss & any(is.na(x)))
        stop(paste("missing values for quantile", arg.name, "not permitted"))
    }
    
    if (length(x) == 0)
      stop(paste("quantile", arg.name, "must be scalar or vector"))
  }
}

#' @export
#' @rdname checking
check.prob <- function(prob, allowvec = TRUE, allownull = FALSE, allowmiss = FALSE) {
  # vector of probability accepted by default

  if (missing(prob)) stop("probability is missing")

  if (!allownull) {
    if (is.null(prob)) stop("probability is NULL")  
  }
  
  if (!is.null(prob)) {
    if ((mode(prob) != "numeric") & !((length(prob) == 1) & is.na(prob[1])))
      stop(paste("probability must be a non-empty numeric", ifelse(allowvec, "vector", "scalar")))
        
    if (!allowvec & (length(prob) > 1))
      stop("probability must be scalar")

    if (!allowmiss & any(!is.finite(prob))) stop("missing values for probability not permitted")
  
    if (any(prob < 0, na.rm = TRUE) | any(prob > 1, na.rm = TRUE))
      stop("probability must be between 0 and 1 (inclusive)")
  }
}

#' @export
#' @rdname checking
check.n <- function(n) {
  # positive integer only permitted

  if (missing(n)) stop("sample size is missing")
  
  if (is.null(n)) stop("sample size is NULL")  
  
  if (length(n) != 1 | (mode(n) != "numeric"))
    stop("sample size must be a positive integer")
        
  if (!is.finite(n)) stop("non-finite sample size not permitted")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
}

#' @export
#' @rdname checking
check.param <- function(param, allowvec = FALSE, allownull = FALSE, allowmiss = FALSE, allowinf = FALSE) {
  args = match.call()
  
  arg.name = args[[2]]

  if (missing(param)) stop(paste("parameter input", arg.name, "is missing"))

  if (!allownull) {
    if (is.null(param)) stop(paste("parameter input", arg.name, "is NULL"))
  }
  
  if (!is.null(param)) {
    if ((mode(param) != "numeric") & !((length(param) == 1) & is.na(param[1])))
      stop(paste("parameter input", arg.name, "must be a non-empty numeric",
        ifelse(allowvec, "vector", "scalar")))
  
    if (!allowvec & (length(param) > 1))
      stop(paste("parameter input", arg.name, "must be scalar"))
      
    if (any(!is.finite(param))) {
      if (!allowinf & any(is.infinite(param)))
        stop(paste("infinite cases for parameter input", arg.name, "not permitted"))
      
      if (!allowmiss & any(is.na(param)))
        stop(paste("missing values for parameter input", arg.name, "not permitted"))
    }
  }
}

#' @export
#' @rdname checking
check.posparam <- function(param, allowvec = FALSE, allownull = FALSE, allowmiss = FALSE, allowinf = FALSE) {
  args = match.call()
  
  arg.name = args[[2]]

  check.param(param, allowvec, allownull, allowmiss, allowinf)
  
  if (any(param < 0, na.rm = TRUE)) stop(paste("non-positive parameter", arg.name, "not permitted"))
}

#' @export
#' @rdname checking
check.logic <- function(logicarg, allowvec = FALSE, allowna = FALSE) {
  # only logical, NA never permitted (NaN not permitted as non-logical)
  args = match.call()
  
  arg.name = args[[2]]

  if (missing(logicarg)) stop(paste("logical input", arg.name, "is missing"))

  if ((length(logicarg) == 0) | (mode(logicarg) != "logical"))
    stop(paste("input", arg.name, "must be logical")) 

  if (!allowvec & (length(logicarg) > 1))
    stop(paste("parameter input", arg.name, "must be logical"))
      
  if (!allowna & any(is.na(logicarg)))
    stop(paste("NA values for logical input", arg.name, "not permitted"))
}

#' @export
#' @rdname checking
check.nparam <- function(param, nparam = 1, allownull = FALSE) {
  args = match.call()
  
  arg.name = args[[2]]

  if (missing(param)) stop(paste("parameter input", arg.name, "is missing"))

  if (!allownull) {
    if (is.null(param)) stop(paste("parameter input", arg.name, "is NULL"))
  }
  
  if (!is.null(param)){
    if ((mode(param) != "numeric") & !((length(param) == 1) & is.na(param[1])))
      stop(paste("parameter input", arg.name, "must be a non-empty numeric of length", nparam))
  
    if (length(param) != nparam)
      stop(paste("parameter input", arg.name, "must be a numeric vector of length", nparam))

    if (any(!is.finite(param))) stop(paste("non-finite values for parameter input", arg.name, "not permitted"))
  }
}

#' @export
#' @rdname checking
check.text <- function(textarg, allowvec = FALSE) {
  # character string only
  args = match.call()
  
  arg.name = args[[2]]

  if (missing(textarg)) stop(paste("character input", arg.name, "is missing"))

  if ((length(textarg) != 1) | (mode(textarg) != "character"))
    stop(paste("input", arg.name, "must be character")) 
}

#' @export
#' @rdname checking
check.inputn <- function(inputn, allowvec = TRUE) {
  # parameter lengths only, vector of inputs permitted by default

  if (missing(inputn)) stop("vector of input lengths is missing")
    
  if ((length(inputn) == 0) | (mode(inputn) != "numeric"))
    stop("input lengths must be integer")
          
  if (!allowvec & any(inputn > 1))
    stop("input lengths must be finite integers")

  if (any(!is.finite(inputn)))
    stop("input lengths must be finite integers")

  n = max(inputn)

  if (n == 0)
    stop("input lengths must be either scalar or vector, with all vectors of same length")

  if (any(inputn[inputn != 1] != n))
    stop("input lengths must be either scalar or vector, with all vectors of same length")
  
  return(n)
}

#' @export
#' @rdname checking
check.phiu <- function(phiu, allowvec = FALSE, allownull = FALSE, allowfalse = FALSE) {
  # missing and infinite values for phiu never permitted
  # TRUE or 0-1 (exclusive) in other functions only, as sample proportion (FALSE) not available
  # in addition FALSE permitted in likelihood only, as full sample data is available

  if (missing(phiu)) stop("phiu is missing")

  if (!allownull) {
    if (is.null(phiu)) stop("phiu is NULL")
  }

  if (allowfalse) {
    phiuerror = "phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach"
  } else{
    phiuerror = "phiu must be either TRUE for bulk model based tail fraction approach, 
      FALSE for parameterised tail fraction approach (estimated using sample proportion)
      or between 0 and 1 (exclusive) for fixed tail fraction"
  }
  
  if (!is.null(phiu)) {
    if (!allowvec & (length(phiu) > 1))
      stop("phiu must be of length one")

    if (any(is.na(phiu))) stop("missing values for phiu not permitted")

    if (is.logical(phiu)) {
      if (!allowfalse) {
        if (any(!phiu)) stop(phiuerror)
      }
    } else {
      if (mode(phiu) != "numeric") stop(phiuerror)
      
      if (any(phiu < 0, na.rm = TRUE) | any(phiu > 1, na.rm = TRUE))
        stop("phiu must between 0 and 1 (inclusive)")
    }

    if (any(!is.finite(phiu))) stop("infinite and missing values for phiu not permitted")
  }
}

#' @export
#' @rdname checking
check.optim <- function(method) {
  # checks optim method specification

  if (missing(method)) stop("optim method is missing")

  alloptim = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")
  if ((mode(method) != "character") | (length(method) != 1))
    stop(paste(c("optim method must be one of", alloptim), collapse = " "))
  
  if (!(method %in% alloptim))
    stop(paste(c("optim method must be one of", alloptim), collapse = " "))
}

#' @export
#' @rdname checking
check.control <- function(control) {
  # checks optim method specification

  if (missing(control)) stop("optim control list is missing")

  if (mode(control) != "list")
    stop("optim control must be list, see optim function for details")

  allcontrol = c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol",
    "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "tmax", "temp")  
  if (any(!(names(control) %in% allcontrol)))
    stop(paste(c("optim control list must be containly only", allcontrol), collapse = " "))
}

#' @export
#' @rdname checking
check.bcmethod <- function(bcmethod) {
  # checks boundary correction method specification

  if (missing(bcmethod)) stop("boundary correction method is missing")

  allbcmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans", 
    "beta1", "beta2", "gamma1", "gamma2", "copula")
  if ((mode(bcmethod) != "character") | (length(bcmethod) != 1))
    stop(paste(c("boundary correction method must be one of", allbcmethods), collapse = " "))
  
  if (!(bcmethod %in% allbcmethods))
    stop(paste(c("boundary correction method must be one of", allbcmethods), collapse = " "))
}

#' @export
#' @rdname checking
check.nn <- function(nn) {
  # checks non-negative correction method specification

  if (missing(nn)) stop("non-negative correction method is missing")

  allnn = c("none", "zero", "jf96")
  if ((mode(nn) != "character") | (length(nn) != 1))
    stop(paste(c("non-negative correction method must be one of", allnn), collapse = " "))
  
  if (!(nn %in% allnn))
    stop(paste(c("non-negative correction method must be one of", allnn), collapse = " "))
}

#' @export
#' @rdname checking
check.offset <- function(offset, bcmethod) {
  # checks offset for different bcmethods

  if (missing(offset)) stop("offset is missing")

  if (missing(bcmethod)) stop("bcmethod is missing")

  if (bcmethod != "logtrans") {
    if (!is.null(offset)) warning("offset only relevant for logtrans method")
  } else {
    if (is.null(offset)) stop("offset must be provided for logtrans method")
    
    if (offset <= 0) stop("offset must be positive")
  }
}


