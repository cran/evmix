#' @export
#' @aliases  lkden nlkden
#' @rdname lkden

# negative cross-validation log-likelihood function for kernel density estimator using normal kernel
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlkden <- function(lambda, x, extracentres = NULL, finitelik = FALSE) {

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

  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  nllh = -lkden(x, lambda = lambda, extracentres = extracentres) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
