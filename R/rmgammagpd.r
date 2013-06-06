#' @export
#' @aliases mgammagpd dmgammagpd pmgammagpd qmgammagpd rmgammagpd
#' @rdname mgammagpd

# random number generation for mixture of gammas bulk with GPD for upper tail
rmgammagpd <- function(n = 1, mgshape = list(1), mgscale = list(1), mgweights = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
  sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0, phiu = TRUE) {

  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")

  if (!is.list(mgshape) | !is.list(mgshape))
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")

  allM = c(length(mgshape), length(mgshape))
  if (!is.null(mgweights)) allM = c(allM, length(mgweights))
  M = unique(allM) # number of components
  if (length(M) != 1)
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(u), length(sigmau), length(xi), length(phiu))
  for (i in 1:M) linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]))
  if (!is.null(mgweights)) {
    for (i in 1:M) linputs = c(linputs, length(mgweights[[i]]))
  }

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")

  if (mode(u) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  for (i in 1:M) {
    if (mode(mgshape[[i]]) != "numeric" | mode(mgscale[[i]]) != "numeric")
      stop(paste("gamma parameters for component", i, "must be numeric"))
  } 
  
  if (any(!is.finite(c(u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  for (i in 1:M) {
    if (any(!is.finite(c(mgshape[[i]], mgscale[[i]]))))
      stop(paste("gamma parameters for component", i, "must be numeric"))
  
    if (min(mgshape[[i]]) < 0)
      stop(paste("gamma shape for component", i, "must be non-negative"))

    if (min(mgscale[[i]]) < 0)
      stop(paste("gamma scale for component", i, "must be non-negative"))
  }

  if (is.null(mgweights)) {
    mgweights = as.list(rep(1/M, M))
  } else {
    for (i in 1:M) {
      if (mode(mgweights[[i]]) != "numeric" | any(!is.finite(c(mgweights[[i]]))))
        stop(paste("gamma weights for component", i, "must be numeric"))

      if (min(mgweights[[i]]) <= 0)
        stop(paste("gamma weights for component", i, "must be non-negative"))
    }
  }
  
  if (min(u) <= 0)
    stop("threshold must be non-negative")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")

  if (any(xi == 1))
    stop("shape cannot be 1")

  if (is.logical(phiu) & any(!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiu < 0) | any(phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }

  qmgammagpd(runif(n), mgshape, mgscale, mgweights, u, sigmau, xi, phiu)
}


