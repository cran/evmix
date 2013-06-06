#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname bckdengpd

# random number generation for boundary corrected kernel density estimators for the bulk
# distribution upto the threshold and conditional GPD above threshold.
rbckdengpd <- function(n = 1, kerncentres, lambda = NULL,
  u = as.vector(quantile(kerncentres, 0.9)), sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, 
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = 0, xmax = Inf) {
  
  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must non-negative integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must non-negative integer")
  
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
  
  if (is.null(lambda))
    stop("bandwidth (lambda) must be specified")
  
  linputs = c(length(lambda), length(u), length(sigmau), length(xi), length(phiu))
  
  if (sum(linputs != 1) > 0)
    stop("parameters must be scalar")
  
  if (mode(lambda) != "numeric" | mode(u) != "numeric" | mode(sigmau) != "numeric" |
      mode(xi) != "numeric")
    stop("parameters must be numeric")

  if (any(!is.finite(c(lambda, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")  
  
  if (lambda <= 0)
    stop("bandwidth must be non-negative")  
  
  if (sigmau <= 0)
    stop("scale must be non-negative")
  
  if (is.logical(phiu) & (!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
         or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if ((phiu < 0) | (phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }
  
  if (is.logical(phiu)) {
    phiu = 1 - pbckden(u, kerncentres, lambda = lambda, bcmethod = bcmethod,
      proper = proper, nn = nn, offset = offset, xmax = xmax)
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pbckden(u, kerncentres, lambda = lambda, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax)
  
  qbckdengpd(runif(n), kerncentres, lambda, u, sigmau, xi, phiu, bcmethod, proper, nn, offset, xmax)
}
