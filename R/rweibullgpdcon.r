#' @export
#' @aliases  weibullgpdcon dweibullgpdcon pweibullgpdcon qweibullgpdcon rweibullgpdcon
#' @rdname weibullgpdcon

# random number generation for weibull bulk with GPD for upper tail
rweibullgpdcon <- function(n = 1, wshape = 1, wscale = 1, u = qweibull(0.9, wshape, wscale),
  xi = 0, phiu = TRUE) {
  
  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(wshape), length(wscale), 
              length(u), length(xi), length(phiu))
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")
  
  if (mode(wshape) != "numeric" | mode(wscale) != "numeric" | mode(u) != "numeric" |
    mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(wshape, wscale, u, xi, phiu))))
    stop("parameters must be numeric")
  
  if (min(wshape) < 0)
    stop("weibull shape must be non-negative")
  
  if (min(wscale) < 0)
    stop("weibull scale must be non-negative")
  
  if (min(u) <= 0)
    stop("threshold must be non-negative")
  
  if (any(xi == 1))
    stop("shape cannot be 1")
  
  if (is.logical(phiu) & any(!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiu < 0) | any(phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }
  
  qweibullgpdcon(runif(n), wshape, wscale, u, xi, phiu)
}
