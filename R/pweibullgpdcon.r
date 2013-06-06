#' @export
#' @aliases  weibullgpdcon dweibullgpdcon pweibullgpdcon qweibullgpdcon rweibullgpdcon
#' @rdname weibullgpdcon

# cumulative distribution function for weibull bulk with GPD for upper tail
pweibullgpdcon <- function(q, wshape = 1, wscale = 1, u = qweibull(0.9, wshape, wscale),
  xi = 0, phiu = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
  
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")
  
  q[is.infinite(q)]=NA # user will have to deal with infinite cases
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(q), length(wshape), length(wscale), 
              length(u), length(xi), length(phiu))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
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
  
  if (is.logical(phiu) & any(!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiu < 0) | any(phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }
  
  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  q = rep(q, length.out = n)
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  if (is.logical(phiu)) {
    phiu = 1 - pweibull(u, shape = wshape, scale = wscale)
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pweibull(u, shape = wshape, scale = wscale)
    
  sigmau = phiu / (phib * dweibull(u, shape = wshape, scale = wscale))
  
  if (any(!is.finite(sigmau)))
    stop("sigmau is not numeric")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  pweibullgpd(q, wshape, wscale, u, sigmau, xi, phiu, lower.tail)
}

