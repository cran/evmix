#' @export
#' @aliases  gammagpdcon dgammagpdcon pgammagpdcon qgammagpdcon rgammagpdcon
#' @rdname gammagpdcon

# cumulative distribution function for gamma bulk with GPD for upper tail
pgammagpdcon <- function(q, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 1/gscale),
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
  linputs = c(length(q), length(gshape), length(gscale), 
              length(u), length(xi), length(phiu))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
  if (mode(gshape) != "numeric" | mode(gscale) != "numeric" | mode(u) != "numeric" |
    mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(gshape, gscale, u, xi, phiu))))
    stop("parameters must be numeric")
  
  if (min(gshape) < 0)
    stop("gamma shape must be non-negative")
  
  if (min(gscale) < 0)
    stop("gamma scale must be non-negative")
  
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
  gshape = rep(gshape, length.out = n)
  gscale = rep(gscale, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  if (is.logical(phiu)) {
    phiu = 1 - pgamma(u, shape = gshape, scale = gscale)
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pgamma(u, shape = gshape, scale = gscale)
  
  sigmau = phiu / (phib * dgamma(u, shape = gshape, scale = gscale))
  
  if (any(!is.finite(sigmau)))
    stop("sigmau is not numeric")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  pgammagpd(q, gshape, gscale, u, sigmau, xi, phiu, lower.tail)
    
}

