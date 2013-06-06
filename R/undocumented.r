#' @name internal
#' 
#' @param x            quantile  
#' @param kerncentres	 kernel centres (typically sample data)
#' @param lambda       bandwidth for KDE
#' @param xmax         upper bound on support, for copula and beta based KDE's only
#' @param offset       offset added to kernel centres, log transform based KDE
#' 
#' @title Internal Functions
#'
#' @description Internal functions not designed to be used directly, but are all exported to make them visible to users.
#'
#' @details Internal functions not designed to be used directly. No error
#' checking of the inputs is carried out, so user must be know what they are doing.
#' They are undocumented, but are made visible to the user.
#' 
#' Mostly, these are used in the kernel density estimation functions.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}. Based on code
#' by Anna MacDonald produced for MATLAB.
#'
#' @seealso \code{\link[evmix:kden]{kden}} and \code{\link[evmix:bckden]{bckden}}.
#' 
NULL

# The following functions are supposed to be used internally for functions
# and are not designed to be used directly by users, there is no error checking


#' @export
#' @rdname internal
kdenx <- function(x, kerncentres, lambda) {
  mean(dnorm(x, kerncentres, lambda))
}

#' @export
#' @rdname internal
pkdenx <- function(x, kerncentres, lambda) {
  mean(pnorm(x, kerncentres, lambda))
}

#' @export
#' @rdname internal
simplebckdenx <- function(x, kerncentres, lambda) {
  # distance of kerncentres to lower boundary (p in paper)
  truncpoint = x/lambda
  
  # Use notation of Jones (1993) to make easier to check
  a0 = pnorm(truncpoint)
  a1 = -dnorm(truncpoint)  
  a2 = a0 + truncpoint * a1
  
  # weights in local linear fitting
  denom = (a2*a0 - a1^2)
  lx = a2/denom
  mx = a1/denom
  
  u = (x - kerncentres)/lambda
  
  mean((lx - mx*u)*dnorm(u))/lambda
}

#' @export
#' @rdname internal
simplepbckdenx <- function(x, kerncentres, lambda) {
  # distance of kerncentres to lower boundary (p in paper)
  truncpoint = x/lambda
  
  # Use notation of Jones (1993) to make easier to check
  a0 = pnorm(truncpoint)
  a1 = -dnorm(truncpoint)  
  a2 = a0 + truncpoint * a1
  
  # weights in local linear fitting
  denom = (a2*a0 - a1^2)
  lx = a2/denom
  mx = a1/denom
  
  u = (x - kerncentres)/lambda
  p = 1 + mean(lx*(pnorm(u) - pnorm(truncpoint)) - 
      mx*(dnorm(truncpoint) - dnorm(u)))
}

#' @export
#' @rdname internal
renormbckdenx <- function(x, kerncentres, lambda) {
  # distance of kernel centres to lower boundary
  truncpoint = kerncentres/lambda
  
  # how much of kernel is in range of support, so (1-a0) gives leakage past boundary
  a0 = pnorm(truncpoint)
  
  u = (x - kerncentres)/lambda
  
  mean(dnorm(u)/a0)/lambda
}

#' @export
#' @rdname internal
renormpbckdenx <- function(x, kerncentres, lambda) {
  # distance of kernel centres to lower boundary
  truncpoint = kerncentres/lambda
  
  # how much of kernel is in range of support, so (1-a0) gives leakage past boundary
  a0 = pnorm(truncpoint)
  
  u = (x - kerncentres)/lambda
  
  mean((pnorm(u) - pnorm(0, kerncentres, lambda))/a0)
}

#' @export
#' @rdname internal
reflectbckdenx <- function(x, kerncentres, lambda) {
  mean(dnorm(x, kerncentres, lambda) + dnorm(x, -kerncentres, lambda))
}

#' @export
#' @rdname internal
reflectpbckdenx <- function(x, kerncentres, lambda) {
  mean(pnorm(x, kerncentres, lambda) - pnorm(0, kerncentres, lambda) 
    + pnorm(0, kerncentres, lambda) - pnorm(-x, kerncentres, lambda) )
}

#' @export
#' @rdname internal
pxb <- function(x, lambda) 2*lambda^2 + 2.5 - sqrt(4*lambda^4 + 6*lambda^2 + 2.25 - x^2 - x/lambda)

#' @export
#' @rdname internal
beta1bckdenx <- function(x, kerncentres, lambda, xmax) {
  
  # rescale x and kernel centres, then treat as bounded on [0, 1] so beta kernel used directly
  x = x/xmax
  kerncentres = kerncentres/xmax
  lambda = lambda/xmax
  
  if ((x >= 0) & (x <= 1)){
    d = mean(dbeta(kerncentres, x/lambda + 1, (1 - x)/lambda + 1))
  } else {
    d = 0
  }
  d/xmax
}

#' @export
#' @rdname internal
beta1pbckdenx <- function(x, kerncentres, lambda, xmax) {
  if (x <= 0) {
    p = 0
  } else if (x > xmax) {
    p = 1
  } else {
    # Re-use density function to do numerical integration
    bckdenint = try(integrate(Vectorize(beta1bckdenx, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, xmax = xmax,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of beta based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
beta2bckdenx <- function(x, kerncentres, lambda, xmax) {
  
  # rescale x and kernel centres, then treat as bounded on [0, 1] so beta kernel used directly
  x = x/xmax
  kerncentres = kerncentres/xmax
  lambda = lambda/xmax

  # split x into cases: middle, near 0, near 1 or outside range
  if ((x >= 2*lambda) & (x <= (1 - 2*lambda))) {
    d = mean(dbeta(kerncentres, x/lambda, (1 - x)/lambda))
  } else if ((x >= 0) & (x < 2*lambda)) {
    d = mean(dbeta(kerncentres, pxb(x, lambda), (1 - x)/lambda))
  } else if ((x > (1 - 2*lambda)) & (x <= 1)) {
    d = mean(dbeta(kerncentres, x/lambda, pxb(1 - x, lambda)))
  } else {
    d = 0
  }
  d/xmax
}

#' @export
#' @rdname internal
beta2pbckdenx <- function(x, kerncentres, lambda, xmax) {
  if (x <= 0) {
    p = 0
  } else if (x > xmax) {
    p = 1
  } else {
    # Re-use density function to do numerical integration
    bckdenint = try(integrate(Vectorize(beta2bckdenx, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, xmax = xmax,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of beta based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
gamma1bckdenx <- function(x, kerncentres, lambda) {
  
  if (x >= 0) {
    d = mean(dgamma(kerncentres, shape = x/lambda + 1, scale = lambda))
  } else {
    d = 0
  }
  d
}

#' @export
#' @rdname internal
gamma1pbckdenx <- function(x, kerncentres, lambda) {

  if (x < 0) {
    p = 0
  } else {
    # Re-use density function to do numerical integration
    bckdenint = try(integrate(Vectorize(gamma1bckdenx, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
  
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of gamma based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
gamma2bckdenx <- function(x, kerncentres, lambda) {

  if ((x >= 0) & (x < 2*lambda)) {
    d = mean(dgamma(kerncentres, shape = (x/lambda)^2/4 + 1, scale = lambda))
  } else if (x >= 2*lambda) {
    d = mean(dgamma(kerncentres, shape = x/lambda, scale = lambda))
  } else {
    d = 0
  }
  d
}

#' @export
#' @rdname internal
gamma2pbckdenx <- function(x, kerncentres, lambda) {

  if (x < 0) {
    p = 0
  } else {
    # Re-use density function to do numerical integration
    bckdenint = try(integrate(Vectorize(gamma2bckdenx, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
  
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of gamma based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
copulabckdenx <- function(x, kerncentres, lambda, xmax) {
  # rescale x and kernel centres, then treat as bounded on [0, 1] so beta kernel used directly
  x = x/xmax
  kerncentres = kerncentres/xmax
  
  # quantities used in Gaussian copula kernels
  lambda2 = lambda^2
  stdinv = ifelse((x >= 0) & (x <= 1), qnorm(x), NA)
  kerninv = ifelse((kerncentres >= 0) & (kerncentres <= 1), qnorm(kerncentres), NA)
  
  # Renormalisation constant and (unscaled) density
  d1 = exp(-((1 - lambda2)*stdinv)^2/2/lambda2/(2 - lambda2))/lambda/sqrt(2 - lambda2)
  d2 = mean(exp(-(1 - lambda2)*((1 - lambda2)*kerninv^2 - 2*stdinv*kerninv)/
      2/lambda2/(2 - lambda2)))
  
  # treat those outside range [0, xmax] as 0 density
  d = ifelse((x >= 0) & (x <= 1), ifelse(is.infinite(d2), 0, d1*d2/xmax), 0)
}

#' @export
#' @rdname internal
copulapbckdenx <- function(x, kerncentres, lambda, xmax) {
  if (x <= 0) {
    p = 0
  } else if (x > xmax) {
    p = 1
  } else {
    # Re-use density function to do numerical integration
    bckdenint = try(integrate(Vectorize(copulabckdenx, "x"), lower = 0, upper = x,
      kerncentres = kerncentres, lambda = lambda, xmax = xmax,
      subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
    
    if (inherits(bckdenint, "try-error")) {
      bckdenint$value = NA
      warning("failed to numerically evaluate cdf of copula based boundary corrected KDE")
    }
    p = bckdenint$value
  }
  p
}

#' @export
#' @rdname internal
logpbckdenx <- function(x, kerncentres, lambda, offset) {

  # Re-use density function to do numerical integration
  bckdenint = try(integrate(dbckden, lower = 0, upper = x,
    kerncentres = kerncentres, lambda = lambda, bcmethod = "logtrans",
    offset = offset, subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
  
  if (inherits(bckdenint, "try-error")) {
    bckdenint$value = NA
    warning("failed to numerically evaluate cdf of log transform based boundary corrected KDE")
  }
  bckdenint$value
}

#' @export
#' @rdname internal
nnbckdenx <- function(x, kerncentres, lambda) {
  mean(dnorm(x, kerncentres, lambda))/pnorm(x/lambda)
}

#' @export
#' @rdname internal
nnpbckdenx <- function(x, kerncentres, lambda) {
  mean(pnorm(x, kerncentres, lambda))/pnorm(x/lambda)
}
