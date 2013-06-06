#' @export
#' @aliases mgammagpd dmgammagpd pmgammagpd qmgammagpd rmgammagpd
#' @rdname mgammagpd

# inverse cumulative distribution function for mixture of gammas bulk with GPD for upper tail
qmgammagpd <- function(p, mgshape = list(1), mgscale = list(1), mgweights = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
  sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
    
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")

  if (!is.list(mgshape) | !is.list(mgshape))
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")

  allM = c(length(mgshape), length(mgshape))
  if (!is.null(mgweights)) allM = c(allM, length(mgweights))
  M = unique(allM) # number of components
  if (length(M) != 1)
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(p), length(u), length(sigmau), length(xi), length(phiu))
  for (i in 1:M) linputs = c(linputs, length(mgshape[[i]]), length(mgscale[[i]]))
  if (!is.null(mgweights)) {
    for (i in 1:M) linputs = c(linputs, length(mgweights[[i]]))
  }
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

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
  
  if (!lower.tail) p = 1 - p

  p = rep(p, length.out = n)
  for (i in 1:M) {
    mgshape[[i]] = rep(mgshape[[i]], length.out = n)
    mgscale[[i]] = rep(mgscale[[i]], length.out = n)
    mgweights[[i]] = rep(mgweights[[i]], length.out = n)
  }
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # Renormalise weights
  totweights = colSums(matrix(unlist(mgweights), nrow = M, ncol = n, byrow = TRUE))
  pmg = rep(0, n)
  for (i in 1:M) {
    mgweights[[i]] = mgweights[[i]]/totweights
    pmg = pmg + pgamma(u, shape = mgshape[[i]], scale = mgscale[[i]])*mgweights[[i]]
  }
  
  if (is.logical(phiu)) {
    phiu = 1 - pmg
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pmg
    
  q = p # this will pass through NA/NaN in p just as they are entered
  
  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)

  # obtain quantile function by interpolation
  # - CDF is quick to calculate, interpolation algorithms also quite fast
  # - an alternative solution is to numerical solve to find quantile, but this is much slower
  # First single parameter valuesthis is much quicker
  if (nb > 0) {
    if (all(linputs[-1] == 1)) {
      qrange = range(0, u[1])
      qk = seq(qrange[1], qrange[2], length.out = 1000)

      pk = sapply(qk, FUN = pmgammagpd, 
        mgshape = lapply(mgshape, FUN = function(x) x[1]),
        mgscale = lapply(mgscale, FUN = function(x) x[1]),
        mweights = lapply(mgweights, FUN = function(x) x[1]),
        u = u[1], sigmau = sigmau[1], xi = xi[1], phiu = phiu[1]) 

      qfun = splinefun(x = pk, y = qk)

      q[whichb] = qfun(p[whichb])
    } else {
      for (i in whichb) {
        qrange = range(0, u[i])
        qk = seq(qrange[1], qrange[2], length.out = 1000)

        pk = phib[i] * sapply(qk, FUN = pmgammagpd, 
          mgshape = lapply(mgshape, FUN = function(x) x[i]),
          mgscale = lapply(mgscale, FUN = function(x) x[i]),
          mweights = lapply(mgweights, FUN = function(x) x[i]),
          u = u[i], sigmau = sigmau[i], xi = xi[i], phiu = phiu[i]) 

        qfun = splinefun(x = pk, y = qk)

        q[i] = qfun(p[i])
      }
    }
  }

  if (nu > 0) q[whichu] = qgpd(p[whichu], u[whichu], sigmau[whichu], xi[whichu], phiu[whichu])

  q
}

