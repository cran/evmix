#' @export
#' @aliases  mgammagpd dmgammagpd pmgammagpd qmgammagpd rmgammagpd
#' @rdname mgammagpd

# cumulative distribution function for mixture of gammas bulk with GPD for upper tail
pmgammagpd <- function(q, mgshape = list(1), mgscale = list(1), mgweights = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
  sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
    
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")

  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")

  q[is.infinite(q)]=NA # user will have to deal with infinite cases

  if (!is.list(mgshape) | !is.list(mgshape))
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")

  allM = c(length(mgshape), length(mgshape))
  if (!is.null(mgweights)) allM = c(allM, length(mgweights))
  M = unique(allM) # number of components
  if (length(M) != 1)
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(q), length(u), length(sigmau), length(xi), length(phiu))
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

  q = rep(q, length.out = n)
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
    
  p = q # this will pass through NA/NaN in q just as they are entered
  
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0){
    pmg = rep(0, nb)
    for (i in 1:M) {
      pmg = pmg + pgamma(q[whichb], shape = mgshape[[i]][whichb],
        scale = mgscale[[i]][whichb])*mgweights[[i]][whichb]
    }
    p[whichb] = phib[whichb] * pmg
  }

  if (nu > 0) p[whichu] = 1 - phiu[whichu] + phiu[whichu] * pgpd(q[whichu], u[whichu], sigmau[whichu], xi[whichu])

  if (!lower.tail) p = 1 - p

  p
}
