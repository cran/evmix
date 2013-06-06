#' @name mgammagpd
#' 
#' @title Mixture of Gammas Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with mixture of gammas for bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the gamma shape \code{gshape} and scale \code{gscale}, threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams gammagpd
#' @param mgshape     mgamma shape (non-negative) as list
#' @param mgscale     mgamma scale (non-negative) as list
#' @param mgweights   mgamma weights (positive) as list or \code{NULL}
#' 
#' @details Extreme value mixture model combining mixture of gammas for the bulk
#' below the threshold and GPD for upper tail. The parameters are input as a list,
#' with one parameter object in the list for each gamma component. There must be the same number of 
#' components in \code{mgshape} and \code{mgscale}. The number of objects in the parameters lists
#' determines the number of components. The parameter object for each gamma component
#' can either be a scalar or vector, consistent with the other mixture models
#' 
#' If \code{mgweights=NULL} then assumes equal weights for each component. Otherwise, 
#' \code{mgweights} must be a list of the same length as \code{mgshape} and 
#' \code{mgscale}, filled with positive values. In the latter case, the weights are rescaled
#' to sum to unity.
#' 
#' The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' gamma bulk model.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the gamma bulk model (\code{phiu=TRUE}), upto the 
#' threshold \eqn{0 < x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the mixture of gammas and conditional GPD
#' cumulative distribution functions respectively.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{0 < x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' The gamma is defined on the non-negative reals, so the threshold must be non-negative.
#' 
#' See \code{\link[evmix:gammagpd]{gammagpd}} for details of simpler parametric mixture model
#' with single gamma for bulk component and GPD for upper tail.
#' 
#' @return \code{\link[evmix:mgammagpd]{dmgammagpd}} gives the density, 
#' \code{\link[evmix:mgammagpd]{pmgammagpd}} gives the cumulative distribution function,
#' \code{\link[evmix:mgammagpd]{qmgammagpd}} gives the quantile function and 
#' \code{\link[evmix:mgammagpd]{rmgammagpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}, and the parameters can 
#' be vectorised within the list.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{rmgammagpd} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:mgammagpd]{rmgammagpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Gamma_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' do Nascimento, F.F., Gamerman, D. and Lopes, H.F. (2011). A semiparametric
#' Bayesian approach to extreme value estimation. Statistical Computing, 22(2), 661-675.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gammagpd]{gammagpd}}, \code{\link[evmix:mgammagpd]{mgammagpd}}, 
#' \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:GammaDist]{dgamma}}
#' @aliases  mgammagpd dmgammagpd pmgammagpd qmgammagpd rmgammagpd
#' @family   mgammagpd
#'   
NULL

#' @export
#' @aliases  mgammagpd dmgammagpd pmgammagpd qmgammagpd rmgammagpd
#' @rdname mgammagpd

# probability density function for mixture of gammas bulk with GPD for upper tail
dmgammagpd <- function(x, mgshape = list(1), mgscale = list(1),  mgweights = NULL,
  u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]), 
  sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0, phiu = TRUE, log = FALSE) {

  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    warning("infinite cases are set to NaN")

  x[is.infinite(x)]=NA # user will have to deal with infinite cases

  if (!is.list(mgshape) | !is.list(mgshape))
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")

  allM = c(length(mgshape), length(mgshape))
  if (!is.null(mgweights)) allM = c(allM, length(mgweights))
  M = unique(allM) # number of components
  if (length(M) != 1)
    stop("gamma mixture parameters must be lists of same length (i.e. one object in list per component)")

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(x), length(u), length(sigmau), length(xi), length(phiu))
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

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  x = rep(x, length.out = n)
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

  d = x # this will pass through NA/NaN in x just as they are entered
  
  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)
  
  if (nb > 0){
    dmg = rep(0, nb)
    for (i in 1:M) {
      dmg = dmg + dgamma(x[whichb], shape = mgshape[[i]][whichb], 
        scale = mgscale[[i]][whichb])*mgweights[[i]][whichb]
    }
    dmg = log(dmg)
    
    d[whichb] = log(phib[whichb]) + dmg
  }
  
  if (nu > 0) d[whichu] = log(phiu[whichu]) + dgpd(x[whichu], u[whichu], sigmau[whichu], xi[whichu], log = TRUE)

  if (!log) d = exp(d)

  d
}
