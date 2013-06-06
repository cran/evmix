#' @name gammagpd
#' 
#' @title Gamma Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with gamma for bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the gamma shape \code{gshape} and scale \code{gscale}, threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams weibullgpd
#' @inheritParams normgpd
#' @inheritParams gpd
#' @param gshape     gamma shape (non-negative)
#' @param gscale     gamma scale (non-negative)
#' 
#' @details Extreme value mixture model combining gamma distribution for the bulk
#' below the threshold and GPD for upper tail. The user can pre-specify \code{phiu} 
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
#' where \eqn{H(x)} and \eqn{G(X)} are the gamma and conditional GPD
#' cumulative distribution functions (i.e. \code{pgamma(x, gshape, scale = gscale)} and
#' \code{pgpd(x, u, sigmau, xi)}).
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
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[stats:GammaDist]{dgamma}} for details of gamma bulk component.
#' 
#' @return \code{\link[evmix:gammagpd]{dgammagpd}} gives the density, 
#' \code{\link[evmix:gammagpd]{pgammagpd}} gives the cumulative distribution function,
#' \code{\link[evmix:gammagpd]{qgammagpd}} gives the quantile function and 
#' \code{\link[evmix:gammagpd]{rgammagpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{rgammagpd} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:gammagpd]{rgammagpd}} is 1.
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
#' Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004). Bayesian analysis of extreme
#' events with threshold estimation. Statistical Modelling. 4(3), 227-244.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:GammaDist]{dgamma}}
#' @aliases  gammagpd dgammagpd pgammagpd qgammagpd rgammagpd
#' @family   gammagpd
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2))
#' x = rgammagpd(1000, gshape = 2)
#' xx = seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgammagpd(xx, gshape = 2))
#' 
#' # three tail behaviours
#' plot(xx, pgammagpd(xx, gshape = 2), type = "l")
#' lines(xx, pgammagpd(xx, gshape = 2, xi = 0.3), col = "red")
#' lines(xx, pgammagpd(xx, gshape = 2, xi = -0.3), col = "blue")
#' legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rgammagpd(1000, gshape = 2, u = 3, phiu = 0.2)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgammagpd(xx, gshape = 2, u = 3, phiu = 0.2))
#' 
#' plot(xx, dgammagpd(xx, gshape = 2, u = 3, xi=0, phiu = 0.2), type = "l")
#' lines(xx, dgammagpd(xx, gshape = 2, u = 3, xi=-0.2, phiu = 0.2), col = "red")
#' lines(xx, dgammagpd(xx, gshape = 2, u = 3, xi=0.2, phiu = 0.2), col = "blue")
#' legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   }
NULL

#' @export
#' @aliases  gammagpd dgammagpd pgammagpd qgammagpd rgammagpd
#' @rdname gammagpd

# probability density function for gamma bulk with GPD for upper tail
dgammagpd <- function(x, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 1/gscale), 
  sigmau = sqrt(gshape) * gscale, xi = 0, phiu = TRUE, log = FALSE) {

  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    warning("infinite cases are set to NaN")

  x[is.infinite(x)]=NA # user will have to deal with infinite cases

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(x), length(gshape), length(gscale), 
    length(u), length(sigmau), length(xi), length(phiu))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

  if (mode(gshape) != "numeric" | mode(gscale) != "numeric" | mode(u) != "numeric" |
    mode(sigmau) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(gshape, gscale, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (min(gshape) < 0)
    stop("gamma shape must be non-negative")

  if (min(gscale) < 0)
    stop("gamma scale must be non-negative")
  
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
  gshape = rep(gshape, length.out = n)
  gscale = rep(gscale, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  if (is.logical(phiu)) {
    phiu = 1 - pgamma(u, shape = gshape, scale = gscale)
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pgamma(u, shape = gshape, scale = gscale)

  d = x # this will pass through NA/NaN in x just as they are entered
  
  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)
  
  if (nb > 0) d[whichb] = log(phib[whichb]) + dgamma(x[whichb], shape = gshape[whichb], scale = gscale[whichb], log = TRUE)
  if (nu > 0) d[whichu] = log(phiu[whichu]) + dgpd(x[whichu], u[whichu], sigmau[whichu], xi[whichu], log = TRUE)

  if (!log) d = exp(d)

  d
}
