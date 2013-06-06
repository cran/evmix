#' @name hpdcon
#' 
#' @title Hybrid Pareto Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the hybrid Pareto extreme value mixture model with a 
#'   single continuiuty constraint. The parameters are the normal mean \code{nmean}
#'   and standard deviation \code{nsd}, threshold \code{u} and GPD shape \code{xi}.  
#'
#' @inheritParams normgpd
#' 
#' @details Extreme value mixture model combining normal distribution for the bulk
#' below the threshold and GPD for upper tail which is continuous at the threshold, 
#' but with one important difference to all the other mixture models.
#' 
#' The hybrid Pareto does not include the usual tail fraction \code{phiu} scaling,
#' i.e. so the GPD is not treated as a conditional model for the exceedances. 
#' The unscaled GPD is simply spliced
#' with the normal truncated at the threshold, with no rescaling to account for the proportion
#' above the threshold being applied. The parameters have to adjust for the lack of tail
#' fraction scaling.
#' 
#' The continuity constraint leads to the GPD scale \code{sigmau} being replaced
#' by a function of the normal mean, standard deviation and threshold parameters. 
#'  
#' The cumulative distribution function defined upto the 
#' threshold \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x) / r }
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = (H(u) +  G(x)) / r }
#' where \eqn{H(x)} and \eqn{G(X)} are the normal and conditional GPD
#' cumulative distribution functions (i.e. \code{pnorm(x, nmean, nsd)},
#' \code{pgpd(x, u, sigmau, xi)}). The normalisation constant \eqn{r} ensures a proper
#' density and is given by\code{r = 1 + pnorm(u, mean = nmean, sd = nsd)}, i.e. the 1 comes from
#' integration of the unscaled GPD and the second term is from the usual normal component.
#' 
#' The continuity constraint on the density at the threshold means that 
#' \eqn{h(u) = g(u)} where \eqn{h(x)} and \eqn{g(x)} are the normal and unscale GPD
#' density functions (i.e. \code{dnorm(x, nmean, nsd)} and
#' \code{dgpd(x, u, sigmau, xi)}).
#' 
#' The continuity constraint on its first derivative at the threshold 
#' means that \eqn{ h'(u) = g'(u)}. Then the Lambert W function is used for replacing
#' the threshold u and GPD scale sigmau in terms of the normal mean, standard deviation
#' and GPD shape xi.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[stats:Normal]{dnorm}} for details of normal bulk component.
#' 
#' @return \code{\link[evmix:hpdcon]{dhpdcon}} gives the density, 
#' \code{\link[evmix:hpdcon]{phpdcon}} gives the cumulative distribution function,
#' \code{\link[evmix:hpdcon]{qhpdcon}} gives the quantile function and 
#' \code{\link[evmix:hpdcon]{rhpdcon}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{rhpd} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:hpd]{rhpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Carreau, J. and Y. Bengio (2008). A hybrid Pareto model for asymmetric fat-tailed data:
#' the univariate case. Extremes 12 (1), 53-76.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}}, \code{\link[stats:Normal]{dnorm}}, 
#' \code{\link[evmix:normgpd]{dnormgpd}} and \code{\link[evmix:normgpdcon]{dnormgpdcon}}. 
#' The \code{\link[condmixt:condmixt-package]{condmixt}} package written by one of the
#' original authors of the hybrid Pareto model (Carreau and Bengio, 2008) also has 
#' similar functions for the hybrid Pareto \code{\link[condmixt:hpareto]{hpareto}} and
#' mixture of hybrid Paretos \code{\link[condmixt:hparetomixt]{hparetomixt}}, which are
#' more flexible as they also permit the model to be truncated at zero.
#' 
#' @aliases  hpdcon dhpdcon phpdcon qhpdcon rhpdcon
#' @family   hpdcon
#' 
#' @examples
#' \dontrun{
#' par(mfrow = c(2, 2))
#' xx = seq(-5, 20, 0.01)
#' f1 = dhpdcon(xx, nmean = 0,1, nsd = 0.4, u=0.5)
#' plot(xx, f1, type = "l")
#' 
#' # three tail behaviours
#' plot(xx, phpdcon(xx), type = "l")
#' lines(xx, phpdcon(xx, xi = 0.3), col = "red")
#' lines(xx, phpdcon(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#'  
#' sim = rhpdcon(1000, nmean = 0, nsd = 1.5, u = 0.5, xi = 0.2)
#' hist(sim, freq = FALSE, 100, xlim = c(-5, 20), ylim = c(0, 0.2))
#' lines(xx, dhpdcon(xx, nmean = 0, nsd = 1.5, u = 0.5, xi = 0.2), col = "blue")
#' 
#' plot(xx, dhpdcon(xx, nmean = 0, nsd = 1.5, u = 0.5, xi = 0), type = "l")
#' lines(xx, dhpdcon(xx, nmean = 0, nsd = 1.5, u = 0.5, xi = 0.2), col = "red")
#' lines(xx, dhpdcon(xx, nmean = 0, nsd = 1.5, u = 0.5, xi = -0.2), col = "blue")
#' legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   }
NULL

#' @export
#' @aliases  hpdcon dhpdcon phpdcon qhpdcon rhpdcon
#' @rdname hpdcon

# probability density function for hybrid Pareto model with a single continuity constraint
dhpdcon <- function(x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), xi = 0, log = FALSE) {
  
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
  linputs = c(length(x), length(nmean), length(nsd), length(u), length(xi))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" | mode(u) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, u, xi))))
    stop("parameters must be numeric")
  
  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")
  
  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")
      
  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  sigmau = 1/dnorm(u, mean = nmean, sd = nsd)
  
  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  r = 1 + pnorm(u, mean = nmean, sd = nsd)
  
  d = x # this will pass through NA/NaN in x just as they are entered
  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)
  
  if (nb > 0) d[whichb] = dnorm(x[whichb], mean = nmean[whichb], sd = nsd[whichb], log = TRUE) - log(r[whichb])
  if (nu > 0) d[whichu] = dgpd(x[whichu], u[whichu], sigmau[whichu], xi[whichu], log = TRUE) - log(r[whichu])
  
  if (!log) d = exp(d)
  
  d
  
}

