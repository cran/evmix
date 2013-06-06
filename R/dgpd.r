#' @name gpd
#' 
#' @title Generalised Pareto Distribution (GPD)
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the GPD conditional on being above a threshold
#'   \code{u} with parameters \code{sigmau} and \code{xi}. Unconditional
#'   quantities are provided when the probability \code{phiu} of being above the
#'   threshold \code{u} is given.
#'
#' @param x          quantile
#' @param q          quantile
#' @param p          cumulative probability
#' @param n          sample size (non-negative integer)
#' @param u          threshold
#' @param sigmau     scale parameter (non-negative)
#' @param xi         shape parameter
#' @param phiu       probability of being above threshold [0,1]
#' @param log        logical, if TRUE then log density
#' @param lower.tail logical, if FALSE then upper tail probabilities
#' 
#' @details The GPD with parameters scale \eqn{\sigma_u} and shape \eqn{\xi} has
#' conditional density given by
#'
#'  \deqn{f(x | X > u) = 1/\sigma_u [1 + \xi(x - u)/\sigma_u]^{-1/\xi - 1}}
#'  
#' for non-zero \eqn{\xi}, \eqn{x > u} and \eqn{\sigma_u > 0}. Further, 
#' \eqn{[1+\xi (x - u) / \sigma_u] > 0} which for \eqn{\xi < 0} implies 
#' \eqn{u < x \le u - \sigma_u/\xi}. In the special case of \eqn{\xi = 0}, which is
#' treated as \eqn{|\xi|<1e-6}, it reduces to the exponential:
#' 
#' \deqn{f(x | X > u) = 1/\sigma_u exp(-(x - u)/\sigma_u).}
#' 
#' The unconditional density is obtained by mutltiplying this by the survival probability
#' (or \emph{tail fraction}) \eqn{\phi_u = P(X > u)} giving \eqn{f(x) = \phi_u f(x | X > u)}.
#' 
#' The syntax of these functions are similar to those of the 
#' \code{\link[evd:gpd]{evd}} package, so most code using these functions can
#' simply be reused. The key difference is the introduction of \code{phiu} to
#' permit output of unconditional quantities.
#' 
#' @return \code{\link[evmix:gpd]{dgpd}} gives the density, \code{\link[evmix:gpd]{pgpd}} 
#' gives the cumulative distribution function,
#' \code{\link[evmix:gpd]{qgpd}} gives the quantile function and 
#' \code{\link[evmix:gpd]{rgpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{rgpd} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default threshold \code{u=0} and tail fraction
#' \code{phiu=1} which assumes the user will default to inputting excesses above 
#' \code{u}, rather than exceedance. The default sample size for 
#' \code{\link[evmix:gpd]{rgpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
#' 
#' Some key differences arise for \code{phiu=1} and \code{phiu<1} (see examples below):
#' 
#' \enumerate{
#' \item For \code{phiu=1} the \code{\link[evmix:gpd]{dgpd}} evaluates as zero for
#' quantiles below the threshold \code{u} and  \code{\link[evmix:gpd]{pgpd}}
#' evaluates over \eqn{[0, 1]}.
#' 
#' \item For \code{phiu=1} then \code{\link[evmix:gpd]{pgpd}} evaluates as zero
#' below the threshold \code{u}. For \code{phiu<1} it evaluates as \eqn{1-\phi_u} at
#' the threshold and \code{NA} below the threshold.
#' 
#' \item For \code{phiu=1} the quantiles from \code{\link[evmix:gpd]{qgpd}} are
#' above threshold and equal to threshold for \code{phiu=0}. For \code{phiu<1} then
#' within upper tail,  \code{p > 1 - phiu}, it will give conditional quantiles
#' above threshold, but when below the threshold, \code{p <= 1 - phiu}, these
#' are set to \code{NA}.
#' 
#' \item When simulating GPD variates using \code{\link[evmix:gpd]{rgpd}} if
#' \code{phiu=1} then all values are above the threshold. For \code{phiu<1} then
#' a standard uniform \eqn{U} is simulated and the variate will be classified as
#' above the threshold if \eqn{U<\phi}, and below the threshold otherwise. This is
#' equivalent to a binomial random variable for simulated number of exceedances. Those
#' above the threshold are then simulated from the conditional GPD and those below
#' the threshold and set to \code{NA}.
#' }
#' 
#' These conditions are intuitive and consistent with \code{\link[evd:gpd]{evd}},
#' which assumes missing data are below threshold.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' Based on GPD functions in the \code{\link[evd:gpd]{evd}} package.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evd:gpd]{evd}} and \code{\link[evd:fpot]{fpot}}
#' @aliases  gpd dgpd pgpd qgpd rgpd
#' @family   gpd
#' 
#' @examples
#' par(mfrow=c(2,2))
#' x = rgpd(1000) # simulate sample from GPD
#' xx = seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgpd(xx))
# 
#' # three tail behaviours
#' plot(xx, pgpd(xx), type = "l")
#' lines(xx, pgpd(xx, xi = 0.3), col = "red")
#' lines(xx, pgpd(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' # GPD when xi=0 is exponential, and demonstrating phiu
#' x = rexp(1000)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgpd(xx, u = 0, sigmau = 1, xi = 0), lwd = 2)
#' lines(xx, dgpd(xx, u = 0.5, phiu = 1 - pexp(0.5)), col = "red", lwd = 2)
#' lines(xx, dgpd(xx, u = 1.5, phiu = 1 - pexp(1.5)), col = "blue", lwd = 2)
#' legend("topright", paste("u =",c(0, 0.5, 1.5)),
#'   col=c("black", "red", "blue"), lty = 1, lwd = 2)
#' 
#' # Quantile function and phiu
#' p = pgpd(xx)
#' plot(qgpd(p), p, type = "l")
#' lines(xx, pgpd(xx, u = 2), col = "red")
#' lines(xx, pgpd(xx, u = 5, phiu = 0.2), col = "blue")
#' legend("bottomright", c("u = 0 phiu = 1","u = 2 phiu = 1","u = 5 phiu = 0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
NULL

#' @export
#' @aliases  gpd dgpd pgpd qgpd rgpd
#' @rdname gpd

# probability density function for GPD
dgpd <- function(x, u = 0, sigmau = 1, xi = 0, phiu = 1, log = FALSE) {

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
  linputs = c(length(x), length(u), length(sigmau), length(xi), length(phiu))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

  if (mode(u) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric" | 
    mode(phiu) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (min(sigmau) <= 0)
    stop("scale must be non-negative")

  if (any(phiu < 0) | any(phiu > 1))
    stop("phiu must between 0 and 1 (inclusive)")

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")
  
  x = rep(x, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)
  
  yu = (x - u)/sigmau # used when shape is zero
  syu = 1 + xi*yu     # used when shape non-zero
  
  # check for x values in range
  yind = ((yu > 0) & (syu > 0)) 

  d = x # this will pass through NA/NaN in x just as they are entered
  d[which(!yind)]= log(0) # zero density is default
  
  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)

  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    d[whichexp] = -log(sigmau[whichexp]) - yu[whichexp]
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    d[whichxi] = -log(sigmau[whichxi]) - (1/xi[whichxi] + 1) * log(syu[whichxi])
  }
  
  d = d + log(phiu) # unconditional density
  
  if (!log) d = exp(d)
  
  d
}
