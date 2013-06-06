#' @name dwm
#' 
#' @title Dynamically Weighted Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the dynamically weighted mixture model. The parameters
#'   are the Weibull shape \code{wshape} and scale \code{wscale}, Cauchy location 
#'   \code{cmu}, Cauchy scale \code{ctau},
#'   GPD scale \code{sigmau}, shape \code{xi} and  initial value for the quantile \code{qinit}.
#'
#' @inheritParams weibullgpd
#' @param cmu       Cauchy location
#' @param ctau      Cauchy scale
#' @param qinit     vector of initial values for the quantile estimate
#' 
#' @details The dynamic weighted mixture model combines a Weibull
#' for the bulk model with GPD for the tail model. However, unlike all the other mixture models
#' the GPD is defined over the entire range of support rather than as a conditional model
#' above some threshold. A transition function is used to apply weights to transition between
#' the bulk and GPD for the upper tail, thus providing the dynamically weighted mixture. They
#' use a Cauchy cumulative distribution function for the transition function.
#' 
#' The density function is then a dynamically weighted mixture given by:
#'  \deqn{f(x) = {[1 - p(x)] h(x) + p(x) g(x)}/r}
#' where \eqn{h(x)} and \eqn{g(x)} are the Weibull and unscaled GPD density functions
#' respectively (i.e. \code{pweibull(x, wshape, wscale)} and \code{pgpd(x, u, sigmau, xi)}).
#' The Cauchy cumulative distribution function used to provide the transition is defined by
#' \eqn{p(x)} (i.e. \code{pcauchy(x, cmu, ctau}. The normalisation constant \eqn{r} ensures a
#' proper density.
#' 
#' The quantile function is not available in closed form, so has to be solved numerically. 
#' The argument \code{qinit} is the initial quantile estimate which is used for numerical
#' optimisation and should be set to a reasonable guess. When the \code{qinit} is \code{NULL},
#' the initial quantile value is given by the midpoint between the Weibull and GPD
#' quantiles. As with the other inputs \code{qinit} is also vectorised, but \code{R} does not
#' permit vectors combining \code{NULL} and numeric entries.
#' 
#' @return \code{\link[evmix:dwm]{ddwm}} gives the density, 
#' \code{\link[evmix:dwm]{pdwm}} gives the cumulative distribution function,
#' \code{\link[evmix:dwm]{qdwm}} gives the quantile function and 
#' \code{\link[evmix:dwm]{rdwm}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:dwm]{rdwm}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:dwm]{rdwm}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Weibull_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Frigessi, A., Haug, O. and Rue, H. (2002). A dynamic mixture model for unsupervised tail
#' estimation without threshold selection. Extremes 5 (3), 219-235
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:weibullgpd]{weibullgpd}}, \code{\link[evmix:gpd]{gpd}}
#' and \code{\link[stats:Weibull]{dweibull}}
#' @aliases  dwm ddwm pdwm qdwm rdwm
#' @family   ldwm nldwm fdwm
#' 
#' @examples
#' \dontrun{
#' par(mfrow = c(2, 2))
#' xx = seq(0.001, 5, 0.01)
#' f = ddwm(xx, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.5)
#' plot(xx, f, ylim = c(0, 1), xlim = c(0, 5), type = 'l', lwd = 2, 
#'   ylab = "density", main = "Plot example in Frigessi et al. (2002)")
#' lines(xx, dgpd(xx, xi = 1, sigmau = 0.5), col = "red", lty = 2, lwd = 2)
#' lines(xx, dweibull(xx, shape = 2, scale = 1/gamma(1.5)), col = "blue", lty = 2, lwd = 2)
#' legend('topright', c('DWM', 'Weibull', 'GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#' 
#' # three tail behaviours
#' plot(xx, pdwm(xx, xi = 0), type = "l")
#' lines(xx, pdwm(xx, xi = 0.3), col = "red")
#' lines(xx, pdwm(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)), col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rdwm(10000, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.1)
#' xx = seq(0, 15, 0.01)
#' hist(x, freq = FALSE, breaks = 100)
#' lines(xx,ddwm(xx, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.1),
#'   lwd = 2, col = 'black')
#'   
#' plot(xx, pdwm(xx, wshape = 2, wscale = 1/gamma(1.5), cmu = 1, ctau = 1, sigmau = 1, xi = 0.1),
#'  xlim = c(0, 15), type = 'l', lwd = 2, 
#'   xlab = "x", ylab = "F(x)")
#' lines(xx, pgpd(xx, sigmau = 1, xi = 0.1), col = "red", lty = 2, lwd = 2)
#' lines(xx, pweibull(xx, shape = 2, scale = 1/gamma(1.5)), col = "blue", lty = 2, lwd = 2)
#' legend('bottomright', c('DWM', 'Weibull', 'GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#' }
NULL

#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# probability density function for dynamically weighted mixture model
ddwm = function(x, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, log = FALSE) {

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
  linputs = c(length(x), length(wshape), length(wscale), 
              length(cmu), length(ctau), length(sigmau), length(xi))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
  if (mode(wshape) != "numeric" | mode(wscale) != "numeric" | mode(cmu) != "numeric" |
    mode(ctau) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(wshape, wscale, cmu, ctau, sigmau, xi))))
    stop("parameters must be numeric")
  
  if (min(wshape) < 0)
    stop("weibull shape must be non-negative")
  
  if (min(wscale) < 0)
    stop("weibull scale must be non-negative")
  
  if (min(ctau) < 0)
    stop("Cauchy scale must be non-negative")
  
  if (min(sigmau) <= 0)
    stop("scale must be non-negative")
  
  x = rep(x, length.out = n)
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
    
  rx <- function(x, wshape, wscale, cmu, ctau, sigmau, xi) {
    (dgpd(x, sigmau = sigmau, xi = xi) - dweibull(x, shape = wshape, scale = wscale))*atan((x - cmu)/ctau)
  }
  
  d = x # this will pass through NA/NaN in x just as they are entered

  whichnonmiss = which(!is.na(x))

  if (max(linputs[-1]) == 1) {
    r = try(integrate(rx, wshape = wshape[1], wscale = wscale[1],
      cmu = cmu[1], ctau = ctau[1], sigmau = sigmau[1], xi = xi[1],
      lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         
  
    if (inherits(r, "try-error")) {
      z = rep(NA, length.out = n)
    } else {              
      z = rep(1 + r/pi, length.out = n)
    }
  } else {

    z = rep(NA, length.out = n)
    for (i in 1:n) {
      r = try(integrate(r, wshape = wshape[i], wscale = wscale[i],
        cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
        lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         

      if (inherits(r, "try-error")) {
        z[i] = NA
      } else {              
        z[i] = 1 + r/pi
      }
    }
  }
  
  d[whichnonmiss] = ((1 - 
      pcauchy(x[whichnonmiss], location = cmu[whichnonmiss], scale = ctau[whichnonmiss]))
    * dweibull(x[whichnonmiss], shape = wshape[whichnonmiss], scale = wscale[whichnonmiss])+
    pcauchy(x[whichnonmiss], location = cmu[whichnonmiss], scale = ctau[whichnonmiss]) *
      dgpd(x, sigmau = sigmau[whichnonmiss], xi = xi[whichnonmiss]))/z[whichnonmiss]
  
  if (log) d = log(d)

  d
}
