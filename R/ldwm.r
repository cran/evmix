#' @export
#' 
#' @title Log-likelihood of dynamically weighted mixture model
#'
#' @description Log-likelihood and negative log-likelihood for the dynamically
#' weighted mixture model
#'
#' @inheritParams fdwm
#' @inheritParams dwm
#' 
#' @details The likelihood functions for the dynamically weighted mixture model
#' \code{\link[evmix:fdwm]{fdwm}}.
#' 
#' Non-positive data are ignored.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:fdwm]{fdwm}}
#' but are available for wider usage, e.g. constructing your own extreme value
#' mixture models.
#' 
#' See \code{\link[evmix:fdwm]{fdwm}} and \code{\link[evmix:fgpd]{fgpd}}
#' for full details.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:ldwm]{ldwm}}, which takes parameters as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:ldwm]{ldwm}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input). 
#' 
#' The function \code{\link[evmix:ldwm]{ldwm}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @return \code{\link[evmix:ldwm]{ldwm}} gives (log-)likelihood and 
#' \code{\link[evmix:ldwm]{nldwm}} gives the negative log-likelihood.
#'  
#' @note Unlike all the distribution functions for this mixture model, the 
#' likelihood functions only permits a scalar value for all the parameters. Only
#' the data is a vector.
#' 
#' Invalid parameters will give \code{0} for likelihood, \code{log(0)=-Inf} for
#' log-likelihood and \code{-log(0)=Inf} for negative log-likelihood. 
#' 
#' See \code{\link[evmix:gpd]{fgpd}} for explanation of \code{finitelik}.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
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
#' Frigessi, A., O. Haug, and H. Rue (2002). A dynamic mixture model for unsupervised tail
#' estimation without threshold selection. Extremes 5 (3), 219-235
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @aliases ldwm nldwm
#' @family  dwm

# log-likelihood function for dynamically weighted mixture model
# will not stop evaluation unless it has to
ldwm = function(x,  wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, log = TRUE) {
  
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    stop("infinite cases must be removed")
  
  if (any(is.na(x)))
    warning("missing values have been removed")
  
  if (any(x < 0))
    warning("negative values have been removed")
  
  # parameter inputs inputs for likelihood must be single values
  lparams = c(length(wshape), length(wscale), length(cmu), length(ctau), length(sigmau), length(xi))
  
  if (sum(lparams != 1) > 0)
    stop("at least one parameter is not a scalar")
  
  if (mode(wshape) != "numeric" | mode(wscale) != "numeric" | mode(cmu) != "numeric" |
        mode(ctau) != "numeric" | mode(sigmau) != "numeric"| mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(wshape, wscale, cmu, ctau, sigmau, xi))))
    stop("parameters must be numeric")

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  x = x[which(!is.na(x))]
  n = length(x)
  
  rx <- function(x, wshape, wscale, cmu, ctau, sigmau, xi) {
    (dgpd(x, sigmau = sigmau, xi = xi) - dweibull(x, shape = wshape, scale = wscale))*atan((x - cmu)/ctau)
  }

  if ((wscale <= 0) | (wshape <= 0) | (cmu <= 0) | (ctau<=0) | (sigmau <= 0) | (cmu >= max(x))) {
    l = -Inf
  } else {
        
    syu = 1 + xi * (x / sigmau)
    
    if ((min(syu) <= 0)) {
      l = -Inf
    } else { 
      px = pcauchy(x, location = cmu, scale = ctau)

      r = try(integrate(rx, wshape = wshape, wscale = wscale, 
        cmu = cmu, ctau = ctau, sigmau = sigmau, xi = xi,
        lower= 0, upper = Inf, subdivisions = 10000, rel.tol = 1.e-10, stop.on.error = FALSE)$value)
      
      if (inherits(r, "try-error")) {
        warning("numerical integration failed, ignore previous messages, optimisation will try again")
        l = -Inf
      } else {
        z = n*log((1 + r/pi))
        
        bulk = sum(log((1 - pcauchy(x, location = cmu, scale = ctau))*dweibull(x, shape = wshape, scale = wscale)+
          pcauchy(x, location = cmu, scale = ctau)*dgpd(x, sigmau = sigmau, xi = xi)))
        l = bulk-z
      }   
    }
  }
  
  l
}
