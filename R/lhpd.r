#' @export
#' 
#' @title Log-likelihood of Hybrid Pareto Extreme Value Mixture Model
#'
#' @description Log-likelihood and negative log-likelihood for the hybrid Pareto extreme
#' value mixture model
#'
#' @inheritParams lnormgpd
#' 
#' @details The likelihood functions for hybrid Pareto extreme value mixture model, 
#' as used in the maximum likelihood fitting function \code{\link[evmix:fhpd]{fhpd}}.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:fhpd]{fhpd}}
#' but are available for wider usage, e.g. constructing your own extreme value
#' mixture models.
#' 
#' See \code{\link[evmix:fhpd]{fhpd}} and \code{\link[evmix:fgpd]{fgpd}}
#' for full details.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:lhpd]{lhpd}}, which takes parameters as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:lhpd]{lhpd}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input).
#' 
#' The function \code{\link[evmix:lhpd]{lhpd}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @return \code{\link[evmix:lhpd]{lhpd}} gives (log-)likelihood and 
#' \code{\link[evmix:lhpd]{nlhpd}} gives the negative log-likelihood.
#'  
#' @note Unlike all the distribution functions for this mixture model, the 
#' likelihood functions only permits a scalar value for all the parameters. Only
#' the data is a vector.
#' 
#' The \code{\link[evmix:lhpd]{lhpd}} also has the usual defaults for
#' the other parameters, but \code{\link[evmix:lhpd]{nlhpd}} has no defaults.
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
#' @seealso \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:gpd]{gpd}}. 
#' The \code{\link[condmixt:condmixt-package]{condmixt}} package written by one of the
#' original authors of the hybrid Pareto model (Carreau and Bengio, 2008) also has 
#' similar functions for the likelihood of the hybrid Pareto 
#' \code{\link[condmixt:hpareto.negloglike]{hpareto.negloglike}} and
#' fitting \code{\link[condmixt:hpareto.negloglike]{hpareto.fit}}.
#' 
#' @aliases  lhpd nlhpd
#' @family   hpd

# log-likelihood function for hybrid Pareto
# will not stop evaluation unless it has to
lhpd <- function(x, nmean = 0, nsd = 1, xi = 0, log = TRUE) {
  
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    stop("infinite cases must be removed")
  
  if (any(is.na(x)))
    warning("missing values have been removed")
  
  # parameter inputs inputs for likelihood must be single values
  lparams = c(length(nmean), length(nsd), length(xi))
  
  if (sum(lparams != 1) > 0)
    stop("at least one parameter is not a scalar")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" | mode(xi) != "numeric" )
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, xi))))
    stop("parameters must be numeric")
  
  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")
  
  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these
 
  z = (1 + xi)^2/(2*pi)
  wz = lambertW(z)
  u = sign(1 + xi) * nsd * sqrt(wz) + nmean

  xu = x[which(x > u)]
  nu = length(xu)
  xb = x[which(x <= u)]
  nb = length(xb)
  n = nb + nu
    
  if ((nsd <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
  
    du = sqrt(wz)
    sigmau = nsd * sqrt((1 + xi)^2)/ du
    
    syu = 1 + xi * (xu - u) / sigmau  
    
    yb = (xb - nmean) / nsd    # used for normal
    
    r = 1 + pnorm(sign(1 + xi) * sqrt(lambertW((1 + xi)^2/(2*pi))))
    phiu = 1/r

    if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps) | (phiu <= 0) | (phiu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l - nb * log(2 * pi * nsd ^ 2) / 2 - sum(yb ^ 2) / 2 + nb * log(phiu)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

