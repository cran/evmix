#' @export
#' 
#' @title Log-likelihood of Generalised Pareto Distribution (GPD)
#'
#' @description Log-likelihood and negative log-likelihood for the GPD conditional on being
#' above a threshold \code{u} with parameters scale \code{sigmau} and shape \code{xi}.
#' Unconditional likelihood also provided when the probability \code{phiu} of being above the
#' threshold \code{u} is given.
#'
#' @inheritParams gpd
#' @inheritParams fgpd
#' 
#' @details The GPD likelihood functions for the exceedances of the threshold \code{u}
#' as used in the maximum likelihood fitting function \code{\link[evmix:fgpd]{fgpd}}.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:fgpd]{fgpd}} but are available for
#' wider usage, e.g. constructing your own extreme value mixture models.
#' 
#' See \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{dgpd}} for full details.
#' 
#' Log-likelihood calculations are carried out in \code{\link[evmix:lgpd]{lgpd}},
#' which takes parameters as inputs in the same form as distribution functions.
#' The negative log-likelihood is a wrapper for \code{\link[evmix:lgpd]{lgpd}}, designed towards
#' making it useable for optimisation (e.g. parameters are given a vector as first input).
#' 
#' Unlike the MLE function \code{\link[evmix:fgpd]{fgpd}}, the \code{phiu} must be
#' in range \eqn{[0, 1]} and cannot be \code{NULL}. Specify \code{phiu=1} for conditional
#' likelihood (default) and \code{phiu<1} for unconditional likelihood.
#' 
#' The function \code{\link[evmix:lgpd]{lgpd}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @return \code{\link[evmix:lgpd]{lgpd}} gives (log-)likelihood and 
#' \code{\link[evmix:lgpd]{nlgpd}} gives the negative log-likelihood.
#'  
#' @note Unlike all the distribution functions for the GPD, the likelihood functions
#' only permits a scalar value for scale and shape parameters, \code{phiu} and
#' threshold \code{u}. Only the data is a vector.
#' 
#' Default values for the threshold \code{u=0} and tail fraction \code{phiu=1} are
#' given in both \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:lgpd]{nlgpd}},
#' assuming the user will default to entering excesses above the threshold,
#' rather than exceeedances. The \code{\link[evmix:lgpd]{lgpd}} has the usual defaults
#' for the GPD scale and shape parameters, but \code{\link[evmix:lgpd]{nlgpd}} has no defaults.
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
#' Based on GPD likelihood function in the \code{\link[evd:gpd]{evd}} package.
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evd:gpd]{dgpd}}, \code{\link[evd:fpot]{fpot}}
#' and \code{\link[MASS:fitdistr]{fitdistr}}
#' @aliases  lgpd nlgpd
#' @family   gpd

# log-likelihood function for GPD
# will not stop evaluation unless it has to
lgpd <- function(x, u = 0, sigmau = 1, xi = 0, phiu = 1, log = TRUE) {

  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    stop("infinite cases must be removed")

  # parameter inputs inputs for likelihood must be single values
  lparams = c(length(u), length(sigmau), length(xi), length(phiu))

  if (sum(lparams != 1) > 0)
    stop("at least one parameter is not a scalar")

  if (mode(u) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric" | 
    mode(phiu) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  # assume NA or NaN are below threshold consistent with evd library
  # hence use which() to ignore these
  
  # check for x values in range
  yind = x > u
  
  xu = x[which(yind)]
  nu = length(xu)
  
  yu = (xu - u) / sigmau # used when shape is zero
  syu = 1 + xi * yu      # used when shape non-zero  
    
  if ((min(syu) <= 0) | (sigmau <= 0) | (phiu < 0)  | (phiu > 1)) {
    l = -Inf
  } else {
    if (abs(xi) < 1e-6) {
      l = - nu * log(sigmau) - sum(yu) + nu * log(phiu)
    } else {
      l = - nu * log(sigmau) - (1/xi + 1) * sum(log(syu)) + nu * log(phiu)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

