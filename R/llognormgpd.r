#' @export
#' 
#' @title Log-likelihood of Log-Normal Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Log-likelihood and negative log-likelihood for the extreme value 
#' mixture model with log-normal for bulk distribution upto the threshold and conditional
#' GPD above threshold.
#'
#' @param phiu probability of being above threshold [0,1] or logical
#' @inheritParams flognormgpd
#' @inheritParams lognormgpd
#' 
#' @details The likelihood functions for the extreme value mixture model with
#' log-normal bulk and GPD tail, as used in the maximum likelihood fitting function 
#' \code{\link[evmix:flognormgpd]{flognormgpd}}.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:flognormgpd]{flognormgpd}}
#' but are available for wider usage, e.g. constructing your own extreme value
#' mixture models.
#' 
#' Negative data are ignored.
#' 
#' See \code{\link[evmix:flognormgpd]{flognormgpd}}, \code{\link[evmix:fnormgpd]{fnormgpd}}
#' and \code{\link[evmix:fgpd]{fgpd}} for full details.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:llognormgpd]{llognormgpd}}, which takes parameters as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:llognormgpd]{llognormgpd}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input). The tail fraction \code{phiu} is treated separately to the other parameters,
#' to allow for all it's representations.
#' 
#' Unlike the distribution functions \code{\link[evmix:lognormgpd]{lognormgpd}} the 
#' \code{phiu} must be either logical (\code{TRUE} or \code{FALSE}) or numerical in
#' range \eqn{(0, 1)}. The default is to specify \code{phiu=TRUE} so that the 
#' tail fraction is specified by log-normal distribution \eqn{\phi_u = 1 - H(u)}, or 
#' \code{phiu=FALSE} to treat the tail fraction as an extra parameter estimated using the 
#' sample proportion. Specify a numeric \code{phiu} as pre-specified probability
#' \eqn{(0, 1)}. Notice that the tail fraction probability cannot be 0 or 1 otherwise 
#' there would be no contribution from either tail or bulk components respectively.
#' 
#' The function \code{\link[evmix:llognormgpd]{llognormgpd}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @return \code{\link[evmix:llognormgpd]{llognormgpd}} gives (log-)likelihood and 
#' \code{\link[evmix:llognormgpd]{nllognormgpd}} gives the negative log-likelihood.
#'  
#' @note Unlike all the distribution functions for this mixture model, the 
#' likelihood functions only permits a scalar value for all the parameters. Only
#' the data is a vector.
#' 
#' A default value for the tail fraction \code{phiu=TRUE} is given in both 
#' \code{\link[evmix:llognormgpd]{llognormgpd}} and \code{\link[evmix:llognormgpd]{nllognormgpd}}. 
#' The \code{\link[evmix:llognormgpd]{llognormgpd}} also has the usual defaults for
#' the other parameters, but \code{\link[evmix:llognormgpd]{nllognormgpd}} has no defaults.
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
#' \url{http://en.wikipedia.org/wiki/Log-normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Solari, S. and Losada, M.A. (2004). A unified statistical model for
#' hydrological variables including the selection of threshold for the peak over
#' threshold method. Water Resources Research. 48, W10541.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @aliases  llognormgpd nllognormgpd
#' @family   lognormgpd

# log-likelihood function for log-normal bulk with GPD for upper tail
# will not stop evaluation unless it has to
llognormgpd <- function(x, lnmean = 0, lnsd = 1, u = qlnorm(0.9, lnmean, lnsd),
  sigmau = lnsd, xi = 0, phiu = TRUE, log = TRUE) {

  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    stop("infinite cases must be removed")

  if (any(is.na(x)))
    warning("missing values have been removed")

  if (any(x <= 0))
    warning("non-positive values have been removed")

  # parameter inputs inputs for likelihood must be single values
  lparams = c(length(lnmean), length(lnsd), length(u), length(sigmau), length(xi), length(phiu))

  if (sum(lparams != 1) > 0)
    stop("at least one parameter is not a scalar")

  if (mode(lnmean) != "numeric" | mode(lnsd) != "numeric" | mode(u) != "numeric" |
      mode(sigmau) != "numeric" | mode(xi) != "numeric" |
      !(mode(phiu) %in% c("logical","numeric")))
    stop("parameters must be numeric, phiu can be numeric or logical")
  
  if (any(!is.finite(c(lnmean, lnsd, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  xu = x[which(x > u)]
  nu = length(xu)
  xb = x[which((x <= u) & (x > 0))]
  nb = length(xb)
  n = nb + nu

  if ((lnsd <= 0) | (sigmau <= 0) | (u <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.logical(phiu)) {
      if (phiu) {
        phiu = 1 - plnorm(u, meanlog = lnmean, sdlog = lnsd)
      } else {
        phiu = nu / n
      }
    }
    phib = (1 - phiu) / plnorm(u, meanlog = lnmean, sdlog = lnsd)
  
    syu = 1 + xi * (xu - u) / sigmau  
    yb = (log(xb) - lnmean) / lnsd    # used for normal
  
    if ((min(syu) <= 0) | (phiu <= 0) | (phiu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l - sum(log(xb)) - nb * log(2 * pi * lnsd ^ 2) / 2 - sum(yb ^ 2) / 2 + nb * log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}
