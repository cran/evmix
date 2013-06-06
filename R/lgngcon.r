#' @export
#' 
#' @title Log-likelihood of Normal Bulk with GPD Upper and Lower Tails Extreme Value Mixture Model with Continuity Constraints
#'
#' @description Log-likelihood and negative log-likelihood for the extreme value 
#' mixture model with normal for bulk distribution between the lower and upper
#' thresholds with conditional GPD for the two tails with continuity constraints.
#'
#' @inheritParams lgng
#' 
#' @details The likelihood functions for the extreme value mixture model with
#' normal bulk and GPD for the two tails, as used in the maximum likelihood
#' fitting function \code{\link[evmix:fgngcon]{fgngcon}}.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:fgngcon]{fgngcon}}
#' but are available for wider usage, e.g. constructing your own extreme value
#' mixture models.
#' 
#' See \code{\link[evmix:fgngcon]{fgngcon}}, \code{\link[evmix:gngcon]{gngcon}} and
#' \code{\link[evmix:fgpd]{fgpd}} for full details.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:lgngcon]{lgngcon}}, which takes parameters as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:lgngcon]{lgngcon}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input). The tail fractions \code{phiul} and \code{phiur} are treated separately
#' to the other parameters, to allow for all it's representations.
#' 
#' Unlike the distribution functions \code{\link[evmix:gngcon]{gngcon}} the 
#' \code{phiu} must be either logical (\code{TRUE} or \code{FALSE}) or numerical in
#' range \eqn{(0, 1)}. The default is to specify \code{phiu=TRUE} so that the 
#' tail fraction is specified by normal distribution \eqn{\phi_u = 1 - H(u)}, or 
#' \code{phiu=FALSE} to treat the tail fraction as an extra parameter estimated using the 
#' sample proportion. Specify a numeric \code{phiu} as pre-specified probability
#' \eqn{(0, 1)}. Notice that the tail fraction probability cannot be 0 or 1 otherwise 
#' there would be no contribution from either tail or bulk components respectively.
#' 
#' The function \code{\link[evmix:lgngcon]{lgngcon}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @return \code{\link[evmix:lgngcon]{lgngcon}} gives (log-)likelihood and 
#' \code{\link[evmix:lgngcon]{nlgngcon}} gives the negative log-likelihood.
#'  
#' @note Unlike all the distribution functions for this mixture model, the 
#' likelihood functions only permits a scalar value for all the parameters. Only
#' the data is a vector.
#' 
#' Default values for the tail fractions \code{phiul=TRUE} and \code{phiur=TRUE}
#' is given in both \code{\link[evmix:lgngcon]{lgngcon}} and \code{\link[evmix:lgngcon]{nlgngcon}}. 
#' The \code{\link[evmix:lgngcon]{lgngcon}} also has the usual defaults for
#' the other parameters, but \code{\link[evmix:lgngcon]{nlgngcon}} has no defaults.
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
#' Zhao, X., Scarrott, C.J. Reale, M. and Oxley, L. (2010). Extreme value modelling
#' for forecasting the market crisis. Applied Financial Econometrics 20(1), 63-72.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:lgng]{lgng}}, \code{\link[evmix:lnormgpd]{lnormgpd}},
#' \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @aliases  lgngcon nlgngcon
#' @family gngcon

# log-likelihood function for normal bulk with GPD's for upper and lower tails
# will not stop evaluation unless it has to
lgngcon <- function(x, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE, log = TRUE) {
  
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
  lparams = c(length(nmean), length(nsd),
              length(ul), length(xil), length(phiul),
              length(ur), length(xir), length(phiur))
  
  if (sum(lparams != 1) > 0)
    stop("at least one parameter is not a scalar")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" |
    mode(ul) != "numeric" | mode(xil) != "numeric" |
    mode(ur) != "numeric" | mode(xir) != "numeric" |
    !(mode(phiul) %in% c("logical","numeric")) | 
    !(mode(phiur) %in% c("logical","numeric")))
    stop("parameters must be numeric, phiu can be numeric or logical")
  
  if (any(!is.finite(c(nmean, nsd, ul, xil, phiul, ur, xir, phiur))))
    stop("parameters must be numeric")
  
  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")
  
  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these
  
  xul = x[which(x < ul)]
  nul = length(xul)
  xb = x[which((x <= ur) & (x >= ul))]
  nb = length(xb)
  xur = x[which(x > ur)]
  nur = length(xur)
  n = nul + nb + nur
  
  if ((nsd <= 0) | (ul <= min(x)) | (ur >= max(x)) | (ul >= ur)) {
    l = -Inf
  } else {
    if (is.logical(phiul)) {
      if (phiul) {
        phiul = pnorm(ul, mean = nmean, sd = nsd)
      } else {
        phiul = nul / n
      }
    }
    if (is.logical(phiur)) {
      if (phiur) {
        phiur = 1 - pnorm(ur, mean = nmean, sd = nsd)
      } else {
        phiur = nur / n
      }
    }
    phib = (1 - phiul - phiur) / (pnorm(ur, mean = nmean, sd = nsd) - pnorm(ul, mean = nmean, sd = nsd))
    
    dul = dnorm(ul, mean = nmean, sd = nsd)
    dur = dnorm(ur, mean = nmean, sd = nsd)

    sigmaul = phiul / (phib * dul)
    sigmaur = phiur / (phib * dur)
        
    syur = 1 + xir * (xur - ur) / sigmaur  
    syul = 1 + xil * (ul - xul) / sigmaul  
    yb = (xb - nmean) / nsd    # used for normal
    
    if ((min(syul) <= 0) | (min(syur) <= 0) | 
      (phib < .Machine$double.eps) | (dul < .Machine$double.eps) | (dur < .Machine$double.eps) | 
      (sigmaul <= 0) | (sigmaur <= 0) | (phiul <= 0) | (phiul >= 1) |
      (phiur <= 0) | (phiur >= 1) | ((phiul + phiur) > 1)) {
      l = -Inf
    } else { 
      l = lgpd(-xul, -ul, sigmaul, xil, phiul)
      l = l - nb * log(2 * pi * nsd ^ 2) / 2 - sum(yb ^ 2) / 2 + nb * log(phib)
      l = l + lgpd(xur, ur, sigmaur, xir, phiur)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}
