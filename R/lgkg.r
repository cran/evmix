#' @export
#' 
#' @title Cross-validation Log-likelihood of Kernel Density Estimation for Bulk and GPD
#' for Both Upper and Lower Tails in Extreme Value Mixture Model
#'
#' @description Cross-validation log-likelihood and negative log-likelihood for the
#' kernel density estimation using normal kernel bulk and GPD upper and lower tails extreme value mixture model.
#'
#' @inheritParams fgkg
#' @inheritParams gkg
#' 
#' @details The cross-validation likelihood functions for the extreme value mixture model
#' with kernel density estimation using normal kernel for bulk distribution between the upper
#' and lower thresholds with conditional GPD's for the two tails. As used in the maximum likelihood fitting function 
#'  \code{\link[evmix:fgkg]{fgkg}}.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:fgkg]{fgkg}}
#' but are available for wider usage, e.g. constructing your own extreme value
#' mixture models.
#' 
#' See \code{\link[evmix:fkdengpd]{fkdengpd}}, \code{\link[evmix:fkden]{fkden}}
#' and \code{\link[evmix:fgpd]{fgpd}} for full details.
#' 
#' Cross-validation likelihood is used for kernel density component, but 
#' standard likelihood is used for GPD components. The cross-validation likelihood for
#' the KDE is obtained by leaving each point out in turn,
#' evaluating the KDE at the point left out:
#'    \deqn{L(\lambda)\prod_{i=1}^{nb} \hat{f}_{-i}(x_i)}
#' where 
#'    \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(n-1)\lambda} \sum_{j=1: j\ne i}^{n} K(\frac{x_i - x_j}{\lambda})}
#' is the KDE obtained when the \eqn{i}th datapoint is dropped out and then 
#' evaluated at that dropped datapoint
#' at \eqn{x_i}. Notice that the KDE sum is indexed over all datapoints (\eqn{j=1, ..., n},
#' except datapoint \eqn{i}) whether they are between the thresholds or in the tails. But the
#' likelihood product is evaluated only for those data between the thresholds 
#' (\eqn{i=1, ..., n_b}). So the \eqn{j = n_b+1, ..., n} datapoint are extra kernel centres
#' from the data in the tails which are used in the KDE but the likelihood is not evaluated there.
#' 
#' Log-likelihood calculations are carried out in \code{\link[evmix:lgkg]{lgkg}}, 
#' which takes bandwidth in the same form as distribution functions. The negative 
#' log-likelihood is a wrapper for \code{\link[evmix:lgkg]{lgkg}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input).
#' 
#' The function \code{\link[evmix:lgkg]{lgkg}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @section Warning:
#' See warning in \code{\link[evmix:fkden]{fkden}}
#'  
#' @return \code{\link[evmix:lgkg]{lgkg}} gives cross-validation (log-)likelihood and 
#' \code{\link[evmix:lgkg]{nlgkg}} gives the negative cross-validation log-likelihood.
#'  
#' @note Invalid bandwidth parameter will give \code{0} for likelihood, \code{log(0)=-Inf} for
#' cross-validation log-likelihood and \code{-log(0)=Inf} for negative cross-validation log-likelihood. 
#' 
#' See \code{\link[evmix:gpd]{fgpd}} for explanation of \code{finitelik}.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Bowman, A.W. (1984). An alternative method of cross-validation for the smoothing of
#' density estimates. Biometrika 71(2), 353-360.
#' 
#' Duin, R.P.W. (1976). On the choice of smoothing parameters for Parzen estimators of
#' probability density functions. IEEE Transactions on Computers C25(11), 1175-1179.
#' 
#' MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B., Reale, M. and Russell, G. (2011).
#' A flexible extreme value mixture model. Computational Statistics and Data Analysis
#' 55(6), 2137-2157.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gkg]{gkg}}, \code{\link[evmix:kdengpd]{kdengpd}}, 
#' \code{\link[evmix:kden]{kden}}, \code{\link[evmix:gpd]{gpd}} and 
#' \code{\link[stats:density]{density}}.
#' @aliases lgkg nlgkg
#' @family  gkg

# cross-validation log-likelihood function for KDE for bulk with GPD's for both upper and lower tails
# will not stop evaluation unless it has to
lgkg <- function(x, lambda = NULL, 
  ul = as.vector(quantile(x, 0.1)), sigmaul = 1, xil = 0, phiul = TRUE, 
  ur = as.vector(quantile(x, 0.9)), sigmaur = 1, xir = 0, phiur = TRUE, log = TRUE) {
  
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    stop("infinite cases must be removed")
  
  if (any(is.na(x)))
    warning("missing values have been removed")
  
  x = x[which(!is.na(x))]
  n = length(x)
  
  if (is.null(lambda)){
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10){
      stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    lambda = bw.nrd0(x) # do not extra kernel centres
  }
  
  # parameter inputs inputs for likelihood must be single values
  lparams = c(length(lambda), length(ul), length(sigmaul), length(xil), length(phiul),
                              length(ur), length(sigmaur), length(xir), length(phiur))
  
  if (sum(lparams != 1) > 0)
    stop("at least one parameter is not a scalar")
  
  if (mode(lambda) != "numeric"| 
      mode(ul) != "numeric" | mode(sigmaul) != "numeric" | mode(xil) != "numeric" |
      mode(ur) != "numeric" | mode(sigmaur) != "numeric" | mode(xir) != "numeric" |
      !(mode(phiul) %in% c("logical","numeric")) | !(mode(phiur) %in% c("logical","numeric")))
    stop("parameters must be numeric, phiu can be numeric or logical")
  
  if (any(!is.finite(c(lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur))))
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
  
  if ((lambda <= 0) | (sigmaul <= 0) | (sigmaur <= 0) |
    (ul <= min(x)) | (ur >= max(x)) | (ul >= ur)) {
    l = -Inf
  } else {
    if (is.logical(phiul)) {
      phiul = pkdenx(ul, x, lambda)
    } else {
      phiul = nul / n
    }
    if (is.logical(phiur)) {
      phiur = 1 - pkdenx(ur, x, lambda)
    } else {
      phiur = nur / n
    }
    phib = (1 - phiul - phiur) / (pkdenx(ur, x, lambda) - pkdenx(ul, x, lambda))
    
    syur = 1 + xir * (xur - ur) / sigmaur  
    syul = 1 + xil * (ul - xul) / sigmaul  
    
    if ((min(syul) <= 0) | (min(syur) <= 0) |
      (phiul <= 0) | (phiul >= 1) |
      (phiur <= 0) | (phiur >= 1) | ((phiul + phiur) > 1)) {
      l = -Inf
    } else {                     
      l = lgpd(-xul, -ul, sigmaul, xil, phiul)
      l = l + lgpd(xur, ur, sigmaur, xir, phiur)
      l = l + lkden(xb, lambda, extracentres = c(xul, xur), log = TRUE) + nb*log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

