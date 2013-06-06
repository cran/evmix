#' @export
#' 
#' @title Cross-validation Log-likelihood of Boundary Corrected Kernel Density Estimation
#'
#' @description Cross-validation log-likelihood and negative log-likelihood for
#' boundary corrected kernel density estimation, by treating it as a mixture model.
#'
#' @inheritParams fbckden
#' @inheritParams bckden
#' @inheritParams fkden
#' 
#' @details The cross-validation likelihood functions for the boundary corrected
#' kernel density estimator, as used in the maximum likelihood fitting function 
#' \code{\link[evmix:fbckden]{fbckden}}.
#' 
#' They are designed to be used for MLE in \code{\link[evmix:fbckden]{fbckden}}
#' but are available for wider usage, e.g. constructing your own extreme value
#' mixture models.
#' 
#' All of the boundary correction methods available in \code{\link[evmix:bckden]{bckden}}
#' are permitted.
#' 
#' See \code{\link[evmix:fbckden]{fkden}} and \code{\link[evmix:fgpd]{fgpd}}
#' for full details.
#' 
#' The cross-validation likelihood is obtained by leaving each point out in turn, obtaining the
#' usual KDE and evaluate at the point left out:
#'    \deqn{L(\lambda)\prod_{i=1}^{n} \hat{f}_{-i}(x_i)}
#' where 
#'    \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(n-1)\lambda} \sum_{j=1: j\ne i}^{n} K(\frac{x_i - x_j}{\lambda})}
#' is the KDE obtained when the \eqn{i}th datapoint is dropped but is evaluated at \eqn{x_i}.
#' 
#' Normally for likelihood estimation of the bandwidth the kernel centres and the data 
#' where the likelihood is evaluated are the same. However, when using KDE for extreme value
#' mixture modelling the likelihood only those data in the bulk of the distribution should contribute
#' to the likelihood, but all the data (including those beyond the threshold) should contribute to
#' the density estimate. The \code{extracentres} option allows the use to specify extra kernel
#' centres used in estimating the density, but not evaluated in the likelihood. The default
#' is to just use the existing data, so \code{extracentres=NULL}.
#' 
#' Log-likelihood calculations are carried out in \code{\link[evmix:lbckden]{lbckden}}, 
#' which takes bandwidth in the same form as distribution functions. The negative 
#' log-likelihood is a wrapper for \code{\link[evmix:lbckden]{lbckden}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input).
#' 
#' The function \code{\link[evmix:lbckden]{lbckden}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' @section Warning:
#' See warning in \code{\link[evmix:fbckden]{fbckden}}
#'  
#' @return \code{\link[evmix:lbckden]{lbckden}} gives cross-validation (log-)likelihood and 
#' \code{\link[evmix:lbckden]{nbclkden}} gives the negative cross-validation log-likelihood.
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
#' @seealso \code{\link[stats:density]{density}}
#' @aliases lbckden nlbckden
#' @family  bckden

# cross-validation log-likelihood function for boundary corrected KDE
# will not stop evaluation unless it has to
lbckden <- function(x, lambda = NULL, extracentres = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, log = TRUE) {

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

  if (!is.null(extracentres)) {
    if (length(extracentres) == 0 | mode(extracentres) != "numeric") 
      stop("extracentres must be a non-empty numeric vector")

    if (any(is.infinite(extracentres)))
      stop("infinite cases in extracentres must be removed")

    if (any(is.na(extracentres)))
      warning("missing values in extracentres have been removed")

    extracentres = extracentres[!is.na(extracentres)]
    kerncentres = c(x, extracentres)
  } else {
    kerncentres = x
  }

  if (any(kerncentres < 0))
    stop("kernel centres cannot be non-positive")

  if (is.null(lambda))
    stop("bandwidth (lambda) must be specified")

  if (length(lambda) != 1)
    stop("bandwidth must be scalar")

  if ((mode(lambda) != "numeric") & !is.finite(lambda))
    stop("bandwidth must be numeric")

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  allmethods = c("simple", "renorm", "reflect", "logtrans", 
    "beta1", "beta2", "gamma1", "gamma2", "copula")
  if (length(bcmethod) != 1)
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if ((mode(bcmethod) != "character"))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))
  
  if (!(bcmethod %in% allmethods))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if (!is.logical(proper))
    stop("proper must be logical")
  
  if (length(proper) != 1)
    stop("proper must be of length 1")

  allnn = c("none", "zero", "jf96")
  if (length(nn) != 1)
    stop(paste("nn must be one of", allnn, collapse = " "))

  if ((mode(nn) != "character"))
    stop(paste("nn must be one of", allnn, collapse = " "))
  
  if (!(nn %in% allnn))
    stop(paste("nn must be one of", allnn, collapse = " "))
  
  if (length(offset) != 1)
    stop("offset must be scalar")

  if ((mode(offset) != "numeric") & !is.finite(offset))
    stop("offset must be numeric")
  
  if (offset < 0)
      stop("offset must be non-negative")  

  if ((offset != 0) & (bcmethod != "logtrans"))
    warning("offset only relevant for logtrans method")
    
  if (length(xmax) != 1)
    stop("xmax must be scalar")

  if ((mode(xmax) != "numeric"))
    stop("xmax must be numeric")
  
  if (xmax <= 0)
      stop("xmax must be positive")  

  upboundmethods = c("beta1", "beta2", "copula")
  if ((!is.infinite(xmax)) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for methods", upboundmethods, collapse = " "))
    
  if ((is.infinite(xmax)) & (bcmethod %in% upboundmethods))
    warning(paste("xmax cannot be infinite for methods", upboundmethods, collapse = " "))

  if ((bcmethod %in% upboundmethods) & (max(kerncentres) > xmax))
    stop("xmax must be higher than largest kernel centre as it is defined as upper bound")
  
  dbckdeni <- function(i, kerncentres, lambda, bcmethod, proper, nn, offset, xmax) {
    di = dbckden(kerncentres[i], kerncentres = kerncentres[-i], lambda = lambda, 
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax, log = TRUE)
  }
  
  if ((lambda <= 0) | ((bcmethod == "copula") & (lambda >= 1)) |
      ((bcmethod == "beta1") & (lambda >= 0.25*xmax)) | 
      ((bcmethod == "beta2") & (lambda >= 0.25*xmax))) {
    l = -Inf
  } else {
    l = sum(sapply(1:n, FUN = dbckdeni, kerncentres = kerncentres, lambda = lambda,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax))
  }
  
  if (!log) l = exp(l)

  l
}
