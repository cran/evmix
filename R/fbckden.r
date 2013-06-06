#' @export
#' 
#' @title Cross-validation MLE Fitting of Boundary Corrected Kernel Density Estimation
#'
#' @description Maximum likelihood estimation for fitting boundary corrected 
#' kernel density estimator, by treating it as a mixture model.
#'
#' @inheritParams lbckden
#' @inheritParams fnormgpd
#' @inheritParams fkden
#' @inheritParams fgpd
#' 
#' @details The boundary corrected kernel density estimator is fitted to the entire
#' dataset using maximum cross-validation likelihood estimation. The estimated
#' bandwidth, variance and standard error are automatically output.
#' 
#' The \code{beta1} and \code{beta2} densities requires renormalisation which is achieved
#' by numerical integration, so is very time consuming. Practically we have found leaving out
#' the renormalisation \code{proper=FALSE} still yields reliable bandwidth estimates.
#' 
#' The cross-validation likelihood estimates of the bandwidth for the \code{simple}, \code{gamma1}
#' and \code{gamma2} methods of boundary correction are biased high (particularly when there is a
#' pole at the boundary) leading to oversmoothing. We have empirically found that leaving off the
#' data from the upper tail in the likelihood appears to help, see examples for an implementation.
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored.
#' 
#' Normally for likelihood estimation of the bandwidth the kernel centres and the data 
#' where the likelihood is evaluated are the same. However, when using KDE for extreme value
#' mixture modelling the likelihood only those data in the bulk of the distribution should contribute
#' to the likelihood, but all the data (including those beyond the threshold) should contribute to
#' the density estimate. The \code{extracentres} option allows the use to specify extra kernel
#' centres used in estimating the density, but not evaluated in the likelihood. The default
#' is to just use the existing data, so \code{extracentres=NULL}.
#' 
#' The default optimisation algorithm is "BFGS", which requires a finite negative 
#' log-likelihood function evaluation \code{finitelik=TRUE}. For invalid 
#' parameters, a zero likelihood is replaced with \code{exp(-1e6)}. The "BFGS" 
#' optimisation algorithms require finite values for likelihood, so any user 
#' input for \code{finitelik} will be overridden and set to \code{finitelik=TRUE} 
#' if either of these optimisation methods is chosen.
#' 
#' It will display a warning for non-zero convergence result comes from 
#' \code{\link[stats:optim]{optim}} function call.
#' 
#' If the hessian is of reduced rank then the variance (from inverse hessian)
#' and standard error of bandwidth parameter cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the bandwidth estimate
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' @section Warning:
#' Two important practical issues arise with MLE for the kernel bandwidth:
#' 1) Cross-validation likelihood is needed for the KDE bandwidth parameter
#' as the usual likelihood degenerates, so that the MLE \eqn{\hat{\lambda} \rightarrow 0} as
#' \eqn{n \rightarrow \infty}, thus giving a negative bias towards a small bandwidth.
#' Leave one out cross-validation essentially ensures that some smoothing between the kernel centres
#' is required (i.e. a non-zero bandwidth), otherwise the resultant density estimates would always
#' be zero if the bandwidth was zero.
#' 
#' This problem occassionally rears its ugly head for data which has been heavily rounded,
#' as even when using cross-validation the density can be non-zero even if the bandwidth is zero.
#' To overcome this issue an option to add a small jitter should be added to the data
#' (\code{x} only) has been included in the fitting inputs, using the 
#' \code{\link[base:jitter]{jitter}} function, to remove the ties. The default options red in the 
#' \code{\link[base:jitter]{jitter}} are specified above, but the user can override these.
#' Notice the default scaling \code{factor=0.1}, which is a tenth of the default value in the
#' \code{\link[base:jitter]{jitter}}
#' function itself.
#' 
#' A warning message is given if the data appear to be rounded
#' (i.e. more than 5% of data are tied). If the estimated bandwidth is too small, then
#' data rounding is the likely culprit. Only use the jittering when the MLE of
#' the bandwidth is far too small. 
#' 
#' 2) For heavy tailed populations the bandwidth is positively biased, giving oversmoothing
#' (see example). The bias is due to the distance between the upper (or lower) order statistics not
#' necessarily decaying to zero as the sample size tends to infinity. Essentially, as the distance
#' between the two largest (or smallest) sample datapoints does not decay to zero, some smoothing between
#' them is required (i.e. bandwidth cannot be zero). One solution to this problem is to splice
#' the GPD at a suitable threshold to remove the problematic tail from the inference for the bandwidth, 
#' using either the \code{kdengpd} function for a single heavy tail or the \code{kdengng} function
#' if both tails are heavy. See MacDonald et al (2013).
#' 
#' @return Returns a simple list with the following elements
#'
#' \tabular{ll}{
#' \code{call}: \tab \code{optim} call\cr
#' \code{x}: \tab (jittered) data vector \code{x}\cr
#' \code{kerncentres}: actual kernel centres used \code{x}\cr
#' \code{init}: \tab \code{linit}\cr
#' \code{optim}: \tab complete \code{optim} output\cr
#' \code{mle}: \tab vector of MLE of bandwidth\cr
#' \code{cov}: \tab variance of MLE of bandwidth\cr
#' \code{se}: \tab standard error of MLE of bandwidth\cr
#' \code{nllh}: \tab minimum negative cross-validation log-likelihood\cr
#' \code{n}: \tab total sample size\cr
#' \code{lambda}: \tab MLE of bandwidth\cr
#' \code{bcmethod}: \tab boundary correction method\cr
#' \code{proper}: \tab logical, whether renormalisation is requested\cr
#' \code{nn}: \tab non-negative correction method\cr
#' \code{offset}: \tab offset for log transformation method\cr
#' \code{xmax}: \tab maximum value of scale beta or copula
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and to make it 
#' as useable as possible.
#'  
#' @note When \code{linit=NULL} then the initial value for the bandwidth is calculated 
#' using \code{\link[stats:bandwidth]{bw.nrd0}} function.
#' 
#' The fitting function will stop if infinite sample values are given.
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
#' MacDonald, A., C. J. Scarrott, and D. S. Lee (2011). Boundary correction, consistency
#' and robustness of kernel densities using extreme value theory. Submitted.
#' Available from: \url{http://www.math.canterbury.ac.nz/~c.scarrott}.

#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[base:jitter]{jitter}}, \code{\link[stats:density]{density}} and
#' \code{\link[stats:bandwidth]{bw.nrd0}}
#' 
#' @family bckden fbckdengpd fkden fgpd
#' 
#' @examples
#' \dontrun{
#' nk=50
#' x = rgamma(nk, shape = 3, scale = 1)
#' xx = seq(-1, 10, 0.01)
#' fit = fbckden(x, linit = 0.2, bcmethod = "renorm")
#' hist(x, nk/5, freq = FALSE) 
#' rug(x)
#' for (i in 1:nk) lines(xx, dbckden(xx, x[i], lambda = fit$lambda, bcmethod = "renorm")*0.05)
#' lines(xx, dgamma(xx, shape = 3, scale = 1), col = "black")
#' lines(xx, dbckden(xx, x, lambda = fit$lambda, bcmethod = "renorm"), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "BC KDE fitted evmix",
#' "KDE Using density, default bandwidth"),
#' lty = c(1, 1, 2), lwd = c(1, 2, 2), col = c("black", "red", "green"))
#' 
#' nk=100
#' x = rgamma(nk, shape = 0.5, scale = 1)
#' q75 = qgamma(0.75, shape = 0.5, scale = 1)
#' xx = seq(-1, 10, 0.01)
#' fit = fbckden(x, linit = 0.2, bcmethod = "simple")
#' fitnotail = fbckden(x[x <= q75], linit = 0.1, bcmethod = "simple", extracentres = x[x > q75])
#' hist(x, nk/5, freq = FALSE, ylim = c(0, 8)) 
#' rug(x)
#' lines(xx, dgamma(xx, shape = 0.5, scale = 1), col = "black")
#' lines(xx, dbckden(xx, x, lambda = fit$lambda, bcmethod = "simple"), lwd = 2, col = "red")
#' lines(xx, dbckden(xx, x, lambda = fitnotail$lambda, bcmethod = "simple"), lwd = 2, col = "blue")
#' legend("topright", c("True Density", "BC KDE (complete dataset)",
#' "BC KDE (upper tail ignored)"),
#' lty = c(1, 1, 2), lwd = c(1, 2, 2), col = c("black", "red", "blue"))
#' }

# maximum cross-validation likelihood fitting for boundary corrected KDE
fbckden <- function(x, linit = NULL, extracentres = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf,
  add.jitter = FALSE, factor = 0.1, amount = NULL,
  std.err = TRUE, method = "BFGS",  control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")

  if (any(is.infinite(x)))
    stop("infinite cases must be removed")

  if (any(is.na(x)))
    warning("missing values have been removed")

  x = x[!is.na(x)]
  n = length(x)
  
  if (!is.null(extracentres)) {
    if (length(extracentres) == 0 | mode(extracentres) != "numeric") 
      stop("extracentres must be a non-empty numeric vector")

    if (any(is.infinite(extracentres)))
      stop("infinite cases in extracentres must be removed")

    if (any(is.na(extracentres)))
      warning("missing values in extracentres have been removed")

    extracentres = extracentres[!is.na(extracentres)]
  }

  if (any(extracentres < 0))
    stop("extra kernel centres must be non-negative")
  
  if (any(x < 0))
    stop("data must be non-negative")
  
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

  if ((mode(xmax) != "numeric") & !is.finite(xmax))
    stop("xmax must be numeric")
  
  if (xmax <= 0)
      stop("xmax must be positive")  

  upboundmethods = c("beta1", "beta2", "copula")
  if ((!is.infinite(xmax)) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for methods", upboundmethods, collapse = " "))
    
  if ((is.infinite(xmax)) & (bcmethod %in% upboundmethods))
    warning(paste("xmax cannot be infinite for methods", upboundmethods, collapse = " "))

  if ((bcmethod %in% upboundmethods) & (max(x) > xmax))
    stop("xmax must be higher than largest kernel centre as it is defined as upper bound")

  if (!is.null(extracentres)) {
    if ((bcmethod %in% upboundmethods) & (max(extracentres) > xmax))
      stop("xmax must be higher than largest kernel centre as it is defined as upper bound")
  }
  
  if (!is.logical(add.jitter))
    stop("add.jitter must be logical")
  
  if (length(add.jitter) != 1)
    stop("add.jitter must be of length 1")

  if (length(factor) != 1)
    stop("jitter factor must be (small) numeric value")

  if ((mode(factor) != "numeric") & !is.finite(factor))
    stop("jitter factor must be (small) numeric value")
  
  if (!is.null(amount)) {
    if (length(amount) != 1)
      stop("jitter amount must be NULL or (small) numeric value")

    if ((mode(amount) != "numeric") & !is.finite(amount))
      stop("jitter amount must be NULL or (small) numeric value")
  }

  if (add.jitter) {
    x = jitter(x, factor, amount)
  }
  
  xuniq = unique(x)
  if (length(xuniq) < (0.95*n)) {
    warning("data could be rounded, as there are more than 5% of ties, so bandwidth could be biased towards zero")
  }

  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  if ((method == "L-BFGS-B") | (method == "BFGS"))
    finitelik = TRUE
  
  if (is.null(linit)) {
    stop("initial guess of bandwidth (linit) must be provided")
  }

  if (length(linit) != 1)
    stop("initial bandwidth must be scalar")

  if ((mode(linit) != "numeric") & !is.finite(linit))
    stop("initial bandwidth must be numeric")

  if ((bcmethod == "copula") & (linit >= 1))
      stop("bandwidth initial value must between (0, 1) for copula method")  

  llhinit = lbckden(x, lambda = linit, extracentres = extracentres, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax)
  tryi = 0
  
  lfirst = linit
  
  while (is.infinite(llhinit) & (tryi < 5)) {
    linit = linit*2
    llhinit = lbckden(x, lambda = linit, extracentres = extracentres, bcmethod = bcmethod,
      proper = proper, nn = nn, offset = offset, xmax = xmax)
    tryi = tryi + 1
  }

  if (is.infinite(llhinit)) {
    linit = lfirst/2
    llhinit = lbckden(x, lambda = linit, extracentres = extracentres, bcmethod = bcmethod,
      proper = proper, nn = nn, offset = offset, xmax = xmax)

    tryi = 0
    while (is.infinite(llhinit) & (tryi < 5)) {
      linit = linit/2
      llhinit = lbckden(x, lambda = linit, extracentres = extracentres, bcmethod = bcmethod,
        proper = proper, nn = nn, offset = offset, xmax = xmax)
      tryi = tryi + 1
    }
  }
  
  if (is.infinite(llhinit))
    stop("likelihood is undefined for initial bandwidth try a value at least 16 times bigger")  

  if (tryi != 0)
    warning(paste("initial bandwidth was too small, so linit=", linit, "is used instead"))
  
  fit = optim(par = linit, fn = nlbckden, x = x, extracentres = extracentres, 
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == linit) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  lambda = fit$par[1]
    
  if (conv & std.err) {
    se = sqrt(1/fit$hessian)
  } else {
    se = NULL
  }

  list(call = call, x = as.vector(x), kerncentres = c(x, extracentres), 
    init = as.vector(linit), optim = fit, conv = conv, cov = 1/fit$hessian, mle = fit$par,
    se = se, nllh = fit$value, n = n, lambda = lambda, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax)
}
