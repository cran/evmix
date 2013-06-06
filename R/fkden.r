#' @export
#' 
#' @title Cross-validation MLE Fitting of Kernel Density Estimator Using normal Kernel
#'
#' @description Maximum likelihood estimation for fitting kernel density estimator
#' using a normal kernels, by treating it as a mixture model.
#'
#' @inheritParams lkden
#' @inheritParams fnormgpd
#' @param linit        initial value for bandwidth parameter or \code{NULL}
#' @param extracentres extra kernel centres used in KDE, but likelihood contribution not evaluated, or \code{NULL}
#' @param add.jitter   logical, whether jitter is needed for rounded data
#' @param factor       see \code{\link[base:jitter]{jitter}}
#' @param amount       see \code{\link[base:jitter]{jitter}}
#' 
#' @details The kernel density estimator with a normal kernel is fitted to the entire
#' dataset using maximum cross-validation likelihood estimation. The estimated
#' bandwidth, variance and standard error are automatically output.
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
#' @family kden
#' 
#' @examples
#' \dontrun{
#' nk=50
#' x = rnorm(nk)
#' xx = seq(-5, 5, 0.01)
#' fit = fkden(x)
#' hist(x, nk/5, freq = FALSE, xlim = c(-5, 7)) 
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit$lambda)*0.05)
#' lines(xx,dnorm(xx), col = "black")
#' lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' lines(density(x, bw = fit$lambda), lwd = 2, lty = 2,  col = "blue")
#' legend("topright", c("True Density", "KDE fitted evmix",
#' "KDE Using density, default bandwidth", "KDE Using density, c-v likelihood bandwidth"),
#' lty = c(1, 1, 2, 2), lwd = c(1, 2, 2, 2), col = c("black", "red", "green", "blue"))
#' 
#' # bandwidth is biased towards oversmoothing for heavy tails
#' nk=100
#' x = rt(nk, df = 2)
#' xx = seq(-8, 8, 0.01)
#' fit = fkden(x)
#' hist(x, seq(floor(min(x)), ceiling(max(x)), 0.5), freq = FALSE, xlim = c(-8, 10)) 
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit$lambda)*0.05)
#' lines(xx,dt(xx , df = 2), col = "black")
#' lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "red")
#' legend("topright", c("True Density", "KDE fitted evmix, c-v likelihood bandwidth"),
#' lty = c(1, 1), lwd = c(1, 2), col = c("black", "red"))
#' 
#' # remove heavy tails from likelihood evaluation, but still include used in KDE within likelihood
#' # often gives better bandwidth
#' nk=100
#' x = rt(nk, df = 2)
#' xx = seq(-8, 8, 0.01)
#' fit2 = fkden(x[(x > -4) & (x < 4)], extracentres = x[(x <= -4) | (x >= 4)])
#' hist(x, seq(floor(min(x)), ceiling(max(x)), 0.5), freq = FALSE, xlim = c(-8, 10)) 
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit2$lambda)*0.05)
#' lines(xx,dt(xx , df = 2), col = "black")
#' lines(xx, dkden(xx, x, lambda = fit2$lambda), lwd = 2, col = "red")
#' lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "blue")
#' legend("topright", c("True Density", "KDE fitted evmix, tails removed",
#' "KDE fitted evmix, tails included"),
#' lty = c(1, 1, 1), lwd = c(1, 2, 2), col = c("black", "red", "blue"))
#' }

# maximum cross-validation likelihood fitting for kernel density estimator using normal kernel
fkden <- function(x, linit = NULL, extracentres = NULL, add.jitter = FALSE, factor = 0.1, amount = NULL,
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
  
  if (is.null(linit)){
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10){
        stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    linit = bw.nrd0(x)
  }

  if (length(linit) != 1)
    stop("initial bandwidth must be scalar")

  if ((mode(linit) != "numeric") & !is.finite(linit))
    stop("initial bandwidth must be numeric")

  nllh = nlkden(linit, x = x, extracentres = extracentres, finitelik = finitelik)
  if (is.infinite(nllh))
    stop("initial bandwidth is invalid")

  fit = optim(par = linit, fn = nlkden, x = x, extracentres = extracentres, finitelik = finitelik,
    method = method, control = control, hessian = TRUE, ...)

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
    se = se, nllh = fit$value, n = n, lambda = lambda)
}
