#' @export
#' 
#' @title Cross-validation MLE Fitting of Kernel Density Estimation for Bulk and GPD
#' for Both Upper and Lower Tails in Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for the extreme value mixture model
#' with kernel density estimation using normal kernel for bulk distribution between
#' the upper and lower thresholds with conditional GPD's for the two tails
#'
#' @inheritParams fkdengpd
#' @inheritParams fgng
#' 
#' @details Extreme value mixture model combining for the extreme value mixture model
#' with kernel density estimation using normal kernel for bulk distribution between the
#' upper and lower thresholds with conditional GPD's for the two tails. Fitted to the
#' entire dataset using maximum cross-validation likelihood
#'   estimation. The estimated parameters, their variance and standard error are
#'   automatically output.
#'   
#'   Cross-validation likelihood is used for kernel density component, but 
#'   standard likelihood is used for GPD components. The default value for
#'   \code{phiul=TRUE} so that the tail fraction is specified by normal
#'   distribution \eqn{\phi_u = 1 - mean(H(ul))}. When \code{phiul=FALSE} then the tail 
#'   fraction is treated as an extra parameter estimated using the MLE which is
#'   the sample proportion below the threshold \code{ul}. In this case the standard error
#'   for \code{phiu} is estimated and output as \code{sephiu}.
#'   
#'   Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so
#'   are ignored, which is inconsistent with the \code{\link[evd:fpot]{evd}}
#'   library which assumes the missing values are below the threshold.
#'   
#'   The default optimisation algorithm is "BFGS", which requires a finite
#'   negative log-likelihood function evaluation \code{finitelik=TRUE}. For
#'   invalid parameters, a zero likelihood is replaced with \code{exp(-1e6)}.
#'   The "BFGS" optimisation algorithms require finite values for likelihood, so
#'   any user input for \code{finitelik} will be overridden and set to
#'   \code{finitelik=TRUE} if either of these optimisation methods is chosen.
#'   
#'   It will display a warning for non-zero convergence result comes from 
#'   \code{\link[stats:optim]{optim}} function call.
#'   
#'   If the hessian is of reduced rank then the variance (from inverse hessian) 
#'   and standard error of bandwidth parameter cannot be calculated, then by
#'   default \code{std.err=TRUE} and the function will stop. If you want the
#'   bandwidth estimate even if the hessian is of reduced rank (e.g. in a
#'   simulation study) then set \code{std.err=FALSE}.
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
#' @return Returns a simple list with the following elements
#'
#' \tabular{ll}{
#' \code{call}: \tab \code{optim} call\cr
#' \code{x}: \tab (jittered) data vector \code{x}\cr
#' \code{kerncentres}: actual kernel centres used \code{x}\cr
#' \code{init}: \tab \code{pvector}\cr
#' \code{optim}: \tab complete \code{optim} output\cr
#' \code{mle}: \tab vector of MLE of parameters\cr
#' \code{cov}: \tab variance of MLE parameters\cr
#' \code{se}: \tab standard error of MLE parameters\cr
#' \code{nllh}: \tab minimum negative cross-validation log-likelihood\cr
#' \code{allparams}: \tab vector of MLE of model parameters, including \code{phiul} and \code{phiur}\cr
#' \code{allse}: \tab vector of standard error of all parameters, including \code{phiul} and \code{phiur}\cr
#' \code{n}: \tab total sample size\cr
#' \code{lambda}: \tab MLE of bandwidth\cr
#' \code{ul}: \tab lower threshold\cr
#' \code{sigmaul}: \tab MLE of lower tail GPD scale\cr
#' \code{xil}: \tab MLE of lower tail GPD shape\cr
#' \code{phiul}: \tab MLE of lower tail fraction\cr
#' \code{ur}: \tab upper threshold\cr
#' \code{sigmaur}: \tab MLE of upper tail GPD scale\cr
#' \code{xir}: \tab MLE of upper tail GPD shape\cr
#' \code{phiur}: \tab MLE of upper tail fraction\cr
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and to make it 
#' as useable as possible.
#'  
#' @note When \code{pvector=NULL} then the initial value for the parameters are calculated 
#' type \code{fkdengpdcon} to see how.
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
#' @seealso \code{\link[evmix:fkdengpd]{fkdengpd}}, \code{\link[evmix:fkden]{fkden}}, 
#' \code{\link[base:jitter]{jitter}},
#' \code{\link[stats:density]{density}} and \code{\link[stats:bandwidth]{bw.nrd0}}
#' 
#' @family gkg
#' 
#' @examples
#' \dontrun{
#' x = rnorm(1000, 0, 1)
#' fit = fgkg(x, phiul = FALSE, phiur = FALSE, std.err = FALSE)
#' hist(x, 100, freq = FALSE, xlim = c(-5, 5)) 
#' xx = seq(-5, 5, 0.01)
#' lines(xx, dgkg(xx, x, fit$lambda, fit$ul, fit$sigmaul, fit$xil, fit$phiul,
#'  fit$ur, fit$sigmaur, fit$xir, fit$phiur), col="blue")
#' abline(v = fit$ul)
#' abline(v = fit$ur)
#' }

# maximum cross-validation likelihood fitting for KDE for bulk with GPD's for both upper and lower tails
fgkg<- function(x, phiul = TRUE, phiur = TRUE, pvector = NULL,
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
  
  if (is.null(pvector)){
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10){
      stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    pvector[1] = bw.nrd0(x)
    pvector[2] = as.vector(quantile(x, 0.1))
    initfgpd = fgpd(-x, -pvector[2], std.err = std.err)
    pvector[3] = initfgpd$sigmau
    pvector[4] = initfgpd$xi 
    pvector[5] = as.vector(quantile(x, 0.9))
    initfgpd = fgpd(x, pvector[5], std.err = std.err)
    pvector[6] = initfgpd$sigmau
    pvector[7] = initfgpd$xi
  }

  if (length(pvector) != 7)
    stop("Initial values for seven parameters must be specified")
  
  if (any(!is.finite(pvector)) | is.logical(pvector))
    stop("initial parameters must be numeric")
  
  nllh = nlgkg(pvector, x = x, phiul = phiul, phiur = phiur, finitelik = finitelik)
  if (is.infinite(nllh))
    stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nlgkg, x = x, 
    phiul = phiul, phiur = phiur, finitelik = finitelik,
    method = method, control = control, hessian = TRUE, ...)
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
  n = length(x)
  lambda = fit$par[1]
  ul = fit$par[2]
  sigmaul = fit$par[3]
  xil = fit$par[4]
  ur = fit$par[5]
  sigmaur = fit$par[6]
  xir = fit$par[7]
  
  if (phiul) {
    phiul = pkdenx(ul, x, lambda)
    sephiul = NA
  } else {
    phiul = mean(x < ul, na.rm = TRUE)
    sephiul = sqrt(phiul * (1 - phiul) / n)
  }
  if (phiur) {
    phiur = 1 - pkdenx(ur, x, lambda)
    sephiur = NA
  } else {
    phiur = mean(x > ur, na.rm = TRUE)
    sephiur = sqrt(phiur * (1 - phiur) / n)
  }
  phib = (1 - phiul - phiur) / (pkdenx(ur, x, lambda) - pkdenx(ul, x, lambda))
  
  if (conv & std.err) {
    qrhess = qr(fit$hessian)
    if (qrhess$rank != ncol(qrhess$qr)) {
      warning("observed information matrix is singular; use std.err = FALSE")
      se = NULL
      invhess = NULL
    } else {
      invhess = solve(qrhess)
      vars = diag(invhess)
      if (any(vars <= 0)) {
        warning("observed information matrix is singular; use std.err = FALSE")
        invhess = NULL
        se = NULL
      } else {
        se = sqrt(vars)
      }  
    }
  } else {
    invhess = NULL
    se = NULL
  }
  
  list(call = call, x = as.vector(x), kerncentres = x,
    init = as.vector(pvector), optim = fit,
    conv = conv, cov = invhess, mle = fit$par, se = se, ratel = phiul, rater = phiur, 
    nllh = fit$value, allparam = c(fit$par, phiul, phiur), allse = c(se, sephiul, sephiur), n = n,
    lambda = lambda, ul = ul, sigmaul = sigmaul, xil = xil, phiul = phiul,
    ur = ur, sigmaur = sigmaur, xir = xir, phiur = phiur)
}
