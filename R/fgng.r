#' @export
#' 
#' @title MLE Fitting of Normal Bulk with GPD Upper and Lower Tails Extreme Value Mixture Model
#' 
#' @description Maximum likelihood estimation for the extreme value 
#' mixture model with normal for bulk distribution between the lower and upper
#' thresholds with conditional GPD for the two tails.
#'
#' @param pvector vector of initial values of mixture model parameters or \code{NULL}
#' @param phiul    logical
#' @param phiur    logical
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with normal bulk and GPD for both tails is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' \code{pvector} is (\code{nmean}, \code{nsd}, \code{ul}, \code{sigmaul}, \code{xil},
#' \code{ur}, \code{sigmaur}, \code{xir})
#' 
#' The default values for \code{phiul=TRUE} and \code{phiur=TRUE} so that the 
#' corresponding tail fractions are specified by normal distribution 
#' \eqn{\phi_{ul} = H(u_l)} and \eqn{\phi_{ur} = 1 - H(u_r)}. When \code{phiul=FALSE}
#' and \code{phiur=FALSE} then the corresponding tail fractions are treated as an
#' extra parameter estimated using the MLE which is the
#' sample proportion beyond the corresponding threshold. In this case the standard error for 
#' \code{phiul} and \code{phiur} are estimated and output as \code{sephiul} and 
#' \code{sephiur}.
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored,
#' which is inconsistent with the \code{\link[evd:fpot]{evd}} library which assumes the 
#' missing values are below the threshold.
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
#' If the hessian is of reduced rank then the variance covariance (from inverse hessian)
#' and standard error of parameters cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the parameter estimates
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' @return Returns a simple list with the following elements
#'
#' \tabular{ll}{
#' \code{x}: \tab data vector \code{x}\cr
#' \code{init}: \tab \code{pvector}\cr
#' \code{optim}: \tab complete \code{optim} output\cr
#' \code{mle}: \tab vector of MLE of model parameters\cr
#' \code{cov}: \tab variance-covariance matrix of MLE of model parameters\cr
#' \code{se}: \tab vector of standard errors of MLE of model parameters\cr
#' \code{rate}: \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#' \code{nllh}: \tab minimum negative log-likelihood\cr
#' \code{allparams}: \tab vector of MLE of model parameters and tail fractions \code{phiul} and \code{phiur}\cr
#' \code{allse}: \tab vector of standard error of model parameters and tail fractions \code{phiul} and \code{phiur}\cr
#' \code{n}: \tab total sample size\cr
#' \code{nmean}: \tab MLE of normal mean\cr
#' \code{nsd}: \tab MLE of normal standard deviation\cr
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
#' @note Unlike all the distribution functions for the extreme value mixture models,
#' the MLE fitting only permits single scalar values for each parameter and 
#' tail fractions \code{phiul} and \code{phiur}. Only the data is a vector.
#' 
#' When \code{pvector=NULL} then the initial values are calculated, type 
#' \code{fgng} to see the default formulae used. The mixture model fitting can be
#' ***extremely*** sensitive to the initial values, so you if you get a poor fit then
#' try some alternatives. Avoid setting the starting value for the shape parameters to
#' \code{xil=0} or \code{xir=0} as depending on the optimisation method it may be get stuck.
#' 
#' If the hessian is of reduced rank then the variance covariance (from inverse hessian)
#' and standard error of parameters cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the parameter estimates
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' The fitting function will stop if infinite sample values are given.
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
#' @seealso \code{\link[evmix:fnormgpd]{fnormgpd}}, \code{\link[evmix:lgpd]{lgpd}}
#' and \code{\link[evmix:gpd]{gpd}}
#' @family gng
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2))
#' x = rnorm(1000)
#' xx = seq(-6, 6, 0.01)
#' y = dnorm(xx)
#' 
#' # Bulk model base tail fraction
#' fit = fgng(x, phiul = TRUE, phiur = TRUE, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6), main = "N(0, 1)")
#' lines(xx, y)
#' lines(xx, dgng(xx, nmean = fit$nmean, nsd = fit$nsd,
#'   ul = fit$ul, sigmaul = fit$sigmaul, xil = fit$xil, phiul = TRUE,
#'   ur = fit$ur, sigmaur = fit$sigmaur, xir = fit$xir, phiur = TRUE), col="red")
#' abline(v = c(fit$ul, fit$ur))
#'   
#' # Parameterised tail fraction
#' fit2 = fgng(x, phiul = TRUE, phiur = TRUE, std.err = FALSE)
#' plot(xx, y, type = "l")
#' lines(xx, dgng(xx, nmean = fit$nmean, nsd = fit$nsd,
#'   ul = fit$ul, sigmaul = fit$sigmaul, xil = fit$xil, phiul = TRUE,
#'   ur = fit$ur, sigmaur = fit$sigmaur, xir = fit$xir, phiur = TRUE), col="red")
#' lines(xx, dgng(xx, nmean = fit2$nmean, nsd = fit2$nsd,
#'   ul = fit2$ul, sigmaul = fit2$sigmaul, xil = fit2$xil, phiul = fit2$phiul,
#'   ur = fit2$ur, sigmaur = fit2$sigmaur, xir = fit2$xir, phiur = fit2$phiur), col="blue")
#' abline(v = c(fit$ul, fit$ur), col = "red")
#' abline(v = c(fit2$ul, fit2$ur), col = "blue")
#' legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#' x = rnorm(1000)
#' xx = seq(-6, 6, 0.01)
#' y = dnorm(xx)
#' 
#' # Two tail is safest if bulk has lower tail which is not normal tail
#' x = rt(1000, df = 3)
#' xx = seq(-10, 10, 0.01)
#' y = dt(xx, df = 3)
#' 
#' # Bulk model base tail fraction
#' fit = fnormgpd(x, phiu = FALSE, std.err = FALSE)
#' fit2 = fgng(x, phiul = FALSE, phiur = FALSE, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-10, 10), main = "t (df=3)")
#' lines(xx, y)
#' lines(xx, dnormgpd(xx, nmean = fit$nmean, nsd = fit$nsd,
#'   u = fit$u, sigmau = fit$sigmau, xi = fit$xi, phiu = fit$phiu), col="red")
#' abline(v = fit$u)
#'   
#' # Bulk model base tail fraction
#' plot(xx, y, type = "l")
#' lines(xx, dnormgpd(xx, nmean = fit$nmean, nsd = fit$nsd,
#'   u = fit$u, sigmau = fit$sigmau, xi = fit$xi, phiu = fit$phiu), col="red")
#' lines(xx, dgng(xx, nmean = fit2$nmean, nsd = fit2$nsd,
#'   ul = fit2$ul, sigmaul = fit2$sigmaul, xil = fit2$xil, phiul = fit2$phiul,
#'   ur = fit2$ur, sigmaur = fit2$sigmaur, xir = fit2$xir, phiur = fit2$phiur), col="blue")
#' abline(v = c(fit$ul, fit$ur), col = "red")
#' abline(v = c(fit2$ul, fit2$ur), col = "blue")
#' legend("topright", c("True Density","GPD Upper Tail Only","GPD Both Tails"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   }

# maximum likelihood fitting for normal bulk with GPD's for upper and lower tails
fgng <- function(x, phiul = TRUE, phiur = TRUE, pvector = NULL, std.err = TRUE,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

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
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  if (!is.logical(phiul) | !is.logical(phiur))
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or FALSE when using parameterised threshold probability approach")

  if ((method == "L-BFGS-B") | (method == "BFGS"))
    finitelik = TRUE

  if (is.null(pvector)) {
    pvector[1] = mean(x, trim = 0.2)
    pvector[2] = sd(x)
    pvector[3] = as.vector(quantile(x, 0.1))
    initfgpd = fgpd(-x, -pvector[3], std.err = std.err)
    pvector[4] = initfgpd$sigmau
    pvector[5] = initfgpd$xi
    pvector[6] = as.vector(quantile(x, 0.9))
    initfgpd = fgpd(x, pvector[6], std.err = std.err)
    pvector[7] = initfgpd$sigmau
    pvector[8] = initfgpd$xi
  } else {
    if (length(pvector) != 8)
      stop("Initial values for eight parameters must be specified")
    if (any(!is.finite(pvector)) | is.logical(pvector))
      stop("initial parameters must be numeric")
  }

  nllh = nlgng(pvector, x = x, phiul = phiul, phiur = phiur, finitelik = finitelik)
  if (is.infinite(nllh))
    stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nlgng, x = x, phiul = phiul, phiur = phiur, finitelik = finitelik,
    method = method, control = control, hessian = TRUE, ...)

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  n = length(x)
  nmean = fit$par[1]
  nsd = fit$par[2]
  ul = fit$par[3]
  sigmaul = fit$par[4]
  xil = fit$par[5]
  ur = fit$par[6]
  sigmaur = fit$par[7]
  xir = fit$par[8]
  
  if (phiul) {
    phiul = pnorm(ul, mean = nmean, sd = nsd)
    sephiul = NA
  } else {
    phiul = mean(x < ul, na.rm = TRUE)
    sephiul = sqrt(phiul * (1 - phiul) / n)
  }
  if (phiur) {
    phiur = 1 - pnorm(ur, mean = nmean, sd = nsd)
    sephiur = NA
  } else {
    phiur = mean(x > ur, na.rm = TRUE)
    sephiur = sqrt(phiur * (1 - phiur) / n)
  }  

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

  list(call = call, x = as.vector(x), init = as.vector(pvector), optim = fit,
    conv = conv, cov = invhess, mle = fit$par, se = se,
    ratel = phiul, rater = phiur, nllh = fit$value,
    allparam = c(fit$par, phiul, phiur), allse = c(se, sephiul, sephiur),
    n = n, nmean = nmean, nsd = nsd,
    ul = ul, sigmaul = sigmaul, xil = xil, phiul = phiul,
    ur = ur, sigmaur = sigmaur, xir = xir, phiur = phiur)
}
