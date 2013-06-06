#' @export
#' 
#' @title MLE Fitting of Dynamically Weighted Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the dynamically weighted mixture model
#'
#' @param pvector vector of initial values of mixture model parameters (\code{wshape}, \code{wscale}, \code{cmu}, \code{ctau}, \code{sigmau}, \code{xi}) or \code{NULL}
#' @inheritParams fnormgpd
#' 
#' @details The dynamically weighted mixture model is fitted to the entire dataset using maximum 
#' likelihood estimation. The estimated parameters, variance-covariance matrix and their standard
#' errors are automatically output.
#' 
#' Non-positive data are ignored.
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
#' \code{call}: \tab \code{optim} call\cr
#' \code{x}: \tab data vector \code{x}\cr
#' \code{init}: \tab \code{pvector}\cr
#' \code{optim}: \tab complete \code{optim} output\cr
#' \code{mle}: \tab vector of MLE of model parameters\cr
#' \code{cov}: \tab variance-covariance matrix of MLE of model parameters\cr
#' \code{se}: \tab vector of standard errors of MLE of model parameters\cr
#' \code{rate}: \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#' \code{nllh}: \tab minimum negative log-likelihood\cr
#' \code{allparams}: \tab vector of MLE of model parameters and \code{phiu}\cr
#' \code{allse}: \tab vector of standard error of all parameters and \code{phiu}\cr
#' \code{n}: \tab total sample size\cr
#' \code{wshape}: \tab MLE of Weibull shape\cr
#' \code{wscale}: \tab MLE of Weibull scale\cr
#' \code{mu}: \tab MLE of Cauchy location\cr
#' \code{tau}: \tab MLE of Cauchy scale\cr
#' \code{sigmau}: \tab MLE of GPD scale\cr
#' \code{xi}: \tab MLE of GPD shape\cr
#' \code{phiu}: \tab MLE of tail fraction\cr
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and to make it 
#' as useable as possible.
#'  
#' @note Unlike all the distribution functions for the extreme value mixture models,
#' the MLE fitting only permits single scalar values for each parameter and 
#' \code{phiu}. Only the data is a vector.
#' 
#' When \code{pvector=NULL} then the initial values are calculated, type 
#' \code{fdwm} to see the default formulae used. The mixture model fitting can be
#' ***extremely*** sensitive to the initial values, so you if you get a poor fit then
#' try some alternatives. Avoid setting the starting value for the shape parameter to
#' \code{xi=0} as depending on the optimisation method it may be get stuck.
#' 
#' The fitting function will stop if infinite sample values are given.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Weibull_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Frigessi, A., O. Haug, and H. Rue (2002). A dynamic mixture model for unsupervised tail
#' estimation without threshold selection. Extremes 5 (3), 219-235
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @family dwm
#' 
#' @examples
#' \dontrun{
#' x = rweibull(1000, shape = 2)
#' xx = seq(-1, 4, 0.01)
#' y = dweibull(xx, shape = 2)
#' 
#' fit = fdwm(x, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 4))
#' lines(xx, y)
#' lines(xx, ddwm(xx, wshape = fit$wshape, wscale = fit$wscale, cmu = fit$cmu, ctau = fit$ctau,
#'   sigmau = fit$sigmau, xi = fit$xi), col="red")
#' }

# maximum likelihood fitting for weibull bulk with GPD for upper tail
fdwm = function(x, pvector = NULL, std.err = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {
  
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
  
  if (any(x < 0))
    warning("negative values have been removed")
  
  x = x[x >= 0]
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")
  
  if ((method == "L-BFGS-B") | (method == "BFGS"))
    finitelik = TRUE
  
  if (is.null(pvector)){
    initfweibull = fitdistr(x, "weibull")
    pvector[1] = initfweibull$estimate[1]
    pvector[2] = initfweibull$estimate[2]    
    pvector[3] = quantile(x, 0.7)
    pvector[4] = sd(x)/10
    pvector[5] = sqrt(6*var(x))/pi
    pvector[6] = 0.1
  } else {
    if (length(pvector) != 6)
      stop("Initial values for six parameters must be specified")
    if (any(!is.finite(pvector)) | is.logical(pvector))
      stop("initial parameters must be numeric")
  }
  
  nllh = nldwm(pvector, x = x, finitelik = finitelik)
  if (is.infinite(nllh))
    stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nldwm, x = x, finitelik = finitelik,
    method = method, control = control, hessian = TRUE, ...)
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
  n = length(x)
  wshape = fit$par[1]
  wscale = fit$par[2]
  cmu = fit$par[3]
  ctau = fit$par[4]
  sigmau = fit$par[5]
  xi = fit$par[6]
  
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
    conv = conv, cov = invhess, mle = fit$par, se = se, nllh = fit$value,
    allparam = fit$par, allse = se, n = n,
    wshape = wshape, wscale = wscale, cmu = cmu, ctau = ctau, sigmau = sigmau, xi = xi)
}
