#' @export
#' 
#' @title MLE Fitting of Normal Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with normal for bulk distribution upto the threshold and conditional
#' GPD above threshold
#'
#' @param pvector vector of initial values of mixture model parameters (\code{nmean}, \code{nsd}, \code{u}, \code{sigmau}, \code{xi}) or \code{NULL}
#' @param phiu    logical
#' @param control    optimisation control list (see \code{\link[stats:optim]{optim}})
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with normal bulk and GPD tail is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' The default value for \code{phiu=TRUE} so that the tail fraction is specified by
#' normal distribution \eqn{\phi_u = 1 - H(u)}. When \code{phiu=FALSE} then the tail
#' fraction is treated as an extra parameter estimated using the MLE which is the
#' sample proportion above the threshold. In this case the standard error for 
#' \code{phiu} is estimated and output as \code{sephiu}.
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
#' \code{nmean}: \tab MLE of normal mean\cr
#' \code{nsd}: \tab MLE of normal standard deviation\cr
#' \code{u}: \tab threshold\cr
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
#' \code{fnormgpd} to see the default formulae used. The mixture model fitting can be
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
#' \url{http://en.wikipedia.org/wiki/Normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004). Bayesian analysis of extreme
#' events with threshold estimation. Statistical Modelling. 4(3), 227-244.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:lgpd]{lgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @family   normgpd
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(2,1))
#' x = rnorm(1000)
#' xx = seq(-4, 4, 0.01)
#' y = dnorm(xx)
#' 
#' # Bulk model base tail fraction
#' fit = fnormgpd(x, phiu = TRUE, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' lines(xx, dnormgpd(xx, nmean = fit$nmean, nsd = fit$nsd, u = fit$u,
#'   sigmau = fit$sigmau, xi = fit$xi, phiu = TRUE), col="red")
#' abline(v = fit$u)
#'   
#' # Parameterised tail fraction
#' fit2 = fnormgpd(x, phiu = FALSE, std.err = FALSE)
#' plot(xx, y, type = "l")
#' lines(xx, dnormgpd(xx, nmean = fit$nmean, nsd = fit$nsd, u = fit$u,
#'   sigmau = fit$sigmau, xi = fit$xi, phiu = TRUE), col="red")
#' lines(xx, dnormgpd(xx, nmean = fit2$nmean, nsd = fit2$nsd, u = fit2$u,
#'   sigmau = fit2$sigmau, xi = fit2$xi, phiu = fit2$phiu), col="blue")
#' abline(v = fit$u, col = "red")
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   }

# maximum likelihood fitting for normal bulk with GPD for upper tail
fnormgpd <- function(x, phiu = TRUE, pvector = NULL, std.err = TRUE, method = "BFGS",
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
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  if (!is.logical(phiu))
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or FALSE when using parameterised threshold probability approach")

  if ((method == "L-BFGS-B") | (method == "BFGS"))
    finitelik = TRUE
  
  if (is.null(pvector)) {
    pvector[1] = mean(x)
    pvector[2] = sd(x)
    pvector[3] = as.vector(quantile(x, 0.9))
    initfgpd = fgpd(x, pvector[3], std.err = std.err)
    pvector[4] = initfgpd$sigmau
    pvector[5] = initfgpd$xi
  } else {
    if (length(pvector) != 5)
      stop("Initial values for five parameters must be specified")
    if (any(!is.finite(pvector)) | is.logical(pvector))
      stop("initial parameters must be numeric")
  }

  nllh = nlnormgpd(pvector, x = x, phiu = phiu, finitelik = finitelik)
  if (is.infinite(nllh))
    stop("initial parameter values are invalid")
  
  fit = optim(par = as.vector(pvector), fn = nlnormgpd, x = x, phiu = phiu, finitelik = finitelik,
    method = method, control = control, hessian = TRUE, ...)

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  n = length(x)
  nmean = fit$par[1]
  nsd = fit$par[2]
  u = fit$par[3]
  sigmau = fit$par[4]
  xi = fit$par[5]
  
  if (phiu) {
    phiu = 1 - pnorm(u, mean = nmean, sd = nsd)
    sephiu = NA
  } else {
    phiu = mean(x > u, na.rm = TRUE)
    sephiu = sqrt(phiu * (1 - phiu) / n)
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
    conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu, nllh = fit$value,
    allparam = c(fit$par, phiu), allse = c(se, sephiu), n = n,
    nmean = nmean, nsd = nsd, u = u, sigmau = sigmau, xi = xi, phiu = phiu)
}
