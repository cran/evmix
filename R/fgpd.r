#' @export
#' @title MLE Fitting of Generalised Pareto Distribution (GPD)
#'
#' @description Maximum likelihood estimation for fitting the GPD with parameters
#' scale \code{sigmau} and shape \code{xi} to the threshold exceedances, conditional
#' on being above a threshold \code{u}. Unconditional likelihood fitting also
#' provided when the probability \code{phiu} of being above the threshold
#' \code{u} is given.
#'
#' @inheritParams gpd
#' @param x          vector of sample data
#' @param pvector    vector of initial values of GPD parameters (\code{sigmau}, \code{xi}) or \code{NULL}
#' @param phiu       probability of being above threshold [0,1] or \code{NULL}
#' @param std.err    logical, should standard errors be calculated
#' @param method     optimisation method (see \code{\link[stats:optim]{optim}})
#' @param finitelik  logical, should log-likelihood return finite value for invalid parameters
#' @param ...        optional inputs passed to \code{\link[stats:optim]{optim}}
#' 
#' @details The GPD is fitted to the exceedances of the threshold \code{u} using
#' maximum likelihood estimation. The estimated parameters, variance-covariance matrix
#' and their standard errors are automatically output.
#' 
#' The default value for \code{phiu} is NULL, which means it will be estimated as the 
#' MLE of the tail fraction which is the sample proportion above the threshold.
#' In this case the standard error for 
#' \code{phiu} is estimated and output as \code{sephiu}. Consistent with the 
#' \code{\link[evd:fpot]{evd}} library the missing values 
#' (\code{NA} and \code{NaN}) are assumed to be below the threshold. 
#' 
#' Otherwise, \code{phiu} can be specified as any value over [0, 1] leading to
#' the unconditional log-likelihood being used for estimation. In this case the
#' standard error will be output as \code{NA}. The value of \code{phiu} does not
#' effect the GPD parameter estimates, only the value of the likelihood, as:
#' 
#' \deqn{L(\sigma_u, \xi; u, \phi_u) = (\phi_u ^ {n_u}) L(\sigma_u, \xi; u, \phi_u=1)}
#' 
#' where the GPD has scale \eqn{\sigma_u} and shape \eqn{\xi}, the threshold is \eqn{u}
#' and \eqn{nu} is the number of exceedances. A non-unit value for \code{phiu} simply
#' scales the likelihood, and shifts the log-likelihood, thus the GPD parameter
#' estimates are invariant to \code{phiu}.
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
#' \code{allparams}: \tab vector of MLE of GPD parameters and (possibly given) \code{phiu}\cr
#' \code{allse}: \tab vector of standard error of GPD parameters and (possibly given) \code{phiu}\cr
#' \code{n}: \tab total sample size\cr
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
#' @note Unlike all the distribution functions for the GPD, the MLE fitting only permits
#' single scalar values for each parameter, \code{phiu} and threshold \code{u}. Only the
#' data is a vector.
#' 
#' When \code{pvector=NULL} then the initial values are calculated, type \code{fgpd} to
#' see the default formulae used. The GPD fitting is not very sensitive to the
#' initial values, so you will rarely have to  give alternatives. Avoid setting the
#' starting value for the shape parameter to \code{xi=0} as depending on the
#' optimisation method it may be get stuck.
#' 
#' Default values for the threshold \code{u=0} and tail fraction \code{phiu=NULL} are given.
#' If the threshold is left as the default \code{u=0} and the tail fraction set to 
#' \code{phiu=1} then MLE assumes tha excesses over the threshold are given, rather
#' than exceedances.
#' 
#' The fitting function will stop if infinite sample values are given.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' 
#' Based on GPD fitting function in the \code{\link[evd:fpot]{evd}} package.
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evd:gpd]{dgpd}}, \code{\link[evd:fpot]{fpot}}
#' and \code{\link[MASS:fitdistr]{fitdistr}}
#' @family   gpd
#' 
#' @examples
#' par(mfrow=c(2,1))
#' x = rgpd(1000, u = 10, sigmau = 5, xi = 0.2)
#' xx = seq(0, 100, 0.1)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(0, 100))
#' lines(xx, dgpd(xx, u = 10, sigmau = 5, xi = 0.2))
#' fit = fgpd(x, u = 10, phiu = NULL, std.err = FALSE)
#' lines(xx, dgpd(xx, u = fit$u, sigmau = fit$sigmau, xi = fit$xi), col="red")
#' 
#' # This time with phiu
#' x = rnorm(10000)
#' xx = seq(-4, 4, 0.01)
#' hist(x, breaks = 200, freq = FALSE, xlim = c(0, 4))
#' lines(xx, dnorm(xx), lwd = 2)
#' fit = fgpd(x, u = 1, phiu = NULL, std.err = FALSE)
#' lines(xx, dgpd(xx, u = fit$u, sigmau = fit$sigmau, xi = fit$xi, phiu = fit$phiu),
#'   col = "red", lwd = 2)
#' legend("topright", c("True Density","Fitted Density"), col=c("black", "red"), lty = 1)

# maximum likelihood fitting for GPD
fgpd <- function(x, u = 0, phiu = NULL, pvector = NULL, std.err = TRUE, method = "BFGS",
  finitelik = TRUE, ...) {

  call <- match.call()
    
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")

  if (any(is.infinite(x)))
    stop("infinite cases must be removed")
  
  if (!is.logical(finitelik))
    stop("finitelik must be logical")

  if ((method == "L-BFGS-B") | (method == "BFGS"))
    finitelik = TRUE
    
  if (any(!is.finite(u) | is.logical(u)))
    stop("parameters must be numeric")
  
  if (any(is.na(x)))
    warning("missing values are treated as below threshold when estimating tail fraction")

  # assume NA or NaN are below threshold consistent with evd library
  # hence use which() to ignore these

  # check for x values in range
  whichexc = which(x > u)
  nu = length(whichexc)
  n = length(x)
  
  if (nu < 1)
    stop("no elements of x are above threshold")
    
  if (is.null(pvector)) {
    yu = x[whichexc] - u
    pvector[1] = sqrt(6 * var(yu)) / pi
    pvector[2] = 0.1
  } else {
    if (length(pvector)!=2)
      stop("Initial values for two GPD parameters must be specified")
    if (any(!is.finite(pvector)) | is.logical(pvector))
      stop("initial parameters must be numeric")
  }
  
  if (is.null(phiu)) {
    # assume NA or NaN are below threshold consistent with evd library
    # hence use which() to ignore these when estimating phiu
    phiu = nu / n
    sephiu = sqrt(phiu * (1 - phiu) / n)
  } else {
    if ((phiu < 0)  | (phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
    sephiu = NA
  }
  
  nllh = nlgpd(pvector, x = x, u = u, phiu = phiu, finitelik = finitelik)
  if (is.infinite(nllh))
    stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nlgpd, x = x, u = u, phiu = phiu, finitelik = finitelik,
    method = method, hessian = TRUE, ...)
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
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
    u = u, sigmau = fit$par[1], xi = fit$par[2], phiu = phiu)
}



