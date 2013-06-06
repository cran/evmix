#' @name kdengpd
#' 
#' @title Kernel Density Estimation Using Normal Kernel and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the kernel density estimation using normal kernel for the bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the bandwidth \code{lambda}, threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams kden
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining  kernel density estimation Using normal kernel for the bulk
#' below the threshold and GPD for upper tail. The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' normal bulk model.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the kernel density estimation using normal kernel (\code{phiu=TRUE}), upto the 
#' threshold \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the kernel density estimator and conditional GPD
#' cumulative distribution functions (i.e. \code{mean(pnorm(x, kerncentres, lambda))} and
#' \code{pgpd(x, u, sigmau, xi)}).
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[stats:Normal]{dnorm}} for details of normal bulk component.
#' 
#' @return \code{\link[evmix:kdengpd]{dkdengpd}} gives the density, 
#' \code{\link[evmix:kdengpd]{pkdengpd}} gives the cumulative distribution function,
#' \code{\link[evmix:kdengpd]{qkdengpd}} gives the quantile function and 
#' \code{\link[evmix:kdengpd]{rkdengpd}} gives a random sample.
#' 
#' @note Unlike all the other extreme value mixture model functions the 
#'   \code{\link[evmix:kdengpd]{kdengpd}} functions have not been vectorised as
#'   this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#'   must be either a scalar or a vector, which also define the output length.
#'   The \code{kerncentres} can also be a scalar or vector.
#' 
#' The kernel centres \code{kerncentres} can either be a single datapoint or a vector
#' of data. The kernel centres (\code{kerncentres}) and locations to evaluate density (\code{x})
#' and cumulative distribution function (\code{q}) would usually be different.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{kerncentres}, \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:kdengpd]{rkdengpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
#' 
#' Due to symmetry, the lower tail can be described by GPD by negating the quantiles. 
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
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
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Normal]{dnorm}}
#' @aliases  kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @family   kdengpd
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2))
#' kerncentres=rnorm(500, 0, 1)
#' xx = seq(-4, 4, 0.01)
#' hist(kerncentres, breaks = 100, freq = FALSE)
#' lines(xx, dkdengpd(xx, kerncentres, u = 1.2, sigmau = 0.56, xi = 0.1))
#' 
#' plot(xx, pkdengpd(xx, kerncentres), type = "l")
#' lines(xx, pkdengpd(xx, kerncentres, xi = 0.3), col = "red")
#' lines(xx, pkdengpd(xx, kerncentres, xi = -0.3), col = "blue")
#' legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
#'       col=c("black", "red", "blue"), lty = 1, cex = 0.5)
#'
#' kerncentres=rnorm(1000, 0, 1)
#' x = rkdengpd(1000, kerncentres, phiu = 0.1, u = 1.2, sigmau = 0.56, xi = 0.1)
#' xx = seq(-4, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
#' lines(xx, dkdengpd(xx, kerncentres, phiu = 0.1))
#'
#' plot(xx, dkdengpd(xx, kerncentres, xi=0, phiu = 0.2), type = "l")
#' lines(xx, dkdengpd(xx, kerncentres, xi=-0.2, phiu = 0.2), col = "red")
#' lines(xx, dkdengpd(xx, kerncentres, xi=0.2, phiu = 0.2), col = "blue")
#' legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'       col=c("black", "red", "blue"), lty = 1)
#'}
NULL

#' @export
#' @aliases kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @rdname kdengpd

# probability density function for kernel density estimation using normal kernel for the bulk
# distribution upto the threshold and conditional GPD above threshold.
dkdengpd <- function(x, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)), 
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, log = FALSE) {
  
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    warning("infinite cases are set to NaN")
  
  x[is.infinite(x)] = NaN # user will have to deal with infinite cases
  
  if (missing(kerncentres))
    stop("kerncentres must be a non-empty numeric vector")
  
  if (length(kerncentres) == 0 | mode(kerncentres) != "numeric") 
    stop("kerncentres must be a non-empty numeric vector")
  
  if (any(!is.finite(kerncentres)))
    warning("non-finite kernel centres are dropped")
  
  kerncentres = kerncentres[is.finite(kerncentres)]
  nk = length(kerncentres)
  
  if (is.null(lambda)){
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10){
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    lambda = bw.nrd0(kerncentres)
  }
  
  linputs = c(length(lambda), length(u), length(sigmau), length(xi), length(phiu))

  if (sum(linputs != 1) > 0)
    stop("parameters must be scalar")
  
  if (mode(lambda) != "numeric" | mode(u) != "numeric" | mode(sigmau) != "numeric" |
      mode(xi) != "numeric")
    stop("parameters must be numeric")

  if (any(!is.finite(c(lambda, u, sigmau, xi, phiu))))
    stop("parameters must be numeric")  
  
  if (lambda <= 0)
    stop("bandwidth must be non-negative")  
  
  if (sigmau <= 0)
    stop("scale must be non-negative")
  
  if (is.logical(phiu) & (!phiu)) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
         or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if ((phiu < 0)| (phiu > 1))
      stop("phiu must between 0 and 1 (inclusive)")
  }
  
  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")
   
  if (is.logical(phiu)) {
    phiu = 1 - pkdenx(u, kerncentres, lambda)
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pkdenx(u, kerncentres, lambda)

  d = x # this will pass through NA/NaN in x just as they are entered

  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)

  if (nb > 0) d[whichb] = log(phib) + log(sapply(x[whichb], FUN = kdenx, kerncentres = kerncentres, lambda = lambda))
  if (nu > 0) d[whichu] = log(phiu) + dgpd(x[whichu], u, sigmau, xi, log = TRUE)
  
  if (!log) d = exp(d)
  
  d
}
