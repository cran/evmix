#' @name kden
#' 
#' @title Kernel Density Estimation Using Normal Kernel
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the kernel density estimation using the normal kernel
#'   with a constant bandwidth \code{lambda}. The kernel centres (typically the data) 
#'   are given by \code{kerncentres}.
#'
#' @inheritParams gpd
#' @param lambda       bandwidth for normal kernel (standard deviation of normal)
#' @param kerncentres  kernel centres (typically sample data)
#' 
#' @details Kernel density estimation using normal density as kernel. 
#' 
#' The density function \code{\link[evmix:kden]{dkden}} produces exactly the
#' same density estimate as \code{\link[stats:density]{density}} when a sequence
#' of \code{x} values are provided, see examples. The latter function is far
#' more efficient in this situation as it takes advantage of the computational
#' savings from doing the kernel smoothing in the spectral domain, where the 
#' convolution becomes a multiplication. So even after accounting for applying
#' the (Fast) Fourier Transform (FFT) and its inverse it is much more efficient
#' especially for a large sample size.
#' 
#' However, this KDE function applies the less efficient convolution using the
#' classic definition:
#' \deqn{\hat{f}_(x) = \frac{1}{n} \sum_{j=1}^{n} K(\frac{x - x_j}{\lambda})}
#' where \eqn{K(.)} is the density function for the standard
#' normal. Thus is no restriction on the values \code{x} can take. Computationally
#' for a particular \code{x} the density is then just
#' \code{mean(dnorm(x, kerncentres, lambda))} for the density and
#' \code{mean(pnorm(x, kerncentres, lambda))} for cumulative distribution
#' function. The random number generation is achieved by treating the KDE as a
#' mixture model with equal probability of coming from each kernel, given by
#' \code{rnorm(rep(1, n), sample(kerncentres, n, replace = TRUE), sd = lambda)}. 
#' The \code{sample()} function decides which kernel each of the \code{n} generated
#' samples comes from and gives the normal random number generator the kernel center
#' as its mean, with \code{lambda} as the kernel standard deviation.
#' 
#' The quantile function is rather more complicated as there is no closed form solution,
#' so is typically obtained by approximation or numerical solution to \eqn{P(X \le x_p) = p}
#' to find \eqn{x_p}. The quantile function \code{\link[evmix:kden]{qkden}} evaluates the
#' KDE cumulative distribution function over the range from 
#' \code{c(min(kerncentre) - 5*lambda, max(kerncentre) - 5*lambda)} as for normal kernel the
#' probability of being outside this range is of the order \code{1e-7}. Outside of the range the
#' quantiles are set to \code{-Inf} for lower tail and \code{Inf} for upper tail. A sequence of values
#' of length fifty times the number of kernels is first calculated. Spline based
#' interpolation using \code{\link[stats:splinefun]{splinefun}}, with default \code{fmm} method,
#' is then used to approximate the quantile function. This is a similar approach to that taken
#' by Matt Wand in the \code{\link[ks:kde.1d]{qkde}} in the \code{\link[ks:kde.1d]{ks}} package.
#' 
#' If no bandwidth is provided \code{lambda=NULL} then the normal reference rule is used,
#' from the function \code{\link[stats:bandwidth]{bw.nrd0}}, which is consistent with the
#' \code{\link[stats:density]{density}} function. At least two kernel centres must be provided
#' as the variance needs to be estimated.
#' 
#' @return \code{\link[evmix:kden]{dkden}} gives the density, 
#' \code{\link[evmix:kden]{pkden}} gives the cumulative distribution function,
#' \code{\link[evmix:kden]{qkden}} gives the quantile function and 
#' \code{\link[evmix:kden]{rkden}} gives a random sample.
#' 
#' @note Unlike all the other extreme value mixture model functions the 
#'   \code{\link[evmix:kden]{kden}} functions have not been vectorised as
#'   this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#'   must be either a scalar or a vector, which also define the output length.
#' 
#' The kernel centres \code{kerncentres} can either be a single datapoint or a vector
#' of data. The kernel centres (\code{kerncentres}) and locations to evaluate density (\code{x})
#' and cumulative distribution function (\code{q}) would usually be different.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{kerncentres}, \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:kden]{rkden}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
#' 
#' Due to symmetry, the lower tail can be described by GPD by negating the quantiles. 
#' The normal mean \code{nmean} and GPD threshold \code{u} will also require negation.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
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
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}. Based on code
#' by Anna MacDonald produced for MATLAB.
#'
#' @seealso \code{\link[stats:density]{density}}, \code{\link[stats:bandwidth]{bw.nrd0}}
#' and \code{\link[ks:kde.1d]{dkde}} in \code{\link[ks:kde.1d]{ks}} package.
#' 
#' @aliases kden dkden pkden qkden rkden
#' @family  kden
#' 
#' @examples
#' \dontrun{
#' nk=50
#' x = rnorm(nk)
#' xx = seq(-5, 5, 0.01)
#' plot(xx, dnorm(xx))
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = bw.nrd0(x))*0.05)
#' lines(xx, dkden(xx, x), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "KDE Using evmix", "KDE Using density function"),
#' lty = c(1, 1, 2), lwd = c(1, 2, 2), col = c("black", "red", "green"))
#' 
#' # Estimate bandwidth using cross-validation likelihood
#' fit = fkden(x)
#' hist(x, nk/5, freq = FALSE, xlim = c(-5, 5)) 
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
#' plot(xx, pnorm(xx), type = "l")
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit$lambda)*0.05)
#' lines(xx, pkden(xx, x), lwd = 2, col = "red")
#' lines(xx, pkden(xx, x, lambda = fit$lambda), lwd = 2, col = "green")
#' # green and blue (quantile) function should be same
#' p = seq(0, 1, 0.001)
#' lines(qkden(p, x, lambda = fit$lambda), p, lwd = 2, lty = 2, col = "blue") 
#' legend("topleft", c("True Density", "KDE using evmix, normal reference rule",
#' "KDE using evmix, c-v likelihood","KDE quantile function, c-v likelihood"),
#' lty = c(1, 1, 1, 2), lwd = c(1, 2, 2, 2), col = c("black", "red", "green", "blue"))
#' 
#' xnew = rkden(1000, x, lambda = fit$lambda)
#' hist(xnew, breaks = 200, freq = FALSE, xlim = c(-5, 5))
#' rug(xnew)
#' lines(xx,dnorm(xx), col = "black")
#' lines(xx, dkden(xx, x), lwd = 2, col = "red")
#' legend("topright", c("True Density", "KDE Using evmix"),
#' lty = c(1, 2), lwd = c(1, 2), col = c("black", "red"))
#'  }
NULL

#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# density function for kernel density estimator using normal kernel
dkden <- function(x, kerncentres, lambda = NULL, log = FALSE) {

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

  if (is.null(lambda)) {
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10){
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    lambda = bw.nrd0(kerncentres)
  }

  if (length(lambda) != 1)
    stop("bandwidth must be scalar")

  if ((mode(lambda) != "numeric") & !is.finite(lambda))
    stop("bandwidth must be numeric")
  
  if (lambda <= 0)
      stop("bandwidth must be non-negative")  

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  d = sapply(x, FUN = kdenx, kerncentres = kerncentres, lambda = lambda) 

  if (log) d = log(d)

  d
}
