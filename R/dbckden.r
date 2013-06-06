#' @name bckden
#' 
#' @title Boundary Corrected Kernel Density Estimation
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the boundary corrected kernel density estimators
#'   with a constant bandwidth \code{lambda}. The kernel
#'   centres (typically the data) are given by \code{kerncentres}.
#'
#' @inheritParams kden
#' @param lambda     scalar value of fixed bandwidth, or \code{NULL} (default)
#' @param bcmethod   boundary correction approach
#' @param proper     logical, should density be renormalised to integrate to unity, 
#'                   simple boundary correction only
#' @param nn         non-negativity correction, so simple boundary correction only
#' @param xmax       upper bound on support, for copula and beta kernels only
#' @param offset     offset added to kernel centres, for logtrans
#'
#' @details Boundary corrected kernel density estimation (KDE) to improve the bias properties
#' near the boundary. A fairly wide range of methods are implemented for the user to choose from
#' to cope with a lower boundary at zero and potentially also both upper and lower boundaries. Some boundary
#' correction methods require a secondary correction for negative density estimates, so an option
#' is available. Further, some methods also need to be normalised
#' to ensure the density estimate is proper (i.e. integrates to one), so an option
#' is provided to renormalise.
#' 
#' It assumes there is a lower boundary at zero, so prior transformation of the data would
#' be required for alternative boundaries (including negation to use it for only an upper boundary).
#' 
#' Renormalisation of the kernel to integrate to unity is assumed by default
#' (\code{proper=TRUE}), but the user can specify if the raw density estimate is provided
#' instead (\code{proper=FALSE}). For the methods implemented thus far, this is needed for
#' \code{bcmethod="simple"} which can be evaluated in closed form, and \code{bcmethod="beta1"}
#' or \code{bcmethod="beta2"} which require numerical integration.
#' 
#' Correction of the density estimate to ensure non-negativity can be applied,
#' which is only relevant for the \code{bcmethod="simple"} approach. The Jones
#' and Foster (1996) method will be applied (\code{nn="jf96"}) by default. The
#' non-negative value can simply be zeroed (\code{nn="zero"}).  By default,
#' correction is applied (\code{nn="none"}). This method can occassionally give
#' an extra boundary bias for certain populations (e.g. Gamma(2, 1)), see their
#' paper for details. Renormlisation should be used after these non-negativity
#' corrections.
#' 
#' The non-negative correction is applied before renormalisation (when either is requested). 
#' 
#' The boundary correction methods implemented are:
#' 
#' \code{bcmethod="simple"} is the default and applies the simple boundary correction method
#' in equation (3.4) of Jones (1993) and is equivalent to the kernel weighted local linear
#' fitting at the boundary. Normal kernels are used.
#' 
#' \code{bcmethod="renorm"} applies the renormalisation method discussed in
#' Diggle (1985), where the kernels are simply truncated at the boundary and 
#' renormalised to unity. But this still exhibits a o(h) boundary bias. Normal kernels are used.
#' 
#' \code{bcmethod="reflect"} applies the reflection method of Boneva, Kendall and Stefanov
#' (1971) which is equivalent to the dataset being supplemented by the same dataset negated. 
#' This method implicitly assumes f'(0)=0, so can causes extra artefacts at the boundary. 
#' Normal kernels are used.
#' 
#' \code{bcmethod="logtrans"} applies KDE on the log-scale and then transforms back (with
#' explicit normalisation), following Marron and Ruppert (1992). This is the approach implmented
#' in the \code{\link[ks:kde.1d]{ks}} package. As the KDE is applied on the log scale, the effective
#' bandwidth on the original scale is not constant. Normal kernels are used on the log-scale. 
#' The \code{offset} option is only used for this method, to offset zero values, to prevent
#' \code{log(0)}.
#' 
#' \code{bcmethod="beta1"} and \code{"beta2"} due to Chen (1999) which uses beta kernels 
#' and modified beta kernels respectively in the KDE. The \code{xmax} argument
#' has been provided so that the user can have the beta kernels rescaled to be appropriate
#' for a support of [0, xmax] rather than [0, 1].
#' 
#' \code{bcmethod="gamma1"} and \code{"gamma2"} due to Chen (2000) which uses gamma kernels 
#' and modified gamma kernels respectively in the KDE.
#' 
#' \code{bcmethod="copula"} due to Jones and Henderson (2007) uses bivariate normal copula based 
#' kernels in the KDE, essentially by taking condition slices thorugh the bivariate copula at
#' each kernel centre. As with the \code{bcmethod="beta"} option the \code{xmax} argument
#' has been provided to rescale the kernels over [0, xmax] rather than [0, 1]. In this case
#' the bandwidth is defined as \eqn{lambda=1-\rho^2}, so is limited to (0, 1).
#' 
#' The examples below show a trick you can use to see the actual kernels used in the chosen
#' boundary correction method.
#' 
#' The quantile function is rather complicated as there is typically no closed form solution,
#' so in these cases it is obtained by approximation or numerical solution to \eqn{P(X \le x_p) = p}
#' to find \eqn{x_p}. The quantile function \code{\link[evmix:bckden]{qbckden}} evaluates the
#' KDE cumulative distribution function over the range from 
#' \code{c(0, max(kerncentre) - 5*lambda)}. Outside of this range the
#' quantiles are set to \code{0} for lower tail and \code{Inf} for upper tail. A sequence of values
#' of length fifty times the number of kernels is first calculated. Spline based
#' interpolation using \code{\link[stats:splinefun]{splinefun}}, with default \code{monoH.FC} method,
#' is then used to approximate the quantile function. This is a similar approach to that taken
#' by Matt Wand in the \code{\link[ks:kde.1d]{qkde}} in the \code{\link[ks:kde.1d]{ks}} package.
#' 
#' Unlike the standard KDE, there is no general rule-of-thumb bandwidth for all these
#' estimators, with only certain methods having a guideline in the literature, so none have been
#' implmented. Hence, the user has to specify a bandwidth, but should consider using 
#' \code{\link[evmix:fbckden]{fbckden}} function for cross-validation likelihood fitting.
#' 
#' Random number generation is slow as inversion sampling using the (numerically evaluated)
#' quantile function is implemented. Users may want to consider alternative approaches instead,
#' like rejection sampling. 
#' 
#' @return \code{\link[evmix:bckden]{dbckden}} gives the density, 
#' \code{\link[evmix:bckden]{pbckden}} gives the cumulative distribution function,
#' \code{\link[evmix:bckden]{qbckden}} gives the quantile function and 
#' \code{\link[evmix:bckden]{rbckden}} gives a random sample.
#' 
#' @note Unlike all the other extreme value mixture model functions the \code{bckden}
#' functions have not been vectorised as this is not appropriate. The main inputs
#' (\code{x}, \code{p} or \code{q}) must be either a scalar or a vector, which also
#' define the output length.
#' 
#' The kernel centres \code{kerncentres} can either be a single datapoint or a vector
#' of data. The kernel centres (\code{kerncentres}) and locations to evaluate density (\code{x})
#' and cumulative distribution function (\code{q}) would usually be different.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{lambda}, \code{kerncentres}, \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:bckden]{rbckden}} is 1.
#' 
#' The \code{xmax} option is only relevant for the beta and copula methods, so a warning is produced
#' if this is changed from the default in other methods.
#' 
#' The renormalisation is only relevant 
#' for the \code{bcmethod="simple"}, \code{"beta1"} and \code{"beta2"} approaches, so will not be
#' applied in any other case (even if the user specifies \code{proper=TRUE}).
#' 
#' Non-negative correction will be applied by default, but this is only relevant for
#' the \code{bcmethod="simple"} approach. It will not be applied in any other cases (even
#' if the user specifies a method for \code{nn}).
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
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
#' Chen, S.X. (1999). Beta kernel estimators for density functions. Computational Statistics
#' and Data Analysis 31, 1310-45.
#' 
#' Chen, S.X. (2000). Probability density function estimation using gamma kernels.
#' Annals of the Institute of Statisical Mathematics 52(3), 471-480.
#' 
#' Boneva, L.I., Kendall, D.G. and Stefanov, I. (1971). Spline transformations: Three new
#' diagnostic aids for the statistical data analyst (with discussion). Journal of the Royal
#' Statistical Society B, 33, 1-70.
#' 
#' Diggle, P.J. (1985). A kernel method for smoothing point process data. Applied Statistics
#' 34, 138-147.
#' 
#' Marron, J.S. and Ruppert, D. (1994) Transformations to reduce boundary bias in kernel
#' density estimation, Journal of the Royal Statistical Society. Series B 56(4), 653-671.
#' 
#' Jones, M.C. and Henderson, D.A. (2007). Kernel-type density estimation on the unit
#' interval. Biometrika 94(4), 977-984.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}. Based on code
#' by Anna MacDonald produced for MATLAB.
#'
#' @seealso \code{\link[evmix:kden]{kden}}, \code{\link[stats:density]{density}}, 
#' \code{\link[logspline:logspline]{logspline}}
#' and \code{\link[ks:kde.1d]{dkde}}.
#' 
#' @aliases bckden dbckden pbckden qbckden rbckden
#' @family  bckdengpd bckden kden
#' 
#' @examples
#' \dontrun{
#' n=100
#' x = rgamma(n, shape = 1, scale = 2)
#' xx = seq(-0.5, 12, 0.01)
#' plot(xx, dgamma(xx, shape = 1, scale = 2), type = "l")
#' rug(x)
#' lines(xx, dbckden(xx, x, lambda = 0.3), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' 
#' #  # Trick to show actual kernels on the plot - just add one kernel centre at a time:
#' for (i in 1:min(n,20)){
#'   lines(xx, dbckden(xx, x[i], lambda = 0.3, proper = FALSE)*0.05, col = "blue")
#' }
#' # Notice the negative weights in the kernels for this approach
#' 
#' legend("topright", c("True Density", "Simple boundary correction", "KDE using density function",
#' "Boundary Corrected Kernels"),
#' lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))
#'
#' n=100
#' x = rbeta(n, shape1 = 3, shape2 = 2)*5
#' xx = seq(-0.5, 5.5, 0.01)
#' plot(xx, dbeta(xx/5, shape1 = 3, shape2 = 2)/5, type = "l")
#' rug(x)
#' lines(xx, dbckden(xx, x, lambda = 0.01, bcmethod = "beta2", xmax = 5),
#'   lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "Modified Beta KDE Using evmix",
#'   "KDE using density function"),
#' lty = c(1, 1, 2), lwd = c(1, 2, 2), col = c("black", "red", "green"))
#'
#' # Demonstrate renormalisation (usually small difference)
#' n=100
#' x = rgamma(n, shape = 2, scale = 2)
#' xx = seq(-0.5, 15, 0.01)
#' plot(xx, dgamma(xx, shape = 2, scale = 2), type = "l")
#' rug(x)
#' lines(xx, dbckden(xx, x, lambda = 0.5, bcmethod = "simple", proper = TRUE),
#'   lwd = 2, col = "red")
#' lines(xx, dbckden(xx, x, lambda = 0.5, bcmethod = "simple", proper = FALSE),
#'   lwd = 2, col = "purple")
#' legend("topright", c("True Density", "Simple BC with renomalisation", 
#' "Simple BC without renomalisation"),
#' lty = 1, lwd = c(1, 2, 2), col = c("black", "red", "purple"))
#' }
NULL

#' @export
#' @aliases bckden dbckden pbckden qbckden
#' @rdname  bckden

# density function for boundary corrected KDE
dbckden <- function(x, kerncentres, lambda = NULL, bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = 0, xmax = Inf, log = FALSE) {

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

  if (any(kerncentres < 0))
    stop("kernel centres cannot be non-positive")
  
  if (is.null(lambda))
    stop("bandwidth (lambda) must be specified")

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
  
  allmethods = c("simple", "renorm", "reflect", "logtrans", 
    "beta1", "beta2", "gamma1", "gamma2", "copula")
  if (length(bcmethod) != 1)
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if ((mode(bcmethod) != "character"))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))
  
  if (!(bcmethod %in% allmethods))
    stop(paste("bcmethod must be one of", allmethods, collapse = " "))

  if ((bcmethod == "copula") & (lambda >= 1))
      stop("bandwidth must between (0, 1) for copula method")  

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
  
  # problems caused in evaluation if no data near boundary
  # so don't evaluate cdf if no kernel centres near threshold
  if (bcmethod %in% upboundmethods) {
    minaccept = 0
    maxaccept = xmax
  } else if (bcmethod == "logtrans") {
    maxaccept = exp(log(max(kerncentres) + offset) + 5*lambda)
    minaccept = max(0, exp(log(min(kerncentres) + offset) - 5*lambda))
  } else {
    maxaccept = max(kerncentres) + 5*lambda
    minaccept = max(0, min(kerncentres) - 5*lambda)
  }

  d = x # this will pass through NA/NaN in x just as they are entered
  
  d[!is.na(d)] = 0 # default for outside of range of support

  # only select non-negative x-values for evaluation
  xok = x[!is.na(x)]
  xok = xok[(xok >= minaccept) & (xok <= maxaccept)]

  if (length(xok) > 0) {
    if (bcmethod == "simple") {
      # simple linear boundary correction method of Jones (1993), eq 3.4
      # which adapts normal kernels so they are sensibly behaved near boundary
      # use similar notation to original paper
      # (truncpoint is their p and lambda is their S_k)
      
      dok = sapply(xok, FUN = simplebckdenx, kerncentres = kerncentres, lambda = lambda)
      
    } else if (bcmethod == "renorm") {
      # really crude renormalisation, based on how much of kernel is above boundary
      
      dok = sapply(xok, FUN = renormbckdenx, kerncentres = kerncentres, lambda = lambda)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "reflect") {
      # really crude reflection, equivalent to using (kerncentres, -kerncentres)
      # only good if f'(0)=0
      
      dok = sapply(xok, FUN = reflectbckdenx, kerncentres = kerncentres, lambda = lambda)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "logtrans") {
      # simple log transformation 
      # (Wand, Marron and Ruppert JASA, 1991, Marron and Ruppert JRSS B, 1994)
      # requires offset in min(x)=0
      if ((min(kerncentres) == 0) & (offset = 0)){
        stop("log transformation requires offset for zero kernel centre")
      }
      
      # transformation KDE is f(x) = g(t(x)) abs(t'(x))
      dok = sapply(log(xok + offset), FUN = kdenx, 
        kerncentres = log(kerncentres), lambda = lambda) / (xok+offset)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta1") {
      # beta kernels by Chen (1999) CSDA
      
      dok = sapply(xok, FUN = beta1bckdenx, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "gamma1") {
      # gamma kernels by Chen (2000) AISM
      
      dok = sapply(xok, FUN = gamma1bckdenx, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    } else if (bcmethod == "gamma2") {
      # modified gamma kernels by Chen (2000) AISM
      
      dok = sapply(xok, FUN = gamma2bckdenx, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta2") {
      # modified beta kernels by Chen (1999) CSDA
      
      dok = sapply(xok, FUN = beta2bckdenx, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "copula") {
      # bivariate normal copula based kernels by Jones and Henderson (2007) Biometrika
      
      dok = sapply(xok, FUN = copulabckdenx, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    }
    
    # Apply non-negative correction first followed by renormalisation (if either requested)
    if (nn == "jf96") {
      dbar = sapply(xok, FUN = nnbckdenx, kerncentres = kerncentres, lambda = lambda)
      dok = ifelse(dbar == 0, dok, dbar*exp(dok/dbar - 1))
    } else if (nn =="zero") {
      dok[which(dok < 0)] = 0
    }

    if (proper) {
      if (bcmethod == "simple") {
        pxmax = simplepbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda)
      } else if (bcmethod == "beta1") {
        pxmax = beta1pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "beta2") {
        pxmax = beta2pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "gamma1") {
        pxmax = gamma1pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda)
      } else if (bcmethod == "gamma2") {
        pxmax = gamma2pbckdenx(maxaccept, kerncentres = kerncentres, lambda = lambda)
      }
      dok = dok/pxmax
    }
    d[ifelse(!is.na(x), (x >= minaccept) & (x <= maxaccept), FALSE)] = dok
  }
  
  if (log) d = log(d)

  d
}
