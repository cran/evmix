#' @name gng
#' 
#' @title Normal Bulk with GPD Upper and Lower Tails Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with normal
#'   for bulk distribution between the upper and lower thresholds with
#'   conditional GPD's for the two tails. The parameters are the normal mean
#'   \code{nmean} and standard deviation \code{nsd}, lower tail (threshold \code{ul}, 
#'   GPD scale \code{sigmaul} and shape \code{xil} and tail fraction \code{phiul})
#'   and upper tail (threshold \code{ur}, GPD scale \code{sigmaur} and shape 
#'   \code{xiR} and tail fraction \code{phiuR}).
#'
#' @inheritParams normgpd
#' @inheritParams gpd
#' @param ul         lower tail threshold
#' @param sigmaul    lower tail GPD scale parameter (non-negative)
#' @param xil        lower tail GPD shape parameter
#' @param phiul      probability of being below lower threshold (0,1)
#' @param ur         upper tail threshold
#' @param sigmaur    upper tail GPD scale parameter (non-negative)
#' @param xir        upper tail GPD shape parameter
#' @param phiur      probability of being above upper threshold (0,1)
#' 
#' @details Extreme value mixture model combining normal distribution for the bulk
#' between the lower and upper thresholds and GPD for upper and lower tails. The
#' user can pre-specify \code{phiul} and \code{phiur} permitting a parameterised
#' value for the lower and upper tail fraction respectively. Alternatively, when
#' \code{phiul=TRUE} or \code{phiur=TRUE} the corresponding tail fraction is
#' estimated as from the normal bulk model.
#' 
#' Notice that the tail fraction cannot be 0 or 1, and the sum of upper and lower tail
#' fractions \code{phiul+phiur<1}, so the lower threshold must be less than the upper, 
#' \code{ul<ur}.
#' 
#' The cumulative distribution function now has three components. The lower tail with 
#' tail fraction \eqn{\phi_{ul}} defined by the normal bulk model (\code{phiul=TRUE})
#' upto the lower threshold \eqn{x < u_l}:
#'   \deqn{F(x) = H(u_l) G_l(x).}
#' where \eqn{H(x)} is the normal cumulative distribution function (i.e. 
#' \code{pnorm(ur, nmean, nsd)}). The 
#' \eqn{G_l(X)} is the conditional GPD cumulative distribution function with negated
#' data and threshold, i.e. \code{dgpd(-x, -ul, sigmaul, xil, phiul)}. The normal
#' bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = H(x).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = H(u_r) + [1 - H(u_r)] G(x)}
#' where \eqn{G(X)}.
#' 
#' The cumulative distribution function for the pre-specified tail fractions 
#' \eqn{\phi_{ul}} and \eqn{\phi_{ur}} is more complicated.  The unconditional GPD
#' is used for the lower tail \eqn{x < u_l}:
#'   \deqn{F(x) = \phi_{ul} G_l(x).}
#' The normal bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = \phi_{ul}+ (1-\phi_{ul}-\phi_{ur}) (H(x) - H(u_l)) / (H(u_r) - H(u_l)).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = (1-\phi_{ur}) + \phi_{ur} G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_{ul} = H(u_l)} and
#' \eqn{\phi_{ur} = 1 - H(u_r)}.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component, 
#' \code{\link[stats:Normal]{dnorm}} for details of normal bulk component and
#' \code{\link[evmix:normgpd]{dnormgpd}} for normal with GPD extreme value
#' mixture model.
#' 
#' @return \code{\link[evmix:gng]{dgng}} gives the density, 
#' \code{\link[evmix:gng]{pgng}} gives the cumulative distribution function,
#' \code{\link[evmix:gng]{qgng}} gives the quantile function and 
#' \code{\link[evmix:gng]{rgng}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main input (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{rgng} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:gng]{rgng}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x} and \code{q}
#' are passed through as is and infinite values are set to \code{NA}.
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
#' Zhao, X., Scarrott, C.J. Reale, M. and Oxley, L. (2010). Extreme value modelling
#' for forecasting the market crisis. Applied Financial Econometrics 20(1), 63-72.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:normgpd]{normgpd}}, \code{\link[evmix:gpd]{gpd}}
#' and \code{\link[stats:Normal]{dnorm}}
#' @aliases  gng dgng pgng qgng rgng
#' @family   gng
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2))
#' x = rgng(1000, phiul = 0.15, phiur = 0.15)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgng(xx, phiul = 0.15, phiur = 0.15))
#' 
#' # three tail behaviours
#' plot(xx, pgng(xx), type = "l")
#' lines(xx, pgng(xx, xil = 0.3, xir = 0.3), col = "red")
#' lines(xx, pgng(xx, xil = -0.3, xir = -0.3), col = "blue")
#' legend("topleft", paste("Symmetric xil=xir=",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rgng(1000, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgng(xx, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2))
#' 
#' plot(xx, dgng(xx, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2), type = "l", ylim = c(0, 0.4))
#' lines(xx, dgng(xx, xil = -0.3, phiul = 0.3, xir = 0.3, phiur = 0.3), col = "red")
#' lines(xx, dgng(xx, xil = -0.3, phiul = TRUE, xir = 0.3, phiur = TRUE), col = "blue")
#' legend("topleft", c("phiul = phiur = 0.2", "phiul = phiur = 0.3", "Bulk Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   }
NULL

#' @export
#' @aliases  gng dgng pgng qgng rgng
#' @rdname gng

# probability density function for normal bulk with GPD's for upper and lower tails
dgng <- function(x, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), sigmaul = nsd, xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), sigmaur = nsd, xir = 0, phiur = TRUE, log = FALSE) {

  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
    
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    warning("infinite cases are set to NaN")

  x[is.infinite(x)]=NA # user will have to deal with infinite cases

  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(x), length(nmean), length(nsd),
    length(ul), length(sigmaul), length(xil), length(phiul),
    length(ur), length(sigmaur), length(xir), length(phiur))
  n = max(linputs)

  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")

  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" |
    mode(ul) != "numeric" | mode(sigmaul) != "numeric" | mode(xil) != "numeric" |
    mode(ur) != "numeric" | mode(sigmaur) != "numeric" | mode(xir) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur))))
    stop("parameters must be numeric")

  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")

  if (any(ul >= ur))
    stop("lower threshold must be below upper threshold")

  if (min(sigmaul) <= 0 | min(sigmaur) <= 0)
    stop("scale must be non-negative")
  
  if ((is.logical(phiul) & any(!phiul)) | (is.logical(phiur) & any(!phiur))) {
    stop("phiu must be either TRUE for bulk parameterised threshold probability approach, 
      or between 0 and 1 (exclusive) when using parameterised threshold probability approach")
  } else {
    if (any(phiul < 0) | any(phiul > 1) | any(phiur < 0) | any(phiur > 1))
      stop("phiu must between 0 and 1 (inclusive)")
    if (!is.logical(phiul) & !is.logical(phiur)) {
      if (any((phiul + phiur) > 1))
        stop("phiu + phiur must be less than 1")
    }
  }

  if (!is.logical(log))
    stop("log must be logical")
  
  if (length(log) != 1)
    stop("log must be of length 1")

  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  ul = rep(ul, length.out = n)
  sigmaul = rep(sigmaul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  sigmaur = rep(sigmaur, length.out = n)
  xir = rep(xir, length.out = n)
  
  if (is.logical(phiul)) {
    phiul = pnorm(ul, mean = nmean, sd = nsd)
  } else {
    phiul = rep(phiul, length.out = n)
  }
  if (is.logical(phiur)) {
    phiur = 1 - pnorm(ur, mean = nmean, sd = nsd)
  } else {
    phiur = rep(phiur, length.out = n)
  }
  phib = (1 - phiul - phiur) / (pnorm(ur, mean = nmean, sd = nsd) - pnorm(ul, mean = nmean, sd = nsd))

  d = x # this will pass through NA/NaN in x just as they are entered
  
  whichul = which(x < ul)
  nul = length(whichul)
  whichb = which((x <= ur) & (x >= ul)) 
  nb = length(whichb)
  whichur = which(x > ur)
  nur = length(whichur)
  
  if (nul > 0) d[whichul] = log(phiul[whichul]) + dgpd(-x[whichul], -ul[whichul], sigmaul[whichul], xil[whichul], log = TRUE)
  if (nb > 0) d[whichb] = log(phib[whichb]) + dnorm(x[whichb], mean = nmean[whichb], sd = nsd[whichb], log = TRUE)
  if (nur > 0) d[whichur] = log(phiur[whichur]) + dgpd(x[whichur], ur[whichur], sigmaur[whichur], xir[whichur], log = TRUE)

  if (!log) d = exp(d)

  d
}
