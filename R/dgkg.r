#' @name gkg
#' 
#' @title Kernel Density Estimation for Bulk and GPD for Both Upper and Lower Tails in
#' Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with kernel density
#'   estimation using normal kernel 
#'   for bulk distribution between the upper and lower thresholds with
#'   conditional GPD's for the two tails. The parameters are the kernel bandwidth
#'  \code{lambda}, lower tail (threshold \code{ul}, 
#'   GPD scale \code{sigmaul} and shape \code{xil} and tail fraction \code{phiul})
#'   and upper tail (threshold \code{ur}, GPD scale \code{sigmaur} and shape 
#'   \code{xiR} and tail fraction \code{phiur}).
#'
#' @inheritParams kdengpd
#' @inheritParams gpd
#' @inheritParams gng
#' 
#' @details Extreme value mixture model combining kernel density estimator (KDE) with normal 
#' kernels to represent the bulk between the lower and upper thresholds and GPD for
#' the upper and lower tails. The
#' user can pre-specify \code{phiul} and \code{phiur} permitting a parameterised
#' value for the lower and upper tail fraction respectively. Alternatively, when
#' \code{phiul=TRUE} or \code{phiur=TRUE} the corresponding tail fraction is
#' estimated as from the normal bulk model.
#' 
#' Notice that the tail fraction cannot be 0 or 1, and the sum of upper and lower tail
#' fractions \code{phiul + phiur < 1}, so the lower threshold must be less than the upper, 
#' \code{ul < ur}.
#' 
#' The cumulative distribution function has three components. The lower tail with 
#' tail fraction \eqn{\phi_{ul}} defined by the KDE bulk model (\code{phiul=TRUE})
#' upto the lower threshold \eqn{x < u_l}:
#'   \deqn{F(x) = H(u_l) [1 - G_l(x)].}
#' where \eqn{H(x)} is the kernel density estimator cumulative distribution function (i.e. 
#' \code{mean(pnorm(x, kerncentres, lambda))} and  
#' \eqn{G_l(X)} is the conditional GPD cumulative distribution function with negated
#' \eqn{x} value and threshold, i.e. \code{pgpd(-x, -ul, sigmaul, xil, phiul)}. The KDE
#' bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = H(x).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = H(u_r) + [1 - H(u_r)] G_r(x)}
#' where \eqn{G_r(X)} is the GPD cumulative distribution function, 
#' i.e. \code{pgpd(x, ur, sigmaur, xir, phiur)}.
#' 
#' The cumulative distribution function for the pre-specified tail fractions 
#' \eqn{\phi_{ul}} and \eqn{\phi_{ur}} is more complicated.  The unconditional GPD
#' is used for the lower tail \eqn{x < u_l}:
#'   \deqn{F(x) = \phi_{ul} [1 - G_l(x)].}
#' The KDE bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = \phi_{ul}+ (1-\phi_{ul}-\phi_{ur}) (H(x) - H(u_l)) / (H(u_r) - H(u_l)).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = (1-\phi_{ur}) + \phi_{ur} G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_{ul} = H(u_l)} and
#' \eqn{\phi_{ur} = 1 - H(u_r)}.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component, 
#' \code{\link[evmix:kden]{dkden}} for details of KDE bulk component and
#' \code{\link[evmix:kdengpd]{dkdengpd}} for KDE with sinlge upper tail GPD extreme value
#' mixture model.
#' 
#' @return \code{\link[evmix:gkg]{dgkg}} gives the density, 
#' \code{\link[evmix:gkg]{pgkg}} gives the cumulative distribution function,
#' \code{\link[evmix:gkg]{qgkg}} gives the quantile function and 
#' \code{\link[evmix:gkg]{rgkg}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main input (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{rgkg} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:gkg]{rgkg}} is 1.
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
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @aliases  kdgkg dgkg pgkg qgkg rgkg
#' @family   gkg
#' 
#' @examples
#' \dontrun{
#' par(mfrow=c(2,2))
#' kerncentres=rnorm(1000,0,1)
#' x = rgkg(1000, kerncentres, phiul = 0.15, phiur = 0.15)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgkg(xx, kerncentres, phiul = 0.15, phiur = 0.15))
#' 
#' # three tail behaviours
#' plot(xx, pgkg(xx, kerncentres), type = "l")
#' lines(xx, pgkg(xx, kerncentres,xil = 0.3, xir = 0.3), col = "red")
#' lines(xx, pgkg(xx, kerncentres,xil = -0.3, xir = -0.3), col = "blue")
#' legend("topleft", paste("Symmetric xil=xir=",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rgkg(1000, kerncentres, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2))
#' 
#' plot(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2),
#'   type = "l", ylim = c(0, 0.4))
#' lines(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = 0.3, xir = 0.3, phiur = 0.3),
#'   col = "red")
#' lines(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = TRUE, xir = 0.3, phiur = TRUE),
#'   col = "blue")
#' legend("topleft", c("phiul = phiur = 0.2", "phiul = phiur = 0.3", "Bulk Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
NULL

#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname  gkg

# probability density function for KDE for bulk with GPD's for both upper and lower tails
dgkg<- function(x, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE, 
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE,  log = FALSE) {
  
  # Check properties of inputs
  if (missing(x))
    stop("x must be a non-empty numeric vector")
  
  if (length(x) == 0 | mode(x) != "numeric") 
    stop("x must be a non-empty numeric vector")
  
  if (any(is.infinite(x)))
    warning("infinite cases are set to NaN")
  
  x[is.infinite(x)]=NA # user will have to deal with infinite cases
  
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
  
  linputs = c(length(lambda), length(ul), length(sigmaul), length(xil), length(phiul),
                              length(ur), length(sigmaur), length(xir), length(phiur)) 
  
  if (sum(linputs != 1) > 0)
    stop("parameters must be scalar")
  
  if (mode(lambda) != "numeric"| 
      mode(ul) != "numeric" | mode(sigmaul) != "numeric" | mode(xil) != "numeric" |
      mode(ur) != "numeric" | mode(sigmaur) != "numeric" | mode(xir) != "numeric" |
      !(mode(phiul) %in% c("logical","numeric")) | !(mode(phiur) %in% c("logical","numeric")))
    stop("parameters must be numeric, phiu can be numeric or logical")

  if (any(!is.finite(c(lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur))))
    stop("parameters must be numeric")  
  
  if (lambda <= 0)
    stop("bandwidth must be non-negative")  
  
  if ((sigmaul <= 0) | (sigmaur <= 0))
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

  if (is.logical(phiul)) {
    phiul = pkdenx(ul, kerncentres, lambda)
  } else {
    phiul = phiul
  }
  if (is.logical(phiur)) {
    phiur = 1 - pkdenx(ur, kerncentres, lambda)
  } else {
    phiur = phiur
  }
  phib = (1 - phiul - phiur) / (pkdenx(ur, kerncentres, lambda) - pkdenx(ul, kerncentres, lambda))
    
  d = x # this will pass through NA/NaN in x just as they are entered
  
  whichul = which(x < ul)
  nul = length(whichul)
  whichb = which((x <= ur) & (x >= ul)) 
  nb = length(whichb)
  whichur = which(x > ur)
  nur = length(whichur)
    
  if (nul > 0) d[whichul] = log(phiul) + dgpd(-x[whichul], -ul, sigmaul, xil, log = TRUE)
  if (nb > 0) d[whichb] =  log(phib) + log(dkden(x[whichb], kerncentres, lambda))
  if (nur > 0) d[whichur] = log(phiur) + dgpd(x[whichur], ur, sigmaur, xir, log = TRUE)
  
  if (!log) d = exp(d)
  
  d
}

