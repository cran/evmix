\name{kdengpd}
\alias{dkdengpd}
\alias{kdengpd}
\alias{pkdengpd}
\alias{qkdengpd}
\alias{rkdengpd}
\title{Kernel Density Estimate and GPD Tail Extreme Value Mixture Model}
\usage{
  dkdengpd(x, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, bw = NULL, kernel = "gaussian",
    log = FALSE)

  pkdengpd(q, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, bw = NULL, kernel = "gaussian",
    lower.tail = TRUE)

  qkdengpd(p, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, bw = NULL, kernel = "gaussian",
    lower.tail = TRUE)

  rkdengpd(n = 1, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, bw = NULL, kernel = "gaussian")
}
\arguments{
  \item{phiu}{probability of being above threshold \eqn{[0,
  1]} or \code{TRUE}}

  \item{x}{quantiles}

  \item{kerncentres}{kernel centres (typically sample data
  vector or scalar)}

  \item{lambda}{bandwidth for kernel (as half-width of
  kernel) or \code{NULL}}

  \item{bw}{bandwidth for kernel (as standard deviations of
  kernel) or \code{NULL}}

  \item{kernel}{kernel name (\code{default = "gaussian"})}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantiles}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probabilities}

  \item{n}{sample size (positive integer)}

  \item{u}{threshold}

  \item{sigmau}{scale parameter (positive)}

  \item{xi}{shape parameter}
}
\value{
  \code{\link[evmix:kdengpd]{dkdengpd}} gives the density,
  \code{\link[evmix:kdengpd]{pkdengpd}} gives the
  cumulative distribution function,
  \code{\link[evmix:kdengpd]{qkdengpd}} gives the quantile
  function and \code{\link[evmix:kdengpd]{rkdengpd}} gives
  a random sample.
}
\description{
  Density, cumulative distribution function, quantile
  function and random number generation for the extreme
  value mixture model with kernel density estimate for bulk
  distribution upto the threshold and conditional GPD above
  threshold. The parameters are the bandwidth
  \code{lambda}, threshold \code{u} GPD scale \code{sigmau}
  and shape \code{xi} and tail fraction \code{phiu}.
}
\details{
  Extreme value mixture model combining kernel density
  estimate (KDE) for the bulk below the threshold and GPD
  for upper tail.

  The user can pre-specify \code{phiu} permitting a
  parameterised value for the tail fraction \eqn{\phi_u}.
  Alternatively, when \code{phiu=TRUE} the tail fraction is
  estimated as the tail fraction from the KDE bulk model.

  The alternate bandwidth definitions are discussed in the
  \code{\link[evmix:kernels]{kernels}}, with the
  \code{lambda} as the default. The \code{bw} specification
  is the same as used in the
  \code{\link[stats:density]{density}} function.

  The possible kernels are also defined in
  \code{\link[evmix:kernels]{kernels}} with the
  \code{"gaussian"} as the default choice.

  The cumulative distribution function with tail fraction
  \eqn{\phi_u} defined by the upper tail fraction of the
  kernel density estimate (\code{phiu=TRUE}), upto the
  threshold \eqn{x \le u}, given by: \deqn{F(x) = H(x)} and
  above the threshold \eqn{x > u}: \deqn{F(x) = H(u) + [1 -
  H(u)] G(x)} where \eqn{H(x)} and \eqn{G(X)} are the KDE
  and conditional GPD cumulative distribution functions
  respectively.

  The cumulative distribution function for pre-specified
  \eqn{\phi_u}, upto the threshold \eqn{x \le u}, is given
  by: \deqn{F(x) = (1 - \phi_u) H(x)/H(u)} and above the
  threshold \eqn{x > u}: \deqn{F(x) = \phi_u + [1 - \phi_u]
  G(x)} Notice that these definitions are equivalent when
  \eqn{\phi_u = 1 - H(u)}.

  If no bandwidth is provided \code{lambda=NULL} and
  \code{bw=NULL} then the normal reference rule is used,
  using the \code{\link[stats:bandwidth]{bw.nrd0}}
  function, which is consistent with the
  \code{\link[stats:density]{density}} function. At least
  two kernel centres must be provided as the variance needs
  to be estimated.

  See \code{\link[evmix:gpd]{gpd}} for details of GPD upper
  tail component and \code{\link[evmix:kden]{dkden}} for
  details of KDE bulk component.
}
\note{
  Unlike all the other extreme value mixture model
  functions the \code{\link[evmix:kdengpd]{kdengpd}}
  functions have not been vectorised as this is not
  appropriate. The main inputs (\code{x}, \code{p} or
  \code{q}) must be either a scalar or a vector, which also
  define the output length. The \code{kerncentres} can also
  be a scalar or vector.

  The kernel centres \code{kerncentres} can either be a
  single datapoint or a vector of data. The kernel centres
  (\code{kerncentres}) and locations to evaluate density
  (\code{x}) and cumulative distribution function
  (\code{q}) would usually be different.

  Default values are provided for all inputs, except for
  the fundamentals \code{kerncentres}, \code{x}, \code{q}
  and \code{p}. The default sample size for
  \code{\link[evmix:kdengpd]{rkdengpd}} is 1.

  Missing (\code{NA}) and Not-a-Number (\code{NaN}) values
  in \code{x}, \code{p} and \code{q} are passed through as
  is and infinite values are set to \code{NA}. None of
  these are not permitted for the parameters or kernel
  centres.

  Due to symmetry, the lower tail can be described by GPD
  by negating the quantiles.

  Error checking of the inputs (e.g. invalid probabilities)
  is carried out and will either stop or give warning
  message as appropriate.
}
\examples{
\dontrun{
par(mfrow=c(2,2))
kerncentres=rnorm(500, 0, 1)
xx = seq(-4, 4, 0.01)
hist(kerncentres, breaks = 100, freq = FALSE)
lines(xx, dkdengpd(xx, kerncentres, u = 1.2, sigmau = 0.56, xi = 0.1))

plot(xx, pkdengpd(xx, kerncentres), type = "l")
lines(xx, pkdengpd(xx, kerncentres, xi = 0.3), col = "red")
lines(xx, pkdengpd(xx, kerncentres, xi = -0.3), col = "blue")
legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
      col=c("black", "red", "blue"), lty = 1, cex = 0.5)

x = rkdengpd(1000, kerncentres, phiu = 0.1, u = 1.2, sigmau = 0.56, xi = 0.1)
xx = seq(-4, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
lines(xx, dkdengpd(xx, kerncentres, phiu = 0.1, u = 1.2, sigmau = 0.56, xi = 0.1))

plot(xx, dkdengpd(xx, kerncentres, xi=0, phiu = 0.1), type = "l")
lines(xx, dkdengpd(xx, kerncentres, xi=0.2, phiu = 0.1), col = "red")
lines(xx, dkdengpd(xx, kerncentres, xi=-0.2, phiu = 0.1), col = "blue")
legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
      col=c("black", "red", "blue"), lty = 1)
}
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}. Based on code by
  Anna MacDonald produced for MATLAB.
}
\references{
  \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}

  \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  Bowman, A.W. (1984). An alternative method of
  cross-validation for the smoothing of density estimates.
  Biometrika 71(2), 353-360.

  Duin, R.P.W. (1976). On the choice of smoothing
  parameters for Parzen estimators of probability density
  functions. IEEE Transactions on Computers C25(11),
  1175-1179.

  MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B.,
  Reale, M. and Russell, G. (2011). A flexible extreme
  value mixture model. Computational Statistics and Data
  Analysis 55(6), 2137-2157.

  Wand, M. and Jones, M.C. (1995). Kernel Smoothing.
  Chapman && Hall.
}
\seealso{
  \code{\link[evmix:kernels]{kernels}},
  \code{\link[evmix:kfun]{kfun}},
  \code{\link[stats:density]{density}},
  \code{\link[stats:bandwidth]{bw.nrd0}} and
  \code{\link[ks:kde.1d]{dkde}} in
  \code{\link[ks:kde.1d]{ks}} package.

  Other kden kdengpd kdengpdcon bckden bckdengpd
  bckdengpdcon fkden fkdengpd fkdengpdcon fbckden
  fbckdengpd fbckdengpdcon: \code{\link{bckden}},
  \code{\link{bckdengpd}}, \code{\link{bckdengpdcon}},
  \code{\link{dbckden}}, \code{\link{dbckdengpd}},
  \code{\link{dbckdengpdcon}}, \code{\link{dkden}},
  \code{\link{dkdengpdcon}}, \code{\link{fbckden}},
  \code{\link{fkden}}, \code{\link{kden}},
  \code{\link{kdengpdcon}}, \code{\link{lbckden}},
  \code{\link{lkden}}, \code{\link{nlbckden}},
  \code{\link{nlkden}}, \code{\link{pbckden}},
  \code{\link{pbckdengpd}}, \code{\link{pbckdengpdcon}},
  \code{\link{pkden}}, \code{\link{pkdengpdcon}},
  \code{\link{qbckden}}, \code{\link{qbckdengpd}},
  \code{\link{qbckdengpdcon}}, \code{\link{qkden}},
  \code{\link{qkdengpdcon}}, \code{\link{rbckden}},
  \code{\link{rbckdengpd}}, \code{\link{rbckdengpdcon}},
  \code{\link{rkden}}, \code{\link{rkdengpdcon}}
}

