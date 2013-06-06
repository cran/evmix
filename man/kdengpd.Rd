\name{kdengpd}
\alias{dkdengpd}
\alias{kdengpd}
\alias{pkdengpd}
\alias{qkdengpd}
\alias{rkdengpd}
\title{Kernel Density Estimation Using Normal Kernel and GPD Tail Extreme Value Mixture Model}
\usage{
  dkdengpd(x, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, log = FALSE)

  pkdengpd(q, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, lower.tail = TRUE)

  qkdengpd(p, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE, lower.tail = TRUE)

  rkdengpd(n = 1, kerncentres, lambda = NULL,
    u = as.vector(quantile(kerncentres, 0.9)),
    sigmau = sqrt(6 * var(kerncentres))/pi, xi = 0,
    phiu = TRUE)
}
\arguments{
  \item{x}{quantile}

  \item{kerncentres}{kernel centres (typically sample
  data)}

  \item{lambda}{bandwidth for normal kernel (standard
  deviation of normal)}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantile}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probability}

  \item{n}{sample size (non-negative integer)}

  \item{u}{threshold}

  \item{sigmau}{scale parameter (non-negative)}

  \item{xi}{shape parameter}

  \item{phiu}{probability of being above threshold [0,1]}
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
  function and random number generation for the kernel
  density estimation using normal kernel for the bulk
  distribution upto the threshold and conditional GPD above
  threshold. The parameters are the bandwidth
  \code{lambda}, threshold \code{u} GPD scale \code{sigmau}
  and shape \code{xi} and tail fraction \code{phiu}.
}
\details{
  Extreme value mixture model combining kernel density
  estimation Using normal kernel for the bulk below the
  threshold and GPD for upper tail. The user can
  pre-specify \code{phiu} permitting a parameterised value
  for the tail fraction \eqn{\phi_u}. Alternatively, when
  \code{phiu=TRUE} the tail fraction is estimated as the
  tail fraction from the normal bulk model.

  The cumulative distribution function with tail fraction
  \eqn{\phi_u} defined by the upper tail fraction of the
  kernel density estimation using normal kernel
  (\code{phiu=TRUE}), upto the threshold \eqn{x \le u},
  given by: \deqn{F(x) = H(x)} and above the threshold
  \eqn{x > u}: \deqn{F(x) = H(u) + [1 - H(u)] G(x)} where
  \eqn{H(x)} and \eqn{G(X)} are the kernel density
  estimator and conditional GPD cumulative distribution
  functions (i.e. \code{mean(pnorm(x, kerncentres,
  lambda))} and \code{pgpd(x, u, sigmau, xi)}).

  The cumulative distribution function for pre-specified
  \eqn{\phi_u}, upto the threshold \eqn{x \le u}, is given
  by: \deqn{F(x) = (1 - \phi_u) H(x)/H(u)} and above the
  threshold \eqn{x > u}: \deqn{F(x) = \phi_u + [1 - \phi_u]
  G(x)} Notice that these definitions are equivalent when
  \eqn{\phi_u = 1 - H(u)}.

  See \code{\link[evmix:gpd]{gpd}} for details of GPD upper
  tail component and \code{\link[stats:Normal]{dnorm}} for
  details of normal bulk component.
}
\note{
  Unlike all the other extreme value mixture model
  functions the \code{\link[evmix:kdengpd]{kdengpd}}
  functions have not been vectorised as this is not
  appropriate. The main inputs (\code{x}, \code{p} or
  \code{q}) must be either a scalar or a vector, which also
  define the output length.  The \code{kerncentres} can
  also be a scalar or vector.

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
  in \code{x} and \code{q} are passed through as is and
  infinite values are set to \code{NA}.

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

kerncentres=rnorm(1000, 0, 1)
x = rkdengpd(1000, kerncentres, phiu = 0.1, u = 1.2, sigmau = 0.56, xi = 0.1)
xx = seq(-4, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
lines(xx, dkdengpd(xx, kerncentres, phiu = 0.1))

plot(xx, dkdengpd(xx, kerncentres, xi=0, phiu = 0.2), type = "l")
lines(xx, dkdengpd(xx, kerncentres, xi=-0.2, phiu = 0.2), col = "red")
lines(xx, dkdengpd(xx, kerncentres, xi=0.2, phiu = 0.2), col = "blue")
legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
      col=c("black", "red", "blue"), lty = 1)
}
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}
}
\references{
  \url{http://en.wikipedia.org/wiki/Normal_distribution}

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
}
\seealso{
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[stats:Normal]{dnorm}}

  Other kdengpd: \code{\link{fkdengpd}}
}

