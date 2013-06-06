\name{lognormgpdcon}
\alias{dlognormgpdcon}
\alias{lognormgpdcon}
\alias{plognormgpdcon}
\alias{qlognormgpdcon}
\alias{rlognormgpdcon}
\title{Log-Normal Bulk and GPD Tail Extreme Value Mixture Model with Continuity Constraint}
\usage{
  dlognormgpdcon(x, lnmean = 0, lnsd = 1,
    u = qlnorm(0.9, lnmean, lnsd), xi = 0, phiu = TRUE,
    log = FALSE)

  plognormgpdcon(q, lnmean = 0, lnsd = 1,
    u = qlnorm(0.9, lnmean, lnsd), xi = 0, phiu = TRUE,
    lower.tail = TRUE)

  qlognormgpdcon(p, lnmean = 0, lnsd = 1,
    u = qlnorm(0.9, lnmean, lnsd), xi = 0, phiu = TRUE,
    lower.tail = TRUE)

  rlognormgpdcon(n = 1, lnmean = 0, lnsd = 1,
    u = qlnorm(0.9, lnmean, lnsd), xi = 0, phiu = TRUE)
}
\arguments{
  \item{x}{quantile}

  \item{lnmean}{mean on log scale}

  \item{lnsd}{standard deviation on log scale
  (non-negative)}

  \item{u}{threshold}

  \item{xi}{shape parameter}

  \item{phiu}{probability of being above threshold [0,1] or
  TRUE}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantile}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probability}

  \item{n}{sample size (non-negative integer)}
}
\value{
  \code{\link[evmix:lognormgpdcon]{dlognormgpdcon}} gives
  the density,
  \code{\link[evmix:lognormgpdcon]{plognormgpdcon}} gives
  the cumulative distribution function,
  \code{\link[evmix:lognormgpdcon]{qlognormgpdcon}} gives
  the quantile function and
  \code{\link[evmix:lognormgpdcon]{rlognormgpdcon}} gives a
  random sample.
}
\description{
  Density, cumulative distribution function, quantile
  function and random number generation for the extreme
  value mixture model with log-normal for bulk distribution
  upto the threshold and conditional GPD above threshold
  with a continuity constraint. The parameters are the
  normal mean \code{lnmean} and standard deviation
  \code{lnsd}, threshold \code{u} GPD scale \code{sigmau}
  and shape \code{xi} and tail fraction \code{phiu}.
}
\details{
  Extreme value mixture model combining log-normal
  distribution for the bulk below the threshold and GPD for
  upper tail with a continuity constraint. The user can
  pre-specify \code{phiu} permitting a parameterised value
  for the tail fraction \eqn{\phi_u}. Alternatively, when
  \code{phiu=TRUE} the tail fraction is estimated as the
  tail fraction from the log-normal bulk model.

  The cumulative distribution function with tail fraction
  \eqn{\phi_u} defined by the upper tail fraction of the
  log-normal bulk model (\code{phiu=TRUE}), upto the
  threshold \eqn{0 < x \le u}, given by: \deqn{F(x) = H(x)}
  and above the threshold \eqn{x > u}: \deqn{F(x) = H(u) +
  [1 - H(u)] G(x)} where \eqn{H(x)} and \eqn{G(X)} are the
  log-normal and conditional GPD cumulative distribution
  functions (i.e. \code{plnorm(x, meanlog = lnmean, sdlog =
  lnsd)} and \code{pgpd(x, u, sigmau, xi)}).

  The cumulative distribution function for pre-specified
  \eqn{\phi_u}, upto the threshold \eqn{0 < x \le u}, is
  given by: \deqn{F(x) = (1 - \phi_u) H(x)/H(u)} and above
  the threshold \eqn{x > u}: \deqn{F(x) = \phi_u + [1 -
  \phi_u] G(x)} Notice that these definitions are
  equivalent when \eqn{\phi_u = 1 - H(u)}.

  The continuity constraint means that \eqn{(1 - \phi_u)
  h(u)/H(u) = \phi_u g(u)} where \eqn{h(x)} and \eqn{g(x)}
  are the log-normal and conditional GPD density functions
  (i.e. \code{dlnorm(x, nmean, nsd)} and \code{dgpd(x, u,
  sigmau, xi)}). The resulting GPD scale parameter is then:
  \deqn{\sigma_u = \phi_u H(u) / [1 - \phi_u] h(u)}. In the
  special case of where the tail fraction is defined by the
  bulk model this reduces to \deqn{\sigma_u = [1 - H(u)] /
  h(u)}.

  The gamma is defined on the non-negative reals, so the
  threshold must be non-negative.

  See \code{\link[evmix:gpd]{gpd}} for details of GPD upper
  tail component and \code{\link[stats:Lognormal]{dlnorm}}
  for details of log-normal bulk component.
}
\note{
  All inputs are vectorised except \code{log} and
  \code{lower.tail}. The main inputs (\code{x}, \code{p} or
  \code{q}) and parameters must be either a scalar or a
  vector. If vectors are provided they must all be of the
  same length, and the function will be evaluated for each
  element of vector. In the case of \code{rlognormgpd} any
  input vector must be of length \code{n}.

  Default values are provided for all inputs, except for
  the fundamentals \code{x}, \code{q} and \code{p}. The
  default sample size for
  \code{\link[evmix:lognormgpdcon]{rlognormgpdcon}} is 1.

  Missing (\code{NA}) and Not-a-Number (\code{NaN}) values
  in \code{x} and \code{q} are passed through as is and
  infinite values are set to \code{NA}.

  Error checking of the inputs (e.g. invalid probabilities)
  is carried out and will either stop or give warning
  message as appropriate.
}
\examples{
\dontrun{
par(mfrow=c(2,2))
x = rlognormgpdcon(1000)
xx = seq(-1, 10, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
lines(xx, dlognormgpdcon(xx))

# three tail behaviours
plot(xx, plognormgpdcon(xx), type = "l")
lines(xx, plognormgpdcon(xx, xi = 0.3), col = "red")
lines(xx, plognormgpdcon(xx, xi = -0.3), col = "blue")
legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
  col=c("black", "red", "blue"), lty = 1)

x = rlognormgpdcon(1000, u = 2, phiu = 0.2)
hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
lines(xx, dlognormgpdcon(xx, u = 2, phiu = 0.2))

plot(xx, dlognormgpdcon(xx, u = 2, xi=0, phiu = 0.2), type = "l")
lines(xx, dlognormgpdcon(xx, u = 2, xi=-0.2, phiu = 0.2), col = "red")
lines(xx, dlognormgpdcon(xx, u = 2, xi=0.2, phiu = 0.2), col = "blue")
legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
  col=c("black", "red", "blue"), lty = 1)
  }
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}
}
\references{
  \url{http://en.wikipedia.org/wiki/Log-normal_distribution}

  \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  Solari, S. and Losada, M.A. (2004). A unified statistical
  model for hydrological variables including the selection
  of threshold for the peak over threshold method. Water
  Resources Research. 48, W10541.
}
\seealso{
  \code{\link[evmix:lognormgpd]{lognormgpd}},
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[stats:Lognormal]{dlnorm}}

  Other lognormgpdcon: \code{\link{flognormgpdcon}},
  \code{\link{llognormgpdcon}},
  \code{\link{nllognormgpdcon}}
}

