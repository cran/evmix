\name{normgpdcon}
\alias{dnormgpdcon}
\alias{normgpdcon}
\alias{pnormgpdcon}
\alias{qnormgpdcon}
\alias{rnormgpdcon}
\title{Normal Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint}
\usage{
  dnormgpdcon(x, nmean = 0, nsd = 1,
    u = qnorm(0.9, nmean, nsd), xi = 0, phiu = TRUE,
    log = FALSE)

  pnormgpdcon(q, nmean = 0, nsd = 1,
    u = qnorm(0.9, nmean, nsd), xi = 0, phiu = TRUE,
    lower.tail = TRUE)

  qnormgpdcon(p, nmean = 0, nsd = 1,
    u = qnorm(0.9, nmean, nsd), xi = 0, phiu = TRUE,
    lower.tail = TRUE)

  rnormgpdcon(n = 1, nmean = 0, nsd = 1,
    u = qnorm(0.9, nmean, nsd), xi = 0, phiu = TRUE)
}
\arguments{
  \item{x}{quantiles}

  \item{nmean}{normal mean}

  \item{nsd}{normal standard deviation (positive)}

  \item{u}{threshold}

  \item{xi}{shape parameter}

  \item{phiu}{probability of being above threshold \eqn{[0,
  1]} or \code{TRUE}}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantiles}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probabilities}

  \item{n}{sample size (positive integer)}
}
\value{
  \code{\link[evmix:normgpdcon]{dnormgpdcon}} gives the
  density, \code{\link[evmix:normgpdcon]{pnormgpdcon}}
  gives the cumulative distribution function,
  \code{\link[evmix:normgpdcon]{qnormgpdcon}} gives the
  quantile function and
  \code{\link[evmix:normgpdcon]{rnormgpdcon}} gives a
  random sample.
}
\description{
  Density, cumulative distribution function, quantile
  function and random number generation for the extreme
  value mixture model with normal for bulk distribution
  upto the threshold and conditional GPD above threshold
  with continuity at threshold. The parameters are the
  normal mean \code{nmean} and standard deviation
  \code{nsd}, threshold \code{u} and GPD shape \code{xi}
  and tail fraction \code{phiu}.
}
\details{
  Extreme value mixture model combining normal distribution
  for the bulk below the threshold and GPD for upper tail
  with continuity at threshold.

  The user can pre-specify \code{phiu} permitting a
  parameterised value for the tail fraction \eqn{\phi_u}.
  Alternatively, when \code{phiu=TRUE} the tail fraction is
  estimated as the tail fraction from the normal bulk
  model.

  The cumulative distribution function with tail fraction
  \eqn{\phi_u} defined by the upper tail fraction of the
  normal bulk model (\code{phiu=TRUE}), upto the threshold
  \eqn{x \le u}, given by: \deqn{F(x) = H(x)} and above the
  threshold \eqn{x > u}: \deqn{F(x) = H(u) + [1 - H(u)]
  G(x)} where \eqn{H(x)} and \eqn{G(X)} are the normal and
  conditional GPD cumulative distribution functions (i.e.
  \code{pnorm(x, nmean, nsd)} and \code{pgpd(x, u, sigmau,
  xi)}) respectively.

  The cumulative distribution function for pre-specified
  \eqn{\phi_u}, upto the threshold \eqn{x \le u}, is given
  by: \deqn{F(x) = (1 - \phi_u) H(x)/H(u)} and above the
  threshold \eqn{x > u}: \deqn{F(x) = \phi_u + [1 - \phi_u]
  G(x)} Notice that these definitions are equivalent when
  \eqn{\phi_u = 1 - H(u)}.

  The continuity constraint means that \eqn{(1 - \phi_u)
  h(u)/H(u) = \phi_u g(u)} where \eqn{h(x)} and \eqn{g(x)}
  are the normal and conditional GPD density functions
  (i.e. \code{dnorm(x, nmean, nsd)} and \code{dgpd(x, u,
  sigmau, xi)}) respectively. The resulting GPD scale
  parameter is then: \deqn{\sigma_u = \phi_u H(u) / [1 -
  \phi_u] h(u)}. In the special case of where the tail
  fraction is defined by the bulk model this reduces to
  \deqn{\sigma_u = [1 - H(u)] / h(u)}.

  See \code{\link[evmix:gpd]{gpd}} for details of GPD upper
  tail component and \code{\link[stats:Normal]{dnorm}} for
  details of normal bulk component.
}
\note{
  All inputs are vectorised except \code{log} and
  \code{lower.tail}. The main inputs (\code{x}, \code{p} or
  \code{q}) and parameters must be either a scalar or a
  vector. If vectors are provided they must all be of the
  same length, and the function will be evaluated for each
  element of vector. In the case of
  \code{\link[evmix:normgpdcon]{rnormgpdcon}} any input
  vector must be of length \code{n}.

  Default values are provided for all inputs, except for
  the fundamentals \code{x}, \code{q} and \code{p}. The
  default sample size for
  \code{\link[evmix:normgpdcon]{rnormgpdcon}} is 1.

  Missing (\code{NA}) and Not-a-Number (\code{NaN}) values
  in \code{x}, \code{p} and \code{q} are passed through as
  is and infinite values are set to \code{NA}. None of
  these are not permitted for the parameters.

  Due to symmetry, the lower tail can be described by GPD
  by negating the quantiles. The normal mean \code{nmean}
  and GPD threshold \code{u} will also require negation.

  Error checking of the inputs (e.g. invalid probabilities)
  is carried out and will either stop or give warning
  message as appropriate.
}
\examples{
\dontrun{
par(mfrow=c(2,2))
x = rnormgpdcon(1000)
xx = seq(-4, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
lines(xx, dnormgpdcon(xx))

# three tail behaviours
plot(xx, pnormgpdcon(xx), type = "l")
lines(xx, pnormgpdcon(xx, xi = 0.3), col = "red")
lines(xx, pnormgpdcon(xx, xi = -0.3), col = "blue")
legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
  col=c("black", "red", "blue"), lty = 1)

x = rnormgpdcon(1000, phiu = 0.2)
xx = seq(-4, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
lines(xx, dnormgpdcon(xx, phiu = 0.2))

plot(xx, dnormgpdcon(xx, xi=0, phiu = 0.2), type = "l")
lines(xx, dnormgpdcon(xx, xi=-0.2, phiu = 0.2), col = "red")
lines(xx, dnormgpdcon(xx, xi=0.2, phiu = 0.2), col = "blue")
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

  Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004).
  Bayesian analysis of extreme events with threshold
  estimation. Statistical Modelling. 4(3), 227-244.
}
\seealso{
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[stats:Normal]{dnorm}}

  Other normgpd normgpdcon gng gngcon fnormgpd fnormgpdcon
  fgng fgngcon: \code{\link{dgng}}, \code{\link{dgngcon}},
  \code{\link{dnormgpd}}, \code{\link{fgng}},
  \code{\link{fgngcon}}, \code{\link{fnormgpd}},
  \code{\link{fnormgpdcon}}, \code{\link{gng}},
  \code{\link{gngcon}}, \code{\link{lgng}},
  \code{\link{lgngcon}}, \code{\link{lnormgpd}},
  \code{\link{lnormgpdcon}}, \code{\link{nlgng}},
  \code{\link{nlgngcon}}, \code{\link{nlnormgpd}},
  \code{\link{nlnormgpdcon}}, \code{\link{nlugng}},
  \code{\link{nlugngcon}}, \code{\link{nlunormgpd}},
  \code{\link{nlunormgpdcon}}, \code{\link{normgpd}},
  \code{\link{pgng}}, \code{\link{pgngcon}},
  \code{\link{pnormgpd}}, \code{\link{proflugng}},
  \code{\link{proflugngcon}}, \code{\link{proflunormgpd}},
  \code{\link{proflunormgpdcon}}, \code{\link{qgng}},
  \code{\link{qgngcon}}, \code{\link{qnormgpd}},
  \code{\link{rgng}}, \code{\link{rgngcon}},
  \code{\link{rnormgpd}}
}

