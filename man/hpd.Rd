\name{hpd}
\alias{dhpd}
\alias{hpd}
\alias{phpd}
\alias{qhpd}
\alias{rhpd}
\title{Hybrid Pareto Extreme Value Mixture Model}
\usage{
  dhpd(x, nmean = 0, nsd = 1, xi = 0, log = FALSE)

  phpd(q, nmean = 0, nsd = 1, xi = 0, lower.tail = TRUE)

  qhpd(p, nmean = 0, nsd = 1, xi = 0, lower.tail = TRUE)

  rhpd(n = 1, nmean = 0, nsd = 1, xi = 0)
}
\arguments{
  \item{x}{quantile}

  \item{nmean}{normal mean}

  \item{nsd}{normal standard deviation (non-negative)}

  \item{xi}{shape parameter}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantile}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probability}

  \item{n}{sample size (non-negative integer)}
}
\value{
  \code{\link[evmix:hpd]{dhpd}} gives the density,
  \code{\link[evmix:hpd]{phpd}} gives the cumulative
  distribution function, \code{\link[evmix:hpd]{qhpd}}
  gives the quantile function and
  \code{\link[evmix:hpd]{rhpd}} gives a random sample.
}
\description{
  Density, cumulative distribution function, quantile
  function and random number generation for the hybrid
  Pareto extreme value mixture model.  The parameters are
  the normal mean \code{nmean} and standard deviation
  \code{nsd} and GPD shape \code{xi}.
}
\details{
  Extreme value mixture model combining normal distribution
  for the bulk below the threshold and GPD for upper tail
  which is continuous in its zeroth and first derivative at
  the threshold, but with one important difference to all
  the other mixture models.

  The hybrid Pareto does not include the usual tail
  fraction \code{phiu} scaling, i.e. so the GPD is not
  treated as a conditional model for the exceedances. The
  unscaled GPD is simply spliced with the normal truncated
  at the threshold, with no rescaling to account for the
  proportion above the threshold being applied. The
  parameters have to adjust for the lack of tail fraction
  scaling.

  The two continuity constraints lead to the threshold
  \code{u} and GPD scale \code{sigmau} being replaced by a
  function of the normal mean, standard deviation and GPD
  shape parameters.

  The continuity constraint on its first derivative at the
  threshold means that \eqn{ h'(u) = g'(u)}. Then the
  Lambert W function is used for replacing the threshold u
  and GPD scale sigmau in terms of the normal mean,
  standard deviation and GPD shape xi. The cumulative
  distribution function defined upto the threshold \eqn{x
  \le u}, given by: \deqn{F(x) = H(x) / r } and above the
  threshold \eqn{x > u}: \deqn{F(x) = (H(u) + G(x)) / r }
  where \eqn{H(x)} and \eqn{G(X)} are the normal and
  conditional GPD cumulative distribution functions (i.e.
  \code{pnorm(x, nmean, nsd)}, \code{pgpd(x, u, sigmau,
  xi)}). The normalisation constant \eqn{r} ensures a
  proper density and is given by\code{r = 1 + pnorm(u, mean
  = nmean, sd = nsd)}, i.e. the 1 comes from integration of
  the unscaled GPD and the second term is from the usual
  normal component.

  The continuity constraint on the density at the threshold
  means that \eqn{ h(u) = g(u)} where \eqn{h(x)} and
  \eqn{g(x)} are the normal and unscaled GPD density
  functions (i.e. \code{dnorm(u, nmean, nsd)} and
  \code{dgpd(u, u, sigmau, xi)}).

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
  element of vector. In the case of \code{rhpdcon} any
  input vector must be of length \code{n}.

  Default values are provided for all inputs, except for
  the fundamentals \code{x}, \code{q} and \code{p}. The
  default sample size for \code{\link[evmix:hpd]{rhpd}} is
  1.

  Missing (\code{NA}) and Not-a-Number (\code{NaN}) values
  in \code{x} and \code{q} are passed through as is and
  infinite values are set to \code{NA}.

  Error checking of the inputs (e.g. invalid probabilities)
  is carried out and will either stop or give warning
  message as appropriate.
}
\examples{
\dontrun{
par(mfrow = c(2, 2))
xx = seq(-5, 20, 0.01)
f1 = dhpd(xx, 0,1, 0.4)
plot(xx, f1, type = "l")
abline(v = 0.4942921)

# three tail behaviours
plot(xx, phpd(xx), type = "l")
lines(xx, phpd(xx, xi = 0.3), col = "red")
lines(xx, phpd(xx, xi = -0.3), col = "blue")
legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
  col=c("black", "red", "blue"), lty = 1)

sim = rhpd(10000, nmean = 0, nsd = 1.5, xi = 0.2)
hist(sim, freq = FALSE, 100, xlim = c(-5, 20), ylim = c(0, 0.2))
lines(xx, dhpd(xx, 0, 1.5, 0.2), col = "blue")

plot(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = 0), type = "l")
lines(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = 0.2), col = "red")
lines(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = -0.2), col = "blue")
legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
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

  Carreau, J. and Y. Bengio (2008). A hybrid Pareto model
  for asymmetric fat-tailed data: the univariate case.
  Extremes 12 (1), 53-76.
}
\seealso{
  \code{\link[evmix:gpd]{gpd}},
  \code{\link[stats:Normal]{dnorm}},
  \code{\link[evmix:normgpd]{dnormgpd}} and
  \code{\link[evmix:normgpdcon]{dnormgpdcon}}. The
  \code{\link[condmixt:condmixt-package]{condmixt}} package
  written by one of the original authors of the hybrid
  Pareto model (Carreau and Bengio, 2008) also has similar
  functions for the hybrid Pareto
  \code{\link[condmixt:hpareto]{hpareto}} and mixture of
  hybrid Paretos
  \code{\link[condmixt:hparetomixt]{hparetomixt}}, which
  are more flexible as they also permit the model to be
  truncated at zero.

  Other hpd: \code{\link{fhpd}}, \code{\link{lhpd}},
  \code{\link{nlhpd}}
}

