% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnormgpd.r
\name{normgpd}
\alias{normgpd}
\alias{dnormgpd}
\alias{pnormgpd}
\alias{qnormgpd}
\alias{rnormgpd}
\title{Normal Bulk and GPD Tail Extreme Value Mixture Model}
\usage{
dnormgpd(x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  sigmau = nsd, xi = 0, phiu = TRUE, log = FALSE)

pnormgpd(q, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  sigmau = nsd, xi = 0, phiu = TRUE, lower.tail = TRUE)

qnormgpd(p, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  sigmau = nsd, xi = 0, phiu = TRUE, lower.tail = TRUE)

rnormgpd(n = 1, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  sigmau = nsd, xi = 0, phiu = TRUE)
}
\arguments{
\item{x}{quantiles}

\item{nmean}{normal mean}

\item{nsd}{normal standard deviation (positive)}

\item{u}{threshold}

\item{sigmau}{scale parameter (positive)}

\item{xi}{shape parameter}

\item{phiu}{probability of being above threshold \eqn{[0, 1]} or \code{TRUE}}

\item{log}{logical, if TRUE then log density}

\item{q}{quantiles}

\item{lower.tail}{logical, if FALSE then upper tail probabilities}

\item{p}{cumulative probabilities}

\item{n}{sample size (positive integer)}
}
\value{
\code{\link[evmix:normgpd]{dnormgpd}} gives the density, 
\code{\link[evmix:normgpd]{pnormgpd}} gives the cumulative distribution function,
\code{\link[evmix:normgpd]{qnormgpd}} gives the quantile function and 
\code{\link[evmix:normgpd]{rnormgpd}} gives a random sample.
}
\description{
Density, cumulative distribution function, quantile function and
  random number generation for the extreme value mixture model with normal for bulk
  distribution upto the threshold and conditional GPD above threshold. The parameters
  are the normal mean \code{nmean} and standard deviation \code{nsd}, threshold \code{u}
  GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
}
\details{
Extreme value mixture model combining normal distribution for the bulk
below the threshold and GPD for upper tail.

The user can pre-specify \code{phiu} 
permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
\code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
normal bulk model.

The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
upper tail fraction of the normal bulk model (\code{phiu=TRUE}), upto the 
threshold \eqn{x \le u}, given by:
\deqn{F(x) = H(x)}
and above the threshold \eqn{x > u}:
\deqn{F(x) = H(u) + [1 - H(u)] G(x)}
where \eqn{H(x)} and \eqn{G(X)} are the normal and conditional GPD
cumulative distribution functions (i.e. \code{pnorm(x, nmean, nsd)} and
\code{pgpd(x, u, sigmau, xi)}) respectively.

The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
threshold \eqn{x \le u}, is given by:
\deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
and above the threshold \eqn{x > u}:
\deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.

See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
\code{\link[stats:Normal]{dnorm}} for details of normal bulk component.
}
\note{
All inputs are vectorised except \code{log} and \code{lower.tail}.
The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
a scalar or a vector. If vectors are provided they must all be of the same length,
and the function will be evaluated for each element of vector. In the case of 
\code{\link[evmix:normgpd]{rnormgpd}} any input vector must be of length \code{n}.

Default values are provided for all inputs, except for the fundamentals 
\code{x}, \code{q} and \code{p}. The default sample size for 
\code{\link[evmix:normgpd]{rnormgpd}} is 1.

Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
\code{p} and \code{q} are passed through as is and infinite values are set to
\code{NA}. None of these are not permitted for the parameters.

Due to symmetry, the lower tail can be described by GPD by negating the quantiles. 
The normal mean \code{nmean} and GPD threshold \code{u} will also require negation.

Error checking of the inputs (e.g. invalid probabilities) is carried out and
will either stop or give warning message as appropriate.
}
\examples{
\dontrun{
set.seed(1)
par(mfrow = c(2, 2))

x = rnormgpd(1000)
xx = seq(-4, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
lines(xx, dnormgpd(xx))

# three tail behaviours
plot(xx, pnormgpd(xx), type = "l")
lines(xx, pnormgpd(xx, xi = 0.3), col = "red")
lines(xx, pnormgpd(xx, xi = -0.3), col = "blue")
legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
  col=c("black", "red", "blue"), lty = 1)

x = rnormgpd(1000, phiu = 0.2)
xx = seq(-4, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
lines(xx, dnormgpd(xx, phiu = 0.2))

plot(xx, dnormgpd(xx, xi=0, phiu = 0.2), type = "l")
lines(xx, dnormgpd(xx, xi=-0.2, phiu = 0.2), col = "red")
lines(xx, dnormgpd(xx, xi=0.2, phiu = 0.2), col = "blue")
legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
  col=c("black", "red", "blue"), lty = 1)
}

}
\references{
\url{http://en.wikipedia.org/wiki/Normal_distribution}

\url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
threshold estimation and uncertainty quantification. REVSTAT - Statistical
Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

Hu Y. and Scarrott, C.J. (2018). evmix: An R Package for Extreme Value Mixture Modeling, 
Threshold Estimation and Boundary Corrected Kernel Density Estimation. Journal of
Statistical Software 84(5), 1-27. doi: 10.18637/jss.v084.i05.

Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004). Bayesian analysis of extreme
events with threshold estimation. Statistical Modelling. 4(3), 227-244.
}
\seealso{
\code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Normal]{dnorm}}

Other normgpd: \code{\link{fgng}}, \code{\link{fhpd}},
  \code{\link{fitmnormgpd}}, \code{\link{flognormgpd}},
  \code{\link{fnormgpdcon}}, \code{\link{fnormgpd}},
  \code{\link{gngcon}}, \code{\link{gng}},
  \code{\link{hpdcon}}, \code{\link{hpd}},
  \code{\link{itmnormgpd}}, \code{\link{lognormgpdcon}},
  \code{\link{lognormgpd}}, \code{\link{normgpdcon}}

Other normgpdcon: \code{\link{fgngcon}},
  \code{\link{fhpdcon}}, \code{\link{flognormgpdcon}},
  \code{\link{fnormgpdcon}}, \code{\link{fnormgpd}},
  \code{\link{gngcon}}, \code{\link{gng}},
  \code{\link{hpdcon}}, \code{\link{hpd}},
  \code{\link{normgpdcon}}

Other gng: \code{\link{fgngcon}}, \code{\link{fgng}},
  \code{\link{fitmgng}}, \code{\link{fnormgpd}},
  \code{\link{gngcon}}, \code{\link{gng}},
  \code{\link{itmgng}}

Other fnormgpd: \code{\link{fnormgpd}}
}
\author{
Yang Hu and Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}
}
\concept{fnormgpd}
\concept{gng}
\concept{normgpd}
\concept{normgpdcon}
