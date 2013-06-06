\name{gngcon}
\alias{dgngcon}
\alias{gngcon}
\alias{pgngcon}
\alias{qgngcon}
\alias{rgngcon}
\title{Normal Bulk with GPD Upper and Lower Tails Extreme Value Mixture Model with Continuity Constraints}
\usage{
  dgngcon(x, nmean = 0, nsd = 1,
    ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE,
    ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE,
    log = FALSE)

  pgngcon(q, nmean = 0, nsd = 1,
    ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE,
    ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE,
    lower.tail = TRUE)

  qgngcon(p, nmean = 0, nsd = 1,
    ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE,
    ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE,
    lower.tail = TRUE)

  rgngcon(n = 1, nmean = 0, nsd = 1,
    ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE,
    ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE)
}
\arguments{
  \item{x}{quantile}

  \item{nmean}{normal mean}

  \item{nsd}{normal standard deviation (non-negative)}

  \item{ul}{lower tail threshold}

  \item{xil}{lower tail GPD shape parameter}

  \item{phiul}{probability of being below lower threshold
  (0,1)}

  \item{ur}{upper tail threshold}

  \item{xir}{upper tail GPD shape parameter}

  \item{phiur}{probability of being above upper threshold
  (0,1)}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantile}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probability}

  \item{n}{sample size (non-negative integer)}
}
\value{
  \code{\link[evmix:gngcon]{dgngcon}} gives the density,
  \code{\link[evmix:gngcon]{pgngcon}} gives the cumulative
  distribution function,
  \code{\link[evmix:gngcon]{qgngcon}} gives the quantile
  function and \code{\link[evmix:gngcon]{rgngcon}} gives a
  random sample.
}
\description{
  Density, cumulative distribution function, quantile
  function and random number generation for the extreme
  value mixture model with normal for bulk distribution
  between the upper and lower thresholds with conditional
  GPD's for the two tails with Continuity Constraints at
  the lower and upper threshold. The parameters are the
  normal mean \code{nmean} and standard deviation
  \code{nsd}, lower tail (threshold \code{ul}, GPD shape
  \code{xil} and tail fraction \code{phiul}) and upper tail
  (threshold \code{ur}, GPD shape \code{xiR} and tail
  fraction \code{phiuR}).
}
\details{
  Extreme value mixture model combining normal distribution
  for the bulk between the lower and upper thresholds and
  GPD for upper and lower tails with Continuity Constraints
  at the lower and upper threshold. The user can
  pre-specify \code{phiul} and \code{phiur} permitting a
  parameterised value for the lower and upper tail fraction
  respectively. Alternatively, when \code{phiul=TRUE} or
  \code{phiur=TRUE} the corresponding tail fraction is
  estimated as from the normal bulk model.

  Notice that the tail fraction cannot be 0 or 1, and the
  sum of upper and lower tail fractions
  \code{phiul+phiur<1}, so the lower threshold must be less
  than the upper, \code{ul<ur}.

  The cumulative distribution function now has three
  components. The lower tail with tail fraction
  \eqn{\phi_{ul}} defined by the normal bulk model
  (\code{phiul=TRUE}) upto the lower threshold \eqn{x <
  u_l}: \deqn{F(x) = H(u_l) G_l(x).} where \eqn{H(x)} is
  the normal cumulative distribution function (i.e.
  \code{pnorm(ur, nmean, nsd)}). The \eqn{G_l(X)} is the
  conditional GPD cumulative distribution function with
  negated data and threshold, i.e. \code{dgpd(-x, -ul,
  sigmaul, xil, phiul)}. The normal bulk model between the
  thresholds \eqn{u_l \le x \le u_r} given by: \deqn{F(x) =
  H(x).} Above the threshold \eqn{x > u_r} the usual
  conditional GPD: \deqn{F(x) = H(u_r) + [1 - H(u_r)] G(x)}
  where \eqn{G(X)}.

  The cumulative distribution function for the
  pre-specified tail fractions \eqn{\phi_{ul}} and
  \eqn{\phi_{ur}} is more complicated.  The unconditional
  GPD is used for the lower tail \eqn{x < u_l}: \deqn{F(x)
  = \phi_{ul} G_l(x).} The normal bulk model between the
  thresholds \eqn{u_l \le x \le u_r} given by: \deqn{F(x) =
  \phi_{ul}+ (1-\phi_{ul}-\phi_{ur}) (H(x) - H(u_l)) /
  (H(u_r) - H(u_l)).} Above the threshold \eqn{x > u_r} the
  usual conditional GPD: \deqn{F(x) = (1-\phi_{ur}) +
  \phi_{ur} G(x)} Notice that these definitions are
  equivalent when \eqn{\phi_{ul} = H(u_l)} and
  \eqn{\phi_{ur} = 1 - H(u_r)}.

  The continuity constraint at the \code{ul} means that:
  \deqn{\phi_{ul} g_l(x) = (1-\phi_{ul}-\phi_{ur}) h(u_l)/
  (H(u_r) - H(u_l)).} The GPD scale parameter
  \code{sigmaul} is replaced by: \deqn{\sigma_ul =
  \phi_{ul} (H(u_r) - H(u_l))/ h(u_l)
  (1-\phi_{ul}-\phi_{ur}).} In the special case of where
  the tail fraction is defined by the bulk model this
  reduces to \deqn{\sigma_ul = H(u_l)/ h(u_l)}.

  The continuity Constraint at the \code{ur} means that:
  \deqn{\phi_{ur} g_r(x) = (1-\phi_{ul}-\phi_{ur}) h(u_l)/
  (H(u_r) - H(u_l)).} By rearrangement, the GPD scale
  parameter \code{sigmaur} is then: \deqn{\sigma_ur =
  \phi_{ur} (H(u_r) - H(u_l))/ h(u_l)
  (1-\phi_{ul}-\phi_{ur}).} where \eqn{h(x)}, \eqn{g_l(x)}
  and \eqn{g_r(x)} are the normal and conditional GPD
  density functions for lower and upper tail (i.e.
  \code{dnorm(x, nmean, nsd)}, \code{dgpd(-x, -ul, sigmaul,
  xil, phiul)}, and \code{dgpd(x, ur, sigmaur, xir,
  phiur)}) In the special case of where the tail fraction
  is defined by the bulk model this reduces to
  \deqn{\sigma_ur = [1-H(u_r)] / h(u_r)}.

  See \code{\link[evmix:gpd]{gpd}} for details of GPD upper
  tail component, \code{\link[stats:Normal]{dnorm}} for
  details of normal bulk component,
  \code{\link[evmix:normgpd]{dnormgpd}} for normal with GPD
  extreme value mixture model and
  \code{\link[evmix:gng]{dgng}} for normal bulk with GPD
  upper and lower tails extreme value mixture model.
}
\note{
  All inputs are vectorised except \code{log} and
  \code{lower.tail}. The main inputs (\code{x}, \code{p} or
  \code{q}) and parameters must be either a scalar or a
  vector. If vectors are provided they must all be of the
  same length, and the function will be evaluated for each
  element of vector. In the case of \code{rgngcon} any
  input vector must be of length \code{n}.

  Default values are provided for all inputs, except for
  the fundamentals \code{x}, \code{q} and \code{p}. The
  default sample size for
  \code{\link[evmix:gngcon]{rgngcon}} is 1.

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
x = rgngcon(1000, phiul = 0.15, phiur = 0.15)
xx = seq(-6, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
lines(xx, dgngcon(xx, phiul = 0.15, phiur = 0.15))

# three tail behaviours
plot(xx, pgngcon(xx), type = "l")
lines(xx, pgngcon(xx, xil = 0.3, xir = 0.3), col = "red")
lines(xx, pgngcon(xx, xil = -0.3, xir = -0.3), col = "blue")
legend("topleft", paste("Symmetric xil=xir=",c(0, 0.3, -0.3)),
  col=c("black", "red", "blue"), lty = 1)

x = rgngcon(1000, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2)
xx = seq(-6, 6, 0.01)
hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
lines(xx, dgngcon(xx, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2))

plot(xx, dgngcon(xx, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2), type = "l", ylim = c(0, 0.4))
lines(xx, dgngcon(xx, xil = -0.3, phiul = 0.3, xir = 0.3, phiur = 0.3), col = "red")
lines(xx, dgngcon(xx, xil = -0.3, phiul = TRUE, xir = 0.3, phiur = TRUE), col = "blue")
legend("topleft", c("phiul = phiur = 0.2", "phiul = phiur = 0.3", "Bulk Tail Fraction"),
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

  Zhao, X., Scarrott, C.J. Reale, M. and Oxley, L. (2010).
  Extreme value modelling for forecasting the market
  crisis. Applied Financial Econometrics 20(1), 63-72.
}
\seealso{
  \code{\link[evmix:gng]{gng}},
  \code{\link[evmix:normgpd]{normgpd}},
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[stats:Normal]{dnorm}}

  Other gngcon: \code{\link{fgngcon}},
  \code{\link{lgngcon}}, \code{\link{nlgngcon}}
}

