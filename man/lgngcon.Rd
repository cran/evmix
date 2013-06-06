\name{lgngcon}
\alias{lgngcon}
\alias{nlgngcon}
\title{Log-likelihood of Normal Bulk with GPD Upper and Lower Tails Extreme Value Mixture Model with Continuity Constraints}
\usage{
  lgngcon(x, nmean = 0, nsd = 1,
    ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE,
    ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE,
    log = TRUE)

  nlgngcon(pvector, x, phiul = TRUE, phiur = TRUE,
    finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{nmean}{normal mean}

  \item{nsd}{normal standard deviation (non-negative)}

  \item{ul}{lower tail threshold}

  \item{xil}{lower tail GPD shape parameter}

  \item{phiul}{probability of being above threshold (0, 1)
  or logical}

  \item{ur}{upper tail threshold}

  \item{xir}{upper tail GPD shape parameter}

  \item{phiur}{probability of being above threshold (0, 1)
  or logical}

  \item{log}{logical, if TRUE then log density}

  \item{pvector}{vector of initial values of mixture model
  parameters or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}
}
\value{
  \code{\link[evmix:lgngcon]{lgngcon}} gives
  (log-)likelihood and
  \code{\link[evmix:lgngcon]{nlgngcon}} gives the negative
  log-likelihood.
}
\description{
  Log-likelihood and negative log-likelihood for the
  extreme value mixture model with normal for bulk
  distribution between the lower and upper thresholds with
  conditional GPD for the two tails with continuity
  constraints.
}
\details{
  The likelihood functions for the extreme value mixture
  model with normal bulk and GPD for the two tails, as used
  in the maximum likelihood fitting function
  \code{\link[evmix:fgngcon]{fgngcon}}.

  They are designed to be used for MLE in
  \code{\link[evmix:fgngcon]{fgngcon}} but are available
  for wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fgngcon]{fgngcon}},
  \code{\link[evmix:gngcon]{gngcon}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:lgngcon]{lgngcon}}, which takes
  parameters as inputs in the same form as distribution
  functions. The negative log-likelihood is a wrapper for
  \code{\link[evmix:lgngcon]{lgngcon}}, designed towards
  making it useable for optimisation (e.g. parameters are
  given a vector as first input). The tail fractions
  \code{phiul} and \code{phiur} are treated separately to
  the other parameters, to allow for all it's
  representations.

  Unlike the distribution functions
  \code{\link[evmix:gngcon]{gngcon}} the \code{phiu} must
  be either logical (\code{TRUE} or \code{FALSE}) or
  numerical in range \eqn{(0, 1)}. The default is to
  specify \code{phiu=TRUE} so that the tail fraction is
  specified by normal distribution \eqn{\phi_u = 1 - H(u)},
  or \code{phiu=FALSE} to treat the tail fraction as an
  extra parameter estimated using the sample proportion.
  Specify a numeric \code{phiu} as pre-specified
  probability \eqn{(0, 1)}. Notice that the tail fraction
  probability cannot be 0 or 1 otherwise there would be no
  contribution from either tail or bulk components
  respectively.

  The function \code{\link[evmix:lgngcon]{lgngcon}} carries
  out the calculations for the log-likelihood directly,
  which can be exponentiated to give actual likelihood
  using (\code{log=FALSE}).
}
\note{
  Unlike all the distribution functions for this mixture
  model, the likelihood functions only permits a scalar
  value for all the parameters. Only the data is a vector.

  Default values for the tail fractions \code{phiul=TRUE}
  and \code{phiur=TRUE} is given in both
  \code{\link[evmix:lgngcon]{lgngcon}} and
  \code{\link[evmix:lgngcon]{nlgngcon}}. The
  \code{\link[evmix:lgngcon]{lgngcon}} also has the usual
  defaults for the other parameters, but
  \code{\link[evmix:lgngcon]{nlgngcon}} has no defaults.

  Invalid parameters will give \code{0} for likelihood,
  \code{log(0)=-Inf} for log-likelihood and
  \code{-log(0)=Inf} for negative log-likelihood.

  See \code{\link[evmix:gpd]{fgpd}} for explanation of
  \code{finitelik}.

  Error checking of the inputs is carried out and will
  either stop or give warning message as appropriate.
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
  \code{\link[evmix:lgng]{lgng}},
  \code{\link[evmix:lnormgpd]{lnormgpd}},
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other gngcon: \code{\link{dgngcon}},
  \code{\link{fgngcon}}, \code{\link{gngcon}},
  \code{\link{pgngcon}}, \code{\link{qgngcon}},
  \code{\link{rgngcon}}
}

