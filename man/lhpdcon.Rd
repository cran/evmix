\name{lhpdcon}
\alias{lhpdcon}
\alias{nlhpdcon}
\title{Log-likelihood of Hybrid Pareto Extreme Value Mixture Model with Single Continuity Constraint}
\usage{
  lhpdcon(x, nmean = 0, nsd = 1,
    u = qnorm(0.9, nmean, nsd), xi = 0, log = TRUE)

  nlhpdcon(pvector, x, finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{nmean}{normal mean}

  \item{nsd}{normal standard deviation (non-negative)}

  \item{u}{threshold}

  \item{xi}{shape parameter}

  \item{log}{logical, if TRUE then log density}

  \item{pvector}{vector of initial values of mixture model
  parameters (\code{nmean}, \code{nsd}, \code{u},
  \code{sigmau}, \code{xi}) or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}
}
\value{
  \code{\link[evmix:lhpdcon]{lhpdcon}} gives
  (log-)likelihood and
  \code{\link[evmix:lhpdcon]{nlhpdcon}} gives the negative
  log-likelihood.
}
\description{
  Log-likelihood and negative log-likelihood for the hybrid
  Pareto extreme value mixture model with a single
  continuity constraint
}
\details{
  The likelihood functions for hybrid Pareto extreme value
  mixture model with single continuity constraint, as used
  in the maximum likelihood fitting function
  \code{\link[evmix:fhpdcon]{fhpdcon}}.

  They are designed to be used for MLE in
  \code{\link[evmix:fhpdcon]{fhpdcon}} but are available
  for wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fhpd]{fhpd}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:lhpdcon]{lhpdcon}}, which takes
  parameters as inputs in the same form as distribution
  functions. The negative log-likelihood is a wrapper for
  \code{\link[evmix:lhpdcon]{lhpdcon}}, designed towards
  making it useable for optimisation (e.g. parameters are
  given a vector as first input).

  The function \code{\link[evmix:lhpdcon]{lhpdcon}} carries
  out the calculations for the log-likelihood directly,
  which can be exponentiated to give actual likelihood
  using (\code{log=FALSE}).
}
\note{
  Unlike all the distribution functions for this mixture
  model, the likelihood functions only permits a scalar
  value for all the parameters. Only the data is a vector.

  The \code{\link[evmix:lhpdcon]{lhpdcon}} also has the
  usual defaults for the other parameters, but
  \code{\link[evmix:lhpdcon]{nlhpdcon}} has no defaults.

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

  Carreau, J. and Y. Bengio (2008). A hybrid Pareto model
  for asymmetric fat-tailed data: the univariate case.
  Extremes 12 (1), 53-76.
}
\seealso{
  \code{\link[evmix:lhpd]{lhpd}},
  \code{\link[evmix:lnormgpdcon]{lnormgpdcon}},
  code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}. The
  \code{\link[condmixt:condmixt-package]{condmixt}} package
  written by one of the original authors of the hybrid
  Pareto model (Carreau and Bengio, 2008) also has similar
  functions for the likelihood of the hybrid Pareto
  \code{\link[condmixt:hpareto.negloglike]{hpareto.negloglike}}
  and fitting
  \code{\link[condmixt:hpareto.negloglike]{hpareto.fit}}.

  Other hpdcon: \code{\link{dhpdcon}},
  \code{\link{fhpdcon}}, \code{\link{hpdcon}},
  \code{\link{phpdcon}}, \code{\link{qhpdcon}},
  \code{\link{rhpdcon}}
}

