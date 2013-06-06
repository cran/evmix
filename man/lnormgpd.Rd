\name{lnormgpd}
\alias{lnormgpd}
\alias{nlnormgpd}
\title{Log-likelihood of Normal Bulk and GPD Tail Extreme Value Mixture Model}
\usage{
  lnormgpd(x, nmean = 0, nsd = 1,
    u = qnorm(0.9, nmean, nsd), sigmau = nsd, xi = 0,
    phiu = TRUE, log = TRUE)

  nlnormgpd(pvector, x, phiu = TRUE, finitelik = FALSE)
}
\arguments{
  \item{phiu}{probability of being above threshold [0,1] or
  logical}

  \item{x}{vector of sample data}

  \item{pvector}{vector of initial values of mixture model
  parameters (\code{nmean}, \code{nsd}, \code{u},
  \code{sigmau}, \code{xi}) or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}

  \item{nmean}{normal mean}

  \item{nsd}{normal standard deviation (non-negative)}

  \item{u}{threshold}

  \item{sigmau}{scale parameter (non-negative)}

  \item{xi}{shape parameter}

  \item{log}{logical, if TRUE then log density}
}
\value{
  \code{\link[evmix:lnormgpd]{lnormgpd}} gives
  (log-)likelihood and
  \code{\link[evmix:lnormgpd]{nlnormgpd}} gives the
  negative log-likelihood.
}
\description{
  Log-likelihood and negative log-likelihood for the
  extreme value mixture model with normal for bulk
  distribution upto the threshold and conditional GPD above
  threshold.
}
\details{
  The likelihood functions for the extreme value mixture
  model with normal bulk and GPD tail, as used in the
  maximum likelihood fitting function
  \code{\link[evmix:fnormgpd]{fnormgpd}}.

  They are designed to be used for MLE in
  \code{\link[evmix:fnormgpd]{fnormgpd}} but are available
  for wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fnormgpd]{fnormgpd}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:lnormgpd]{lnormgpd}}, which takes
  parameters as inputs in the same form as distribution
  functions. The negative log-likelihood is a wrapper for
  \code{\link[evmix:lnormgpd]{lnormgpd}}, designed towards
  making it useable for optimisation (e.g. parameters are
  given a vector as first input). The tail fraction
  \code{phiu} is treated separately to the other
  parameters, to allow for all it's representations.

  Unlike the distribution functions
  \code{\link[evmix:normgpd]{normgpd}} the \code{phiu} must
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

  The function \code{\link[evmix:lnormgpd]{lnormgpd}}
  carries out the calculations for the log-likelihood
  directly, which can be exponentiated to give actual
  likelihood using (\code{log=FALSE}).
}
\note{
  Unlike all the distribution functions for this mixture
  model, the likelihood functions only permits a scalar
  value for all the parameters. Only the data is a vector.

  A default value for the tail fraction \code{phiu=TRUE} is
  given in both \code{\link[evmix:lnormgpd]{lnormgpd}} and
  \code{\link[evmix:lnormgpd]{nlnormgpd}}. The
  \code{\link[evmix:lnormgpd]{lnormgpd}} also has the usual
  defaults for the other parameters, but
  \code{\link[evmix:lnormgpd]{nlnormgpd}} has no defaults.

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

  Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004).
  Bayesian analysis of extreme events with threshold
  estimation. Statistical Modelling. 4(3), 227-244.
}
\seealso{
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other normgpd: \code{\link{dnormgpd}},
  \code{\link{fnormgpd}}, \code{\link{normgpd}},
  \code{\link{pnormgpd}}, \code{\link{qnormgpd}},
  \code{\link{rnormgpd}}
}

