\name{lbetagpd}
\alias{lbetagpd}
\alias{nlbetagpd}
\title{Log-likelihood of beta Bulk and GPD Tail Extreme Value Mixture Model}
\usage{
  lbetagpd(x, bshape1 = 1, bshape2 = 1,
    u = qbeta(0.9, bshape1, bshape2),
    sigmau = sqrt(bshape1 * bshape2/(bshape1 + bshape2)^2/(bshape1 + bshape2 + 1)),
    xi = 0, phiu = TRUE, log = TRUE)

  nlbetagpd(pvector, x, phiu = TRUE, finitelik = FALSE)
}
\arguments{
  \item{phiu}{probability of being above threshold [0,1] or
  logical}

  \item{x}{vector of sample data}

  \item{pvector}{vector of initial values mixture model
  parameters (\code{bshape1}, \code{bshape2}, \code{u},
  \code{sigmau}, \code{xi}) or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}

  \item{bshape1}{beta shape 1 (non-negative)}

  \item{bshape2}{beta shape 2 (non-negative)}

  \item{u}{threshold over \eqn{(0,1)}}

  \item{sigmau}{scale parameter (non-negative)}

  \item{xi}{shape parameter}

  \item{log}{logical, if TRUE then log density}
}
\value{
  \code{\link[evmix:lbetagpd]{lbetagpd}} gives
  (log-)likelihood and
  \code{\link[evmix:lbetagpd]{nlbetagpd}} gives the
  negative log-likelihood.
}
\description{
  Log-likelihood and negative log-likelihood for the
  extreme value mixture model with beta for bulk
  distribution upto the threshold and conditional GPD above
  threshold.
}
\details{
  The likelihood functions for the extreme value mixture
  model with beta bulk and GPD tail, as used in the maximum
  likelihood fitting function
  \code{\link[evmix:fbetagpd]{fbetagpd}}.

  Non-positive data are ignored. Values above 1 must come
  from GPD component, as threshold \code{u<1}..

  They are designed to be used for MLE in
  \code{\link[evmix:fbetagpd]{fbetagpd}} but are available
  for wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fbetagpd]{fbetagpd}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:lbetagpd]{lbetagpd}}, which takes
  parameters as inputs in the same form as distribution
  functions. The negative log-likelihood is a wrapper for
  \code{\link[evmix:lbetagpd]{lbetagpd}}, designed towards
  making it useable for optimisation (e.g. parameters are
  given a vector as first input). The tail fraction
  \code{phiu} is treated separately to the other
  parameters, to allow for all it's representations.

  Unlike the distribution functions
  \code{\link[evmix:betagpd]{betagpd}} the \code{phiu} must
  be either logical (\code{TRUE} or \code{FALSE}) or
  numerical in range \eqn{(0, 1)}. The default is to
  specify \code{phiu=TRUE} so that the tail fraction is
  specified by beta distribution \eqn{\phi_u = 1 - H(u)},
  or \code{phiu=FALSE} to treat the tail fraction as an
  extra parameter estimated using the sample proportion.
  Specify a numeric \code{phiu} as pre-specified
  probability \eqn{(0, 1)}. Notice that the tail fraction
  probability cannot be 0 or 1 otherwise there would be no
  contribution from either tail or bulk components
  respectively.

  The function \code{\link[evmix:lbetagpd]{lbetagpd}}
  carries out the calculations for the log-likelihood
  directly, which can be exponentiated to give actual
  likelihood using (\code{log=FALSE}).
}
\note{
  Unlike all the distribution functions for this mixture
  model, the likelihood functions only permits a scalar
  value for all the parameters. Only the data is a vector.

  A default value for the tail fraction \code{phiu=TRUE} is
  given in both \code{\link[evmix:lbetagpd]{lbetagpd}} and
  \code{\link[evmix:lbetagpd]{nlnormgpd}}. The
  \code{\link[evmix:lbetagpd]{lbetagpd}} also has the usual
  defaults for the other parameters, but
  \code{\link[evmix:lbetagpd]{nlbetagpd}} has no defaults.

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
  \url{http://en.wikipedia.org/wiki/Beta_distribution}

  \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  MacDonald, A. (2012). Extreme value mixture modelling
  with medical and industrial applications. PhD thesis,
  University of Canterbury, New Zealand.
  \url{http://ir.canterbury.ac.nz/bitstream/10092/6679/1/thesis_fulltext.pdf}
}
\seealso{
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other betagpd: \code{\link{betagpd}},
  \code{\link{dbetagpd}}, \code{\link{fbetagpd}},
  \code{\link{pbetagpd}}, \code{\link{qbetagpd}},
  \code{\link{rbetagpd}}
}

