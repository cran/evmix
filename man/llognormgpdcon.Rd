\name{llognormgpdcon}
\alias{llognormgpdcon}
\alias{nllognormgpdcon}
\title{Log-likelihood of Log-Normal Bulk and GPD Tail Extreme Value Mixture Model with Continuity Constraint}
\usage{
  llognormgpdcon(x, lnmean = 0, lnsd = 1,
    u = qlnorm(0.9, lnmean, lnsd), xi = 0, phiu = TRUE,
    log = TRUE)

  nllognormgpdcon(pvector, x, phiu = TRUE,
    finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{lnmean}{mean on log scale}

  \item{lnsd}{standard deviation on log scale
  (non-negative)}

  \item{u}{threshold}

  \item{xi}{shape parameter}

  \item{phiu}{probability of being above threshold [0,1] or
  logical}

  \item{log}{logical, if TRUE then log density}

  \item{pvector}{vector of initial values of mixture model
  parameters (\code{lnmean}, \code{lnsd}, \code{u},
  \code{sigmau}, \code{xi}) or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}
}
\value{
  \code{\link[evmix:llognormgpdcon]{llognormgpdcon}} gives
  (log-)likelihood and
  \code{\link[evmix:llognormgpdcon]{nllognormgpdcon}} gives
  the negative log-likelihood.
}
\description{
  Log-likelihood and negative log-likelihood for the
  extreme value mixture model with log-normal for bulk
  distribution upto the threshold and conditional GPD above
  threshold and with a continuity constraint.
}
\details{
  The likelihood functions for the extreme value mixture
  model with log-normal bulk and GPD tail with continuity
  constraint, as used in the maximum likelihood fitting
  function
  \code{\link[evmix:flognormgpdcon]{flognormgpdcon}}.

  They are designed to be used for MLE in
  \code{\link[evmix:flognormgpdcon]{flognormgpdcon}} but
  are available for wider usage, e.g. constructing your own
  extreme value mixture models.

  Negative data are ignored.

  See \code{\link[evmix:flognormgpdcon]{flognormgpdcon}},
  \code{\link[evmix:flognormgpd]{flognormgpd}},
  \code{\link[evmix:fnormgpd]{fnormgpd}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:llognormgpdcon]{llognormgpdcon}}, which
  takes parameters as inputs in the same form as
  distribution functions. The negative log-likelihood is a
  wrapper for
  \code{\link[evmix:llognormgpdcon]{llognormgpdcon}},
  designed towards making it useable for optimisation (e.g.
  parameters are given a vector as first input). The tail
  fraction \code{phiu} is treated separately to the other
  parameters, to allow for all it's representations.

  Unlike the distribution functions
  \code{\link[evmix:lognormgpdcon]{lognormgpdcon}} the
  \code{phiu} must be either logical (\code{TRUE} or
  \code{FALSE}) or numerical in range \eqn{(0, 1)}. The
  default is to specify \code{phiu=TRUE} so that the tail
  fraction is specified by log-normal distribution
  \eqn{\phi_u = 1 - H(u)}, or \code{phiu=FALSE} to treat
  the tail fraction as an extra parameter estimated using
  the sample proportion. Specify a numeric \code{phiu} as
  pre-specified probability \eqn{(0, 1)}. Notice that the
  tail fraction probability cannot be 0 or 1 otherwise
  there would be no contribution from either tail or bulk
  components respectively.

  The function
  \code{\link[evmix:llognormgpdcon]{llognormgpdcon}}
  carries out the calculations for the log-likelihood
  directly, which can be exponentiated to give actual
  likelihood using (\code{log=FALSE}).
}
\note{
  Unlike all the distribution functions for this mixture
  model, the likelihood functions only permits a scalar
  value for all the parameters. Only the data is a vector.

  A default value for the tail fraction \code{phiu=TRUE} is
  given in both
  \code{\link[evmix:llognormgpdcon]{llognormgpdcon}} and
  \code{\link[evmix:llognormgpdcon]{nllognormgpdcon}}. The
  \code{\link[evmix:llognormgpdcon]{llognormgpdcon}} also
  has the usual defaults for the other parameters, but
  \code{\link[evmix:llognormgpdcon]{nllognormgpdcon}} has
  no defaults.

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
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other lognormgpdcon: \code{\link{dlognormgpdcon}},
  \code{\link{flognormgpdcon}},
  \code{\link{lognormgpdcon}},
  \code{\link{plognormgpdcon}},
  \code{\link{qlognormgpdcon}},
  \code{\link{rlognormgpdcon}}
}

