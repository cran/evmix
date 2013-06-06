\name{ldwm}
\alias{ldwm}
\alias{nldwm}
\title{Log-likelihood of dynamically weighted mixture model}
\usage{
  ldwm(x, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
    sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
    xi = 0, log = TRUE)

  nldwm(pvector, x, finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{pvector}{vector of initial values of mixture model
  parameters (\code{wshape}, \code{wscale}, \code{cmu},
  \code{ctau}, \code{sigmau}, \code{xi}) or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}

  \item{wshape}{Weibull shape (non-negative)}

  \item{wscale}{Weibull scale (non-negative)}

  \item{cmu}{Cauchy location}

  \item{ctau}{Cauchy scale}

  \item{sigmau}{scale parameter (non-negative)}

  \item{xi}{shape parameter}

  \item{log}{logical, if TRUE then log density}
}
\value{
  \code{\link[evmix:ldwm]{ldwm}} gives (log-)likelihood and
  \code{\link[evmix:ldwm]{nldwm}} gives the negative
  log-likelihood.
}
\description{
  Log-likelihood and negative log-likelihood for the
  dynamically weighted mixture model
}
\details{
  The likelihood functions for the dynamically weighted
  mixture model \code{\link[evmix:fdwm]{fdwm}}.

  Non-positive data are ignored.

  They are designed to be used for MLE in
  \code{\link[evmix:fdwm]{fdwm}} but are available for
  wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fdwm]{fdwm}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:ldwm]{ldwm}}, which takes parameters as
  inputs in the same form as distribution functions. The
  negative log-likelihood is a wrapper for
  \code{\link[evmix:ldwm]{ldwm}}, designed towards making
  it useable for optimisation (e.g. parameters are given a
  vector as first input).

  The function \code{\link[evmix:ldwm]{ldwm}} carries out
  the calculations for the log-likelihood directly, which
  can be exponentiated to give actual likelihood using
  (\code{log=FALSE}).
}
\note{
  Unlike all the distribution functions for this mixture
  model, the likelihood functions only permits a scalar
  value for all the parameters. Only the data is a vector.

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
  \url{http://en.wikipedia.org/wiki/Weibull_distribution}

  \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  Frigessi, A., O. Haug, and H. Rue (2002). A dynamic
  mixture model for unsupervised tail estimation without
  threshold selection. Extremes 5 (3), 219-235
}
\seealso{
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other dwm: \code{\link{fdwm}}
}

