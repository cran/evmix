\name{mgammagpd}
\alias{dmgammagpd}
\alias{mgammagpd}
\alias{pmgammagpd}
\alias{qmgammagpd}
\alias{rmgammagpd}
\title{Mixture of Gammas Bulk and GPD Tail Extreme Value Mixture Model}
\usage{
  dmgammagpd(x, mgshape = list(1), mgscale = list(1),
    mgweights = NULL,
    u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
    sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0,
    phiu = TRUE, log = FALSE)

  pmgammagpd(q, mgshape = list(1), mgscale = list(1),
    mgweights = NULL,
    u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
    sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0,
    phiu = TRUE, lower.tail = TRUE)

  qmgammagpd(p, mgshape = list(1), mgscale = list(1),
    mgweights = NULL,
    u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
    sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0,
    phiu = TRUE, lower.tail = TRUE)

  rmgammagpd(n = 1, mgshape = list(1), mgscale = list(1),
    mgweights = NULL,
    u = qgamma(0.9, mgshape[[1]], 1/mgscale[[1]]),
    sigmau = sqrt(mgshape[[1]]) * mgscale[[1]], xi = 0,
    phiu = TRUE)
}
\arguments{
  \item{mgshape}{mgamma shape (non-negative) as list}

  \item{mgscale}{mgamma scale (non-negative) as list}

  \item{mgweights}{mgamma weights (positive) as list or
  \code{NULL}}

  \item{x}{quantile}

  \item{u}{threshold (non-negative)}

  \item{sigmau}{scale parameter (non-negative)}

  \item{xi}{shape parameter}

  \item{phiu}{probability of being above threshold [0,1] or
  TRUE}

  \item{log}{logical, if TRUE then log density}

  \item{q}{quantile}

  \item{lower.tail}{logical, if FALSE then upper tail
  probabilities}

  \item{p}{cumulative probability}

  \item{n}{sample size (non-negative integer)}
}
\value{
  \code{\link[evmix:mgammagpd]{dmgammagpd}} gives the
  density, \code{\link[evmix:mgammagpd]{pmgammagpd}} gives
  the cumulative distribution function,
  \code{\link[evmix:mgammagpd]{qmgammagpd}} gives the
  quantile function and
  \code{\link[evmix:mgammagpd]{rmgammagpd}} gives a random
  sample.
}
\description{
  Density, cumulative distribution function, quantile
  function and random number generation for the extreme
  value mixture model with mixture of gammas for bulk
  distribution upto the threshold and conditional GPD above
  threshold. The parameters are the gamma shape
  \code{gshape} and scale \code{gscale}, threshold \code{u}
  GPD scale \code{sigmau} and shape \code{xi} and tail
  fraction \code{phiu}.
}
\details{
  Extreme value mixture model combining mixture of gammas
  for the bulk below the threshold and GPD for upper tail.
  The parameters are input as a list, with one parameter
  object in the list for each gamma component. There must
  be the same number of components in \code{mgshape} and
  \code{mgscale}. The number of objects in the parameters
  lists determines the number of components. The parameter
  object for each gamma component can either be a scalar or
  vector, consistent with the other mixture models

  If \code{mgweights=NULL} then assumes equal weights for
  each component. Otherwise, \code{mgweights} must be a
  list of the same length as \code{mgshape} and
  \code{mgscale}, filled with positive values. In the
  latter case, the weights are rescaled to sum to unity.

  The user can pre-specify \code{phiu} permitting a
  parameterised value for the tail fraction \eqn{\phi_u}.
  Alternatively, when \code{phiu=TRUE} the tail fraction is
  estimated as the tail fraction from the gamma bulk model.

  The cumulative distribution function with tail fraction
  \eqn{\phi_u} defined by the upper tail fraction of the
  gamma bulk model (\code{phiu=TRUE}), upto the threshold
  \eqn{0 < x \le u}, given by: \deqn{F(x) = H(x)} and above
  the threshold \eqn{x > u}: \deqn{F(x) = H(u) + [1 - H(u)]
  G(x)} where \eqn{H(x)} and \eqn{G(X)} are the mixture of
  gammas and conditional GPD cumulative distribution
  functions respectively.

  The cumulative distribution function for pre-specified
  \eqn{\phi_u}, upto the threshold \eqn{0 < x \le u}, is
  given by: \deqn{F(x) = (1 - \phi_u) H(x)/H(u)} and above
  the threshold \eqn{x > u}: \deqn{F(x) = \phi_u + [1 -
  \phi_u] G(x)} Notice that these definitions are
  equivalent when \eqn{\phi_u = 1 - H(u)}.

  The gamma is defined on the non-negative reals, so the
  threshold must be non-negative.

  See \code{\link[evmix:gammagpd]{gammagpd}} for details of
  simpler parametric mixture model with single gamma for
  bulk component and GPD for upper tail.
}
\note{
  All inputs are vectorised except \code{log} and
  \code{lower.tail}, and the parameters can be vectorised
  within the list. The main inputs (\code{x}, \code{p} or
  \code{q}) and parameters must be either a scalar or a
  vector. If vectors are provided they must all be of the
  same length, and the function will be evaluated for each
  element of vector. In the case of \code{rmgammagpd} any
  input vector must be of length \code{n}.

  Default values are provided for all inputs, except for
  the fundamentals \code{x}, \code{q} and \code{p}. The
  default sample size for
  \code{\link[evmix:mgammagpd]{rmgammagpd}} is 1.

  Missing (\code{NA}) and Not-a-Number (\code{NaN}) values
  in \code{x} and \code{q} are passed through as is and
  infinite values are set to \code{NA}.

  Error checking of the inputs (e.g. invalid probabilities)
  is carried out and will either stop or give warning
  message as appropriate.
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}
}
\references{
  \url{http://en.wikipedia.org/wiki/Gamma_distribution}

  \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  do Nascimento, F.F., Gamerman, D. and Lopes, H.F. (2011).
  A semiparametric Bayesian approach to extreme value
  estimation. Statistical Computing, 22(2), 661-675.
}
\seealso{
  \code{\link[evmix:gammagpd]{gammagpd}},
  \code{\link[evmix:mgammagpd]{mgammagpd}},
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[stats:GammaDist]{dgamma}}
}

