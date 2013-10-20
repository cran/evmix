\name{fbetagpdcon}
\alias{fbetagpdcon}
\alias{lbetagpdcon}
\alias{nlbetagpdcon}
\alias{nlubetagpdcon}
\alias{proflubetagpdcon}
\title{MLE Fitting of beta Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint}
\usage{
  fbetagpdcon(x, phiu = TRUE, useq = NULL, fixedu = FALSE,
    pvector = NULL, std.err = TRUE, method = "BFGS",
    control = list(maxit = 10000), finitelik = TRUE, ...)

  lbetagpdcon(x, bshape1 = 1, bshape2 = 1,
    u = qbeta(0.9, bshape1, bshape2), xi = 0, phiu = TRUE,
    log = TRUE)

  nlbetagpdcon(pvector, x, phiu = TRUE, finitelik = FALSE)

  proflubetagpdcon(u, pvector, x, phiu = TRUE,
    method = "BFGS", control = list(maxit = 10000),
    finitelik = FALSE, ...)

  nlubetagpdcon(pvector, u, x, phiu = TRUE,
    finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{phiu}{probability of being above threshold \eqn{(0,
  1)} or logical, see Details in help for
  \code{\link[evmix:fnormgpd]{fnormgpd}}}

  \item{useq}{vector of thresholds (or scalar) to be
  considered in profile likelihood or \code{NULL} for no
  profile likelihood}

  \item{fixedu}{logical, should threshold be fixed (at
  either scalar value in \code{useq}, or estimated from
  maximum of profile likelihood evaluated at sequence of
  thresholds in \code{useq})}

  \item{pvector}{vector of initial values of parameters or
  \code{NULL} for default values, see below}

  \item{std.err}{logical, should standard errors be
  calculated}

  \item{method}{optimisation method (see
  \code{\link[stats:optim]{optim}})}

  \item{control}{optimisation control list (see
  \code{\link[stats:optim]{optim}})}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}

  \item{...}{optional inputs passed to
  \code{\link[stats:optim]{optim}}}

  \item{bshape1}{scalar beta shape 1 (positive)}

  \item{bshape2}{scalar beta shape 2 (positive)}

  \item{u}{scalar threshold over \eqn{(0, 1)}}

  \item{xi}{scalar shape parameter}

  \item{log}{logical, if \code{TRUE} then log-likelihood
  rather than likelihood is output}
}
\value{
  Log-likelihood is given by
  \code{\link[evmix:fbetagpdcon]{lbetagpdcon}} and it's
  wrappers for negative log-likelihood from
  \code{\link[evmix:fbetagpdcon]{nlbetagpdcon}} and
  \code{\link[evmix:fbetagpdcon]{nlubetagpdcon}}. Profile
  likelihood for single threshold given by
  \code{\link[evmix:fbetagpdcon]{proflubetagpdcon}}.
  Fitting function
  \code{\link[evmix:fbetagpdcon]{fbetagpdcon}} returns a
  simple list with the following elements

  \tabular{ll}{ \code{call}: \tab \code{optim} call\cr
  \code{x}: \tab data vector \code{x}\cr \code{init}: \tab
  \code{pvector}\cr \code{fixedu}: \tab fixed threshold,
  logical\cr \code{useq}: \tab threshold vector for profile
  likelihood or scalar for fixed threshold\cr \code{optim}:
  \tab complete \code{optim} output\cr \code{mle}: \tab
  vector of MLE of parameters\cr \code{cov}: \tab
  variance-covariance matrix of MLE of parameters\cr
  \code{se}: \tab vector of standard errors of MLE of
  parameters\cr \code{rate}: \tab \code{phiu} to be
  consistent with \code{\link[evd:fpot]{evd}}\cr
  \code{nllh}: \tab minimum negative log-likelihood\cr
  \code{n}: \tab total sample size\cr \code{bshape1}: \tab
  MLE of beta shape1\cr \code{bshape2}: \tab MLE of beta
  shape2\cr \code{u}: \tab threshold (fixed or MLE)\cr
  \code{sigmau}: \tab MLE of GPD scale (estimated from
  other parameters)\cr \code{xi}: \tab MLE of GPD shape\cr
  \code{phiu}: \tab MLE of tail fraction (bulk model or
  parameterised approach)\cr \code{se.phiu}: \tab standard
  error of MLE of tail fraction\cr }
}
\description{
  Maximum likelihood estimation for fitting the extreme
  value mixture model with beta for bulk distribution upto
  the threshold and conditional GPD above threshold with
  continuity at threshold. With options for profile
  likelihood estimation for threshold and fixed threshold
  approach.
}
\details{
  The extreme value mixture model with beta bulk and GPD
  tail with continuity at threshold is fitted to the entire
  dataset using maximum likelihood estimation. The
  estimated parameters, variance-covariance matrix and
  their standard errors are automatically output.

  See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for
  details, type \code{help fnormgpd}. Only the different
  features are outlined below for brevity.

  The GPD \code{sigmau} parameter is now specified as
  function of other parameters, see help for
  \code{\link[evmix:betagpdcon]{dbetagpdcon}} for details,
  type \code{help betagpdcon}. Therefore, \code{sigmau}
  should not be included in the parameter vector if initial
  values are provided, making the full parameter vector
  (\code{bshape1}, \code{bshape2}, \code{u}, \code{xi}) if
  threshold is also estimated and (\code{bshape1},
  \code{bshape2}, \code{xi}) for profile likelihood or
  fixed threshold approach.

  Negative data are ignored. Values above 1 must come from
  GPD component, as threshold \code{u<1}.
}
\note{
  When \code{pvector=NULL} then the initial values are:
  \itemize{ \item method of moments estimator of beta
  parameters assuming entire population is beta; and \item
  threshold 90\% quantile (not relevant for profile
  likelihood for threshold or fixed threshold approaches);
  \item MLE of GPD shape parameter above threshold. }
}
\section{Acknowledgments}{
  See Acknowledgments in
  \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help
  fnormgpd}. Based on MATLAB code written by Anna
  MacDonald.
}
\examples{
\dontrun{
par(mfrow=c(2,1))
x = rbeta(1000, shape1 = 2, shape2 = 4)
xx = seq(-0.1, 2, 0.01)
y = dbeta(xx, shape1 = 2, shape2 = 4)

# Continuity constraint
fit = fbetagpdcon(x)
hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 2))
lines(xx, y)
with(fit, lines(xx, dbetagpdcon(xx, bshape1, bshape2, u, xi), col="red"))
abline(v = fit$u, col = "red")

# No continuity constraint
fit2 = fbetagpd(x, phiu = FALSE)
with(fit2, lines(xx, dbetagpd(xx, bshape1, bshape2, u, sigmau, xi, phiu), col="blue"))
abline(v = fit2$u, col = "blue")
legend("topright", c("True Density","No continuity constraint","With continuty constraint"),
  col=c("black", "blue", "red"), lty = 1)

# Profile likelihood for initial value of threshold and fixed threshold approach
fitu = fbetagpdcon(x, useq = seq(0.5, 0.95, length = 20))
fitfix = fbetagpdcon(x, useq = seq(0.5, 0.95, length = 20), fixedu = TRUE)

hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 2))
lines(xx, y)
with(fit, lines(xx, dbetagpdcon(xx, bshape1, bshape2, u, xi), col="red"))
abline(v = fit$u, col = "red")
with(fitu, lines(xx, dbetagpdcon(xx, bshape1, bshape2, u, xi), col="purple"))
abline(v = fitu$u, col = "purple")
with(fitfix, lines(xx, dbetagpdcon(xx, bshape1, bshape2, u, xi), col="darkgreen"))
abline(v = fitfix$u, col = "darkgreen")
legend("topright", c("True Density","Default initial value (90\% quantile)",
 "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
 col=c("black", "red", "purple", "darkgreen"), lty = 1)
}
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}
}
\references{
  \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}

  \url{http://en.wikipedia.org/wiki/Beta_distribution}

  \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  Hu, Y. (2013). Extreme value mixture modelling: An R
  package and simulation study. MSc (Hons) thesis,
  University of Canterbury, New Zealand.
  \url{http://ir.canterbury.ac.nz/simple-search?query=extreme&submit=Go}

  MacDonald, A. (2012). Extreme value mixture modelling
  with medical and industrial applications. PhD thesis,
  University of Canterbury, New Zealand.
  \url{http://ir.canterbury.ac.nz/bitstream/10092/6679/1/thesis_fulltext.pdf}
}
\seealso{
  \code{\link[stats:Beta]{dbeta}},
  \code{\link[evmix:fgpd]{fgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other betagpd betagpdcon fbetagpd fbetagpdcon normgpd
  fnormgpd: \code{\link{fbetagpd}}, \code{\link{lbetagpd}},
  \code{\link{nlbetagpd}}, \code{\link{nlubetagpd}},
  \code{\link{proflubetagpd}}
}

