\name{fhpd}
\alias{fhpd}
\title{MLE Fitting of Hybrid Pareto Extreme Value Mixture Model}
\usage{
  fhpd(x, pvector = NULL, std.err = TRUE, method = "BFGS",
    control = list(maxit = 10000), finitelik = TRUE, ...)
}
\arguments{
  \item{pvector}{vector of initial values of mixture model
  parameters (\code{nmean}, \code{nsd}, \code{xi}) or
  \code{NULL}}

  \item{x}{vector of sample data}

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
}
\value{
  Returns a simple list with the following elements

  \tabular{ll}{ \code{call}: \tab \code{optim} call\cr
  \code{x}: \tab data vector \code{x}\cr \code{init}: \tab
  \code{pvector}\cr \code{optim}: \tab complete
  \code{optim} output\cr \code{mle}: \tab vector of MLE of
  model parameters\cr \code{cov}: \tab variance-covariance
  matrix of MLE of model parameters\cr \code{se}: \tab
  vector of standard errors of MLE of model parameters\cr
  \code{rate}: \tab \code{phiu} to be consistent with
  \code{\link[evd:fpot]{evd}}\cr \code{nllh}: \tab minimum
  negative log-likelihood\cr \code{allparams}: \tab vector
  of MLE of model parameters, including \code{u},
  \code{sigmau} and \code{phiu}\cr \code{allse}: \tab
  vector of standard error of all parameters, including
  \code{u}, \code{sigmau} and \code{phiu}\cr \code{n}: \tab
  total sample size\cr \code{nmean}: \tab MLE of normal
  mean\cr \code{nsd}: \tab MLE of normal standard
  deviation\cr \code{u}: \tab threshold\cr \code{sigmau}:
  \tab MLE of GPD scale\cr \code{xi}: \tab MLE of GPD
  shape\cr }

  The output list has some duplicate entries and repeats
  some of the inputs to both provide similar items to those
  from \code{\link[evd:fpot]{fpot}} and to make it as
  useable as possible.
}
\description{
  Maximum likelihood estimation for fitting the hybrid
  Pareto extreme value mixture model
}
\details{
  The hybrid Pareto model is fitted to the entire dataset
  using maximum likelihood estimation. The estimated
  parameters, variance-covariance matrix and their standard
  errors are automatically output.

  Missing values (\code{NA} and \code{NaN}) are assumed to
  be invalid data so are ignored, which is inconsistent
  with the \code{\link[evd:fpot]{evd}} library which
  assumes the missing values are below the threshold.

  The default optimisation algorithm is "BFGS", which
  requires a finite negative log-likelihood function
  evaluation \code{finitelik=TRUE}. For invalid parameters,
  a zero likelihood is replaced with \code{exp(-1e6)}. The
  "BFGS" optimisation algorithms require finite values for
  likelihood, so any user input for \code{finitelik} will
  be overridden and set to \code{finitelik=TRUE} if either
  of these optimisation methods is chosen.

  It will display a warning for non-zero convergence result
  comes from \code{\link[stats:optim]{optim}} function
  call.

  If the hessian is of reduced rank then the variance
  covariance (from inverse hessian) and standard error of
  parameters cannot be calculated, then by default
  \code{std.err=TRUE} and the function will stop. If you
  want the parameter estimates even if the hessian is of
  reduced rank (e.g. in a simulation study) then set
  \code{std.err=FALSE}.
}
\note{
  Unlike all the distribution functions for the extreme
  value mixture models, the MLE fitting only permits single
  scalar values for each parameter. Only the data is a
  vector.

  When \code{pvector=NULL} then the initial values are
  calculated, type \code{fhpd} to see the default formulae
  used. The mixture model fitting can be ***extremely***
  sensitive to the initial values, so you if you get a poor
  fit then try some alternatives. Avoid setting the
  starting value for the shape parameter to \code{xi=0} as
  depending on the optimisation method it may be get stuck.

  The fitting function will stop if infinite sample values
  are given.

  Error checking of the inputs is carried out and will
  either stop or give warning message as appropriate.
}
\examples{
\dontrun{
par(mfrow=c(2,1))
x = rnorm(1000)
xx = seq(-4, 4, 0.01)
y = dnorm(xx)

# Hybrid Pareto provides reasonable fit for asymmetric heavy tailed distribution
# but not for cases such as the normal distribution
fit = fhpd(x, std.err = FALSE)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
lines(xx, y)
lines(xx, dhpd(xx, nmean = fit$nmean, nsd = fit$nsd,
  xi = fit$xi), col="red")
abline(v = fit$u)

# Notice that if tail fraction is included a better fit is obtained
fit2 = fnormgpdcon(x, std.err = FALSE)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
lines(xx, y)
lines(xx, dnormgpdcon(xx, nmean = fit2$nmean, nsd = fit2$nsd, u = fit2$u,
  xi = fit2$xi), col="blue")
abline(v = fit2$u)
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

  Carreau, J. and Y. Bengio (2008). A hybrid Pareto model
  for asymmetric fat-tailed data: the univariate case.
  Extremes 12 (1), 53-76.
}
\seealso{
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}} The
  \code{\link[condmixt:condmixt-package]{condmixt}} package
  written by one of the original authors of the hybrid
  Pareto model (Carreau and Bengio, 2008) also has similar
  functions for the likelihood of the hybrid Pareto
  \code{\link[condmixt:hpareto.negloglike]{hpareto.negloglike}}
  and fitting
  \code{\link[condmixt:hpareto.negloglike]{hpareto.fit}}.

  Other hpd: \code{\link{dhpd}}, \code{\link{hpd}},
  \code{\link{lhpd}}, \code{\link{nlhpd}},
  \code{\link{phpd}}, \code{\link{qhpd}},
  \code{\link{rhpd}}
}

