\name{fdwm}
\alias{fdwm}
\title{MLE Fitting of Dynamically Weighted Mixture Model}
\usage{
  fdwm(x, pvector = NULL, std.err = TRUE, method = "BFGS",
    control = list(maxit = 10000), finitelik = TRUE, ...)
}
\arguments{
  \item{pvector}{vector of initial values of mixture model
  parameters (\code{wshape}, \code{wscale}, \code{cmu},
  \code{ctau}, \code{sigmau}, \code{xi}) or \code{NULL}}

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
  of MLE of model parameters and \code{phiu}\cr
  \code{allse}: \tab vector of standard error of all
  parameters and \code{phiu}\cr \code{n}: \tab total sample
  size\cr \code{wshape}: \tab MLE of Weibull shape\cr
  \code{wscale}: \tab MLE of Weibull scale\cr \code{mu}:
  \tab MLE of Cauchy location\cr \code{tau}: \tab MLE of
  Cauchy scale\cr \code{sigmau}: \tab MLE of GPD scale\cr
  \code{xi}: \tab MLE of GPD shape\cr \code{phiu}: \tab MLE
  of tail fraction\cr }

  The output list has some duplicate entries and repeats
  some of the inputs to both provide similar items to those
  from \code{\link[evd:fpot]{fpot}} and to make it as
  useable as possible.
}
\description{
  Maximum likelihood estimation for fitting the dynamically
  weighted mixture model
}
\details{
  The dynamically weighted mixture model is fitted to the
  entire dataset using maximum likelihood estimation. The
  estimated parameters, variance-covariance matrix and
  their standard errors are automatically output.

  Non-positive data are ignored.

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
  scalar values for each parameter and \code{phiu}. Only
  the data is a vector.

  When \code{pvector=NULL} then the initial values are
  calculated, type \code{fdwm} to see the default formulae
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
x = rweibull(1000, shape = 2)
xx = seq(-1, 4, 0.01)
y = dweibull(xx, shape = 2)

fit = fdwm(x, std.err = FALSE)
hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 4))
lines(xx, y)
lines(xx, ddwm(xx, wshape = fit$wshape, wscale = fit$wscale, cmu = fit$cmu, ctau = fit$ctau,
  sigmau = fit$sigmau, xi = fit$xi), col="red")
}
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

  Other dwm: \code{\link{ldwm}}, \code{\link{nldwm}}
}

