\name{fnormgpdcon}
\alias{fnormgpdcon}
\title{MLE Fitting of Normal Bulk and GPD Tail Extreme Value Mixture Model with Continuity Constraint}
\usage{
  fnormgpdcon(x, phiu = TRUE, pvector = NULL,
    std.err = TRUE, method = "BFGS",
    control = list(maxit = 10000), finitelik = TRUE, ...)
}
\arguments{
  \item{pvector}{vector of initial values mixture model
  parameters (\code{nmean}, \code{nsd}, \code{u},
  \code{xi}) or \code{NULL}}

  \item{x}{vector of sample data}

  \item{phiu}{logical}

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
  of MLE of model parameters, including \code{phiu} and
  \code{sigmau}\cr \code{allse}: \tab vector of standard
  error of all parameters, including \code{phiu} and
  \code{sigmau}\cr \code{n}: \tab total sample size\cr
  \code{nmean}: \tab MLE of normal mean\cr \code{nsd}: \tab
  MLE of normal standard deviation\cr \code{u}: \tab
  threshold\cr \code{sigmau}: \tab MLE of GPD scale\cr
  \code{xi}: \tab MLE of GPD shape\cr \code{phiu}: \tab MLE
  of tail fraction\cr }

  The output list has some duplicate entries and repeats
  some of the inputs to both provide similar items to those
  from \code{\link[evd:fpot]{fpot}} and to make it as
  useable as possible.
}
\description{
  Maximum likelihood estimation for fitting the extreme
  value mixture model with normal for bulk distribution
  upto the threshold and conditional GPD above threshold
  with a continuity constraint
}
\details{
  The extreme value mixture model with normal bulk and GPD
  tail with a continuity constraint is fitted to the entire
  dataset using maximum likelihood estimation. The
  estimated parameters, variance-covariance matrix and
  their standard errors are automatically output.

  The default value for \code{phiu=TRUE} so that the tail
  fraction is specified by normal distribution \eqn{\phi_u
  = 1 - H(u)}. When \code{phiu=FALSE} then the tail
  fraction is treated as an extra parameter estimated using
  the MLE which is the sample proportion above the
  threshold. In this case the standard error for
  \code{phiu} is estimated and output as \code{sephiu}.

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
  calculated, type \code{fnormgpdcon} to see the default
  formulae used. The mixture model fitting can be
  ***extremely*** sensitive to the initial values, so you
  if you get a poor fit then try some alternatives. Avoid
  setting the starting value for the shape parameter to
  \code{xi=0} as depending on the optimisation method it
  may be get stuck.

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

# Bulk model base tail fraction
fit = fnormgpdcon(x, phiu = TRUE, std.err = FALSE)
hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
lines(xx, y)
lines(xx, dnormgpdcon(xx, nmean = fit$nmean, nsd = fit$nsd, u = fit$u,
  xi = fit$xi, phiu = TRUE), col="red")
abline(v = fit$u)

# Parameterised tail fraction
fit2 = fnormgpdcon(x, phiu = FALSE, std.err = FALSE)
plot(xx, y, type = "l")
lines(xx, dnormgpdcon(xx, nmean = fit$nmean, nsd = fit$nsd, u = fit$u,
  xi = fit$xi, phiu = TRUE), col="red")
lines(xx, dnormgpdcon(xx, nmean = fit2$nmean, nsd = fit2$nsd, u = fit2$u,
  xi = fit2$xi, phiu = fit2$phiu), col="blue")
abline(v = fit$u, col = "red")
abline(v = fit2$u, col = "blue")
legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
  col=c("black", "red", "blue"), lty = 1)
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

  Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004).
  Bayesian analysis of extreme events with threshold
  estimation. Statistical Modelling. 4(3), 227-244.
}
\seealso{
  \code{\link[evmix:lgpd]{lgpd}} and
  \code{\link[evmix:gpd]{gpd}}

  Other normgpdcon: \code{\link{dnormgpdcon}},
  \code{\link{lnormgpdcon}}, \code{\link{nlnormgpdcon}},
  \code{\link{normgpdcon}}, \code{\link{pnormgpdcon}},
  \code{\link{qnormgpdcon}}, \code{\link{rnormgpdcon}}
}

