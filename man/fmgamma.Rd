% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmgamma.r
\name{fmgamma}
\alias{fmgamma}
\alias{lmgamma}
\alias{nlmgamma}
\alias{nlEMmgamma}
\title{MLE Fitting of Mixture of Gammas Using EM Algorithm}
\usage{
fmgamma(x, M, pvector = NULL, std.err = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...)

lmgamma(x, mgshape, mgscale, mgweight, log = TRUE)

nlmgamma(pvector, x, M, finitelik = FALSE)

nlEMmgamma(pvector, tau, mgweight, x, M, finitelik = FALSE)
}
\arguments{
\item{x}{vector of sample data}

\item{M}{number of gamma components in mixture}

\item{pvector}{vector of initial values of GPD parameters (\code{sigmau}, \code{xi}) or \code{NULL}}

\item{std.err}{logical, should standard errors be calculated}

\item{method}{optimisation method (see \code{\link[stats:optim]{optim}})}

\item{control}{optimisation control list (see \code{\link[stats:optim]{optim}})}

\item{finitelik}{logical, should log-likelihood return finite value for invalid parameters}

\item{...}{optional inputs passed to \code{\link[stats:optim]{optim}}}

\item{mgshape}{mgamma shape (positive) as vector of length \code{M}}

\item{mgscale}{mgamma scale (positive) as vector of length \code{M}}

\item{mgweight}{mgamma weights (positive) as vector of length \code{M}}

\item{log}{logical, if \code{TRUE} then log-likelihood rather than likelihood is output}

\item{tau}{matrix of posterior probability of being in each component
(\code{nxM} where \code{n} is \code{length(x)})}
}
\value{
Log-likelihood is given by \code{\link[evmix:fmgamma]{lmgamma}} and it's
  wrapper for negative log-likelihood from \code{\link[evmix:fmgamma]{nlmgamma}}. 
  The conditional negative log-likelihood
  using the posterior probabilities is given by \code{\link[evmix:fmgamma]{nlEMmgamma}}.
  Fitting function \code{\link[evmix:fmgammagpd]{fmgammagpd}} using EM algorithm returns
  a simple list with the following elements
  
\tabular{ll}{
 \code{call}:      \tab \code{optim} call\cr
 \code{x}:         \tab data vector \code{x}\cr
 \code{init}:      \tab \code{pvector}\cr
 \code{optim}:     \tab complete \code{optim} output\cr
 \code{mle}:       \tab vector of MLE of parameters\cr
 \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
 \code{se}:        \tab vector of standard errors of MLE of parameters\cr
 \code{nllh}:      \tab minimum negative log-likelihood\cr
 \code{n}:         \tab total sample size\cr
 \code{M}:         \tab number of gamma components\cr
 \code{mgshape}:   \tab MLE of gamma shapes\cr
 \code{mgscale}:   \tab MLE of gamma scales\cr
 \code{mgweight}:  \tab MLE of gamma weights\cr
 \code{EMresults}: \tab EM results giving complete negative log-likelihood, estimated parameters
                        and conditional "maximisation step" negative log-likelihood and convergence result\cr
 \code{posterior}: \tab posterior probabilites\cr
}
}
\description{
Maximum likelihood estimation for fitting the mixture of gammas distribution
using the EM algorithm.
}
\details{
The weighted mixture of gammas distribution is fitted to the entire
dataset by maximum likelihood estimation using the EM algorithm. The estimated parameters,
variance-covariance matrix and their standard errors are automatically output.

The expectation step estimates the expected probability of being in each component
conditional on gamma component parameters. The maximisation step optimizes the
negative log-likelihood conditional on posterior probabilities of each observation
being in each component.

The optimisation of the likelihood for these mixture models can be very sensitive to
the initial parameter vector, as often there are numerous local modes. This is an
inherent feature of such models and the EM algorithm. The EM algorithm is guaranteed
to reach the maximum of the local mode. Multiple initial values should be considered
to find the global maximum. If the \code{pvector} is input as \code{NULL} then 
random component probabilities are simulated as the initial values, so multiple such runs
should be run to check the sensitivity to initial values. Alternatives to black-box
likelihood optimisers (e.g. simulated annealing), or moving to computational Bayesian
inference, are also worth considering.

The log-likelihood functions are provided for wider usage, e.g. constructing profile
likelihood functions. The parameter vector \code{pvector} must be specified in the
negative log-likelihood functions \code{\link[evmix:fmgamma]{nlmgamma}} and
\code{\link[evmix:fmgamma]{nlEMmgamma}}.

Log-likelihood calculations are carried out in \code{\link[evmix:fmgamma]{lmgamma}},
which takes parameters as inputs in the same form as the distribution functions. The
negative log-likelihood function \code{\link[evmix:fmgamma]{nlmgamma}} is a wrapper
for \code{\link[evmix:fmgamma]{lmgamma}} designed towards making it useable for optimisation,
i.e. \code{\link[evmix:fmgamma]{nlmgamma}} has complete parameter vector as first input.
Similarly, for the maximisation step negative log-likelihood
\code{\link[evmix:fmgamma]{nlEMmgamma}}, which also has the second input as the component
probability vector \code{mgweight}.

Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored.

The function \code{\link[evmix:fnormgpd]{lnormgpd}} carries out the calculations
for the log-likelihood directly, which can be exponentiated to give actual
likelihood using (\code{log=FALSE}).

The default optimisation algorithm in the "maximisation step" is "BFGS", which
requires a finite negative 
log-likelihood function evaluation \code{finitelik=TRUE}. For invalid 
parameters, a zero likelihood is replaced with \code{exp(-1e6)}. The "BFGS" 
optimisation algorithms require finite values for likelihood, so any user 
input for \code{finitelik} will be overridden and set to \code{finitelik=TRUE} 
if either of these optimisation methods is chosen.

It will display a warning for non-zero convergence result comes from 
\code{\link[stats:optim]{optim}} function call or for common indicators of lack
of convergence (e.g. any estimated parameters same as initial values).

If the hessian is of reduced rank then the variance covariance (from inverse hessian)
and standard error of parameters cannot be calculated, then by default 
\code{std.err=TRUE} and the function will stop. If you want the parameter estimates
even if the hessian is of reduced rank (e.g. in a simulation study) then
set \code{std.err=FALSE}. 

Suppose there are \eqn{M} gamma components with (scalar) shape and scale parameters and
weight for each component. Only \eqn{M-1} are to be provided in the initial parameter
vector, as the \eqn{M}th components weight is uniquely determined from the others.

For the fitting function \code{\link[evmix:fmgamma]{fmgamma}} and negative log-likelihood
functions the parameter vector \code{pvector} is a \code{3*M-1} length vector
containing all \eqn{M} gamma component shape parameters first, 
followed by the corresponding \eqn{M} gamma scale parameters,
then all the corresponding \eqn{M-1} probability weight parameters. The full parameter vector
is then \code{c(mgshape, mgscale, mgweight[1:(M-1)])}.

For the maximisation step negative log-likelihood functions the parameter vector
\code{pvector} is a \code{2*M} length vector containing all \eqn{M} gamma component
shape parameters first followed by the corresponding \eqn{M} gamma scale parameters. The
partial parameter vector is then \code{c(mgshape, mgscale)}.

For identifiability purposes the mean of each gamma component must be in ascending in order. 
If the initial parameter vector does not satisfy this constraint then an error is given. 

Non-positive data are ignored as likelihood is infinite, except for \code{gshape=1}.
}
\note{
In the fitting and profile likelihood functions, when \code{pvector=NULL} then the default initial values
are obtained under the following scheme:
\itemize{
 \item number of sample from each component is simulated from symmetric multinomial distribution;
 \item sample data is then sorted and split into groups of this size (works well when components
       have modes which are well separated);
 \item for data within each component approximate MLE's for the
       gamma shape and scale parameters are estimated.
}
The \code{\link[evmix:fmgamma]{lmgamma}}, \code{\link[evmix:fmgamma]{nlmgamma}} and
\code{\link[evmix:fmgamma]{nlEMmgamma}} have no defaults.

If the hessian is of reduced rank then the variance covariance (from inverse hessian)
and standard error of parameters cannot be calculated, then by default 
\code{std.err=TRUE} and the function will stop. If you want the parameter estimates
even if the hessian is of reduced rank (e.g. in a simulation study) then
set \code{std.err=FALSE}. 

Invalid parameter ranges will give \code{0} for likelihood, \code{log(0)=-Inf} for
log-likelihood and \code{-log(0)=Inf} for negative log-likelihood. 

Infinite and missing sample values are dropped.

Error checking of the inputs is carried out and will either stop or give warning message
as appropriate.
}
\section{Acknowledgments}{
 Thanks to Daniela Laas, University of St Gallen, Switzerland for reporting various bugs in these functions.
}

\examples{
\dontrun{
set.seed(1)
par(mfrow = c(1, 1))

x = c(rgamma(1000, shape = 1, scale = 1), rgamma(3000, shape = 6, scale = 2))
xx = seq(-1, 40, 0.01)
y = (dgamma(xx, shape = 1, scale = 1) + 3 * dgamma(xx, shape = 6, scale = 2))/4

# Fit by EM algorithm
fit = fmgamma(x, M = 2)
hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
lines(xx, y)
with(fit, lines(xx, dmgamma(xx, mgshape, mgscale, mgweight), col="red"))
}

}
\references{
\url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}

\url{http://en.wikipedia.org/wiki/Gamma_distribution}

\url{http://en.wikipedia.org/wiki/Mixture_model}

McLachlan, G.J. and Peel, D. (2000). Finite Mixture Models. Wiley.
}
\seealso{
\code{\link[stats:GammaDist]{dgamma}} and \code{\link[mixtools:gammamixEM]{gammamixEM}}
 in \code{mixtools} package

Other gammagpd: \code{\link{fgammagpdcon}},
  \code{\link{fgammagpd}}, \code{\link{fmgammagpd}},
  \code{\link{gammagpdcon}}, \code{\link{gammagpd}},
  \code{\link{mgammagpd}}

Other mgamma: \code{\link{fmgammagpdcon}},
  \code{\link{fmgammagpd}}, \code{\link{mgammagpdcon}},
  \code{\link{mgammagpd}}, \code{\link{mgamma}}

Other mgammagpd: \code{\link{fgammagpd}},
  \code{\link{fmgammagpdcon}}, \code{\link{fmgammagpd}},
  \code{\link{gammagpd}}, \code{\link{mgammagpdcon}},
  \code{\link{mgammagpd}}, \code{\link{mgamma}}

Other mgammagpdcon: \code{\link{fgammagpdcon}},
  \code{\link{fmgammagpdcon}}, \code{\link{fmgammagpd}},
  \code{\link{gammagpdcon}}, \code{\link{mgammagpdcon}},
  \code{\link{mgammagpd}}, \code{\link{mgamma}}

Other fmgamma: \code{\link{mgamma}}
}
\author{
Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}
}
\concept{fmgamma}
\concept{gammagpd}
\concept{mgamma}
\concept{mgammagpd}
\concept{mgammagpdcon}
