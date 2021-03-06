% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmgamma.r
\name{mgamma}
\alias{mgamma}
\alias{dmgamma}
\alias{pmgamma}
\alias{qmgamma}
\alias{rmgamma}
\title{Mixture of Gammas Distribution}
\usage{
dmgamma(x, mgshape = 1, mgscale = 1, mgweight = NULL, log = FALSE)

pmgamma(q, mgshape = 1, mgscale = 1, mgweight = NULL,
  lower.tail = TRUE)

qmgamma(p, mgshape = 1, mgscale = 1, mgweight = NULL,
  lower.tail = TRUE)

rmgamma(n = 1, mgshape = 1, mgscale = 1, mgweight = NULL)
}
\arguments{
\item{x}{quantiles}

\item{mgshape}{mgamma shape (positive) as list or vector}

\item{mgscale}{mgamma scale (positive) as list or vector}

\item{mgweight}{mgamma weights (positive) as list or vector (\code{NULL} for equi-weighted)}

\item{log}{logical, if TRUE then log density}

\item{q}{quantiles}

\item{lower.tail}{logical, if FALSE then upper tail probabilities}

\item{p}{cumulative probabilities}

\item{n}{sample size (positive integer)}
}
\value{
\code{\link[evmix:mgamma]{dmgamma}} gives the density, 
\code{\link[evmix:mgamma]{pmgamma}} gives the cumulative distribution function,
\code{\link[evmix:mgamma]{qmgamma}} gives the quantile function and 
\code{\link[evmix:mgamma]{rmgamma}} gives a random sample.
}
\description{
Density, cumulative distribution function, quantile function and
  random number generation for the mixture of gammas distribution. The parameters
  are the multiple gamma shapes \code{mgshape} scales \code{mgscale} and weights \code{mgweights}.
}
\details{
Distribution functions for weighted mixture of gammas.

Suppose there are \eqn{M>=1} gamma components in the mixture model. If you 
wish to have a single (scalar) value for each parameter within each of the
\eqn{M} components then these can be input as a vector of length \eqn{M}. If
you wish to input a vector of values for each parameter within each of the
\eqn{M} components, then they are input as a list with each entry the
parameter object for each component (which can either be a scalar or
vector as usual). No matter whether they are input as a vector or list there
must be \eqn{M} elements in \code{mgshape} and \code{mgscale}, one for each
gamma mixture component. Further, any vectors in the list of parameters must
of the same length of the \code{x, q, p} or equal to the sample size \code{n}, where
relevant.

If \code{mgweight=NULL} then equal weights for each component are assumed. Otherwise, 
\code{mgweight} must be a list of the same length as \code{mgshape} and 
\code{mgscale}, filled with positive values. In the latter case, the weights are rescaled
to sum to unity.

The gamma is defined on the non-negative reals. Though behaviour at zero depends on
the shape (\eqn{\alpha}):
\itemize{
 \item \eqn{f(0+)=\infty} for \eqn{0<\alpha<1};
 \item \eqn{f(0+)=1/\beta} for \eqn{\alpha=1} (exponential);
 \item \eqn{f(0+)=0} for \eqn{\alpha>1};
}
where \eqn{\beta} is the scale parameter.
}
\note{
All inputs are vectorised except \code{log} and \code{lower.tail}, and
the gamma mixture parameters can be vectorised within the list. The main
inputs (\code{x}, \code{p} or \code{q}) and parameters must be either a
scalar or a vector. If vectors are provided they must all be of the same
length, and the function will be evaluated for each element of vector. In
the case of \code{\link[evmix:mgamma]{rmgamma}} any input vector must be of
length \code{n}. The only exception is when the parameters are single scalar
values, input as vector of length \eqn{M}.

Default values are provided for all inputs, except for the fundamentals 
\code{x}, \code{q} and \code{p}. The default sample size for 
\code{\link[evmix:mgamma]{rmgamma}} is 1.

Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
\code{p} and \code{q} are passed through as is and infinite values are set to
\code{NA}. None of these are not permitted for the parameters.

Error checking of the inputs (e.g. invalid probabilities) is carried out and
will either stop or give warning message as appropriate.
}
\section{Acknowledgments}{
 Thanks to Daniela Laas, University of St Gallen, Switzerland for reporting various bugs in these functions.
}

\examples{
\dontrun{
set.seed(1)
par(mfrow = c(2, 1))

n = 1000
x = rmgamma(n, mgshape = c(1, 6), mgscale = c(1,2), mgweight = c(1, 2))
xx = seq(-1, 40, 0.01)

hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
lines(xx, dmgamma(xx, mgshape = c(1, 6), mgscale = c(1, 2), mgweight = c(1, 2)))

# By direct simulation
n1 = rbinom(1, n, 1/3) # sample size from population 1
x = c(rgamma(n1, shape = 1, scale = 1), rgamma(n - n1, shape = 6, scale = 2))

hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
lines(xx, dmgamma(xx, mgshape = c(1, 6), mgscale = c(1, 2), mgweight = c(1, 2)))
}

}
\references{
\url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}

\url{http://en.wikipedia.org/wiki/Gamma_distribution}

\url{http://en.wikipedia.org/wiki/Mixture_model}

McLachlan, G.J. and Peel, D. (2000). Finite Mixture Models. Wiley.
}
\seealso{
\code{\link[evmix:gammagpd]{gammagpd}}, \code{\link[evmix:gpd]{gpd}}
and \code{\link[stats:GammaDist]{dgamma}}

Other mgamma: \code{\link{fmgammagpdcon}},
  \code{\link{fmgammagpd}}, \code{\link{fmgamma}},
  \code{\link{mgammagpdcon}}, \code{\link{mgammagpd}}

Other mgammagpd: \code{\link{fgammagpd}},
  \code{\link{fmgammagpdcon}}, \code{\link{fmgammagpd}},
  \code{\link{fmgamma}}, \code{\link{gammagpd}},
  \code{\link{mgammagpdcon}}, \code{\link{mgammagpd}}

Other mgammagpdcon: \code{\link{fgammagpdcon}},
  \code{\link{fmgammagpdcon}}, \code{\link{fmgammagpd}},
  \code{\link{fmgamma}}, \code{\link{gammagpdcon}},
  \code{\link{mgammagpdcon}}, \code{\link{mgammagpd}}

Other fmgamma: \code{\link{fmgamma}}
}
\author{
Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}
}
\concept{fmgamma}
\concept{mgamma}
\concept{mgammagpd}
\concept{mgammagpdcon}
