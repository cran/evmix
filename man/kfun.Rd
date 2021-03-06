% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kfun.r
\name{kfun}
\alias{kfun}
\alias{klambda}
\alias{kbw}
\alias{check.kinputs}
\alias{check.kernel}
\alias{check.kbw}
\alias{ka0}
\alias{ka1}
\alias{ka2}
\title{Various subsidiary kernel function, conversion of bandwidths and evaluating certain
kernel integrals.}
\usage{
check.kinputs(x, lambda, bw, kerncentres, allownull = FALSE)

check.kernel(kernel)

check.kbw(lambda, bw, allownull = FALSE)

klambda(bw = NULL, kernel = "gaussian", lambda = NULL)

kbw(lambda = NULL, kernel = "gaussian", bw = NULL)

ka0(truncpoint, kernel = "gaussian")

ka1(truncpoint, kernel = "gaussian")

ka2(truncpoint, kernel = "gaussian")
}
\arguments{
\item{x}{location to evaluate KDE (single scalar or vector)}

\item{lambda}{bandwidth for kernel (as half-width of kernel) or \code{NULL}}

\item{bw}{bandwidth for kernel (as standard deviations of kernel) or \code{NULL}}

\item{kerncentres}{kernel centres (typically sample data vector or scalar)}

\item{allownull}{logical, where TRUE permits NULL values}

\item{kernel}{kernel name (\code{default = "gaussian"})}

\item{truncpoint}{upper endpoint as standardised location \code{x/lambda}}
}
\value{
\code{\link[evmix:kfun]{klambda}} and \code{\link[evmix:kfun]{kbw}} return the
\code{lambda} and \code{bw} bandwidths respectively.

The checking functions \code{\link[evmix:kfun]{check.kinputs}},
\code{\link[evmix:kfun]{check.klambda}} and \code{\link[evmix:kfun]{check.kernel}}
will stop on errors and return no value.

\code{\link[evmix:kfun]{ka0}}, \code{\link[evmix:kfun]{ka1}} and \code{\link[evmix:kfun]{ka2}}
return the partial moment integrals specified above.
}
\description{
Functions for checking the inputs to the kernel functions, evaluating 
integrals \eqn{\int u^l K*(u) du} for \eqn{l = 0, 1, 2} and conversion between the two bandwidth
definitions.
}
\details{
Various boundary correction methods require integral of (partial moments of)
kernel within the range of support, over the range \eqn{[-1, p]} where \eqn{p}
is the \code{truncpoint} determined by the standardised distance of location \eqn{x}
where KDE is being evaluated to the lower bound of zero, i.e. \code{truncpoint = x/lambda}.
The exception is the normal kernel which has unbounded support so the \eqn{[-5*\lambda, p]} where
\code{lambda} is the standard deviation bandwidth. There is a function for each partial moment
of degree (0, 1, 2):
\itemize{
 \item \code{ka0} - \eqn{\int_{-1}^{p} K*(z) dz}
 \item \code{ka1} - \eqn{\int_{-1}^{p} u K*(z) dz}
 \item \code{ka2} - \eqn{\int_{-1}^{p} u^2 K*(z) dz}
}
Notice that when evaluated at the upper endpoint on the support \eqn{p = 1}
(or \eqn{p = \infty} for normal) these are the zeroth, first and second moments. In the
normal distribution case the lower bound on the region of integration is \eqn{\infty} but
implemented here as \eqn{-5*\lambda}. 
These integrals are all specified in closed form, there is no need for numerical integration
(except normal which uses the \code{\link[stats:Normal]{pnorm}} function). 

See \code{\link[evmix:kernels]{kpu}} for list of kernels and discussion of bandwidth 
definitions (and their default values):
\enumerate{
 \item \code{bw} - in terms of number of standard deviations of the kernel, consistent
   with the defined values in the \code{\link[stats:density]{density}} function in
   the \code{R} base libraries
 \item \code{lambda} - in terms of half-width of kernel
}
The \code{\link[evmix:kfun]{klambda}} function converts the \code{bw} to the \code{lambda}
equivalent, and \code{\link[evmix:kfun]{kbw}} applies converse. These conversions are
kernel specific as they depend on the kernel standard deviations. If both \code{bw} and
\code{lambda} are provided then the latter is used by default. If neither are provided 
(\code{bw=NULL} and \code{lambda=NULL}) then default is \code{lambda=1}.

\code{\link[evmix:kfun]{check.kinputs}} checks all the kernel function inputs,
\code{\link[evmix:kfun]{check.klambda}} checks the pair of inputted bandwidths and
\code{\link[evmix:kfun]{check.kernel}} checks the kernel names.
}
\examples{
xx = seq(-2, 2, 0.01)
plot(xx, kdgaussian(xx), type = "l", col = "black",ylim = c(0, 1.2))
lines(xx, kduniform(xx), col = "grey")
lines(xx, kdtriangular(xx), col = "blue")
lines(xx, kdepanechnikov(xx), col = "darkgreen")
lines(xx, kdbiweight(xx), col = "red")
lines(xx, kdtriweight(xx), col = "purple")
lines(xx, kdtricube(xx), col = "orange")
lines(xx, kdparzen(xx), col = "salmon")
lines(xx, kdcosine(xx), col = "cyan")
lines(xx, kdoptcosine(xx), col = "goldenrod")
legend("topright", c("Gaussian", "uniform", "triangular", "Epanechnikov",
"biweight", "triweight", "tricube", "Parzen", "cosine", "optcosine"), lty = 1,
col = c("black", "grey", "blue", "darkgreen", "red", "purple",
  "salmon", "orange", "cyan", "goldenrod"))

}
\references{
\url{http://en.wikipedia.org/wiki/Kernel_density_estimation}

\url{http://en.wikipedia.org/wiki/Kernel_(statistics)}

Wand and Jones (1995). Kernel Smoothing. Chapman & Hall.
}
\seealso{
\code{\link[evmix:kernels]{kernels}}, \code{\link[stats:density]{density}}, 
\code{\link[evmix:kden]{kden}} and \code{\link[evmix:bckden]{bckden}}.

Other kernels: \code{\link{kernels}}
}
\author{
Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}.
}
\concept{kernels}
