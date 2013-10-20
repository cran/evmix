\name{internal}
\alias{bckdenxbeta1}
\alias{bckdenxbeta2}
\alias{bckdenxcopula}
\alias{bckdenxcutnorm}
\alias{bckdenxgamma1}
\alias{bckdenxgamma2}
\alias{bckdenxreflect}
\alias{bckdenxrenorm}
\alias{bckdenxsimple}
\alias{internal}
\alias{kdenx}
\alias{pbckdenxbeta1}
\alias{pbckdenxbeta2}
\alias{pbckdenxcopula}
\alias{pbckdenxcutnorm}
\alias{pbckdenxgamma1}
\alias{pbckdenxgamma2}
\alias{pbckdenxlog}
\alias{pbckdenxnn}
\alias{pbckdenxreflect}
\alias{pbckdenxrenorm}
\alias{pbckdenxsimple}
\alias{pkdenx}
\alias{pxb}
\title{Internal Functions}
\usage{
  kdenx(x, kerncentres, lambda, kernel = "gaussian")

  pkdenx(x, kerncentres, lambda, kernel = "gaussian")

  bckdenxsimple(x, kerncentres, lambda,
    kernel = "gaussian")

  pbckdenxsimple(x, kerncentres, lambda,
    kernel = "gaussian")

  bckdenxcutnorm(x, kerncentres, lambda,
    kernel = "gaussian")

  pbckdenxcutnorm(x, kerncentres, lambda,
    kernel = "gaussian")

  bckdenxrenorm(x, kerncentres, lambda,
    kernel = "gaussian")

  pbckdenxrenorm(x, kerncentres, lambda,
    kernel = "gaussian")

  bckdenxreflect(x, kerncentres, lambda,
    kernel = "gaussian")

  pbckdenxreflect(x, kerncentres, lambda,
    kernel = "gaussian")

  pxb(x, lambda)

  bckdenxbeta1(x, kerncentres, lambda, xmax)

  pbckdenxbeta1(x, kerncentres, lambda, xmax)

  bckdenxbeta2(x, kerncentres, lambda, xmax)

  pbckdenxbeta2(x, kerncentres, lambda, xmax)

  bckdenxgamma1(x, kerncentres, lambda)

  pbckdenxgamma1(x, kerncentres, lambda)

  bckdenxgamma2(x, kerncentres, lambda)

  pbckdenxgamma2(x, kerncentres, lambda)

  bckdenxcopula(x, kerncentres, lambda, xmax)

  pbckdenxcopula(x, kerncentres, lambda, xmax)

  pbckdenxlog(x, kerncentres, lambda, offset,
    kernel = "gaussian")

  pbckdenxnn(x, kerncentres, lambda, kernel = "gaussian",
    nn)
}
\arguments{
  \item{x}{quantiles}

  \item{kerncentres}{kernel centres (typically sample data
  vector or scalar)}

  \item{lambda}{bandwidth for kernel (as half-width of
  kernel) or \code{NULL}}

  \item{kernel}{kernel name (\code{default = "gaussian"})}

  \item{xmax}{upper bound on support (copula and beta
  kernels only) or \code{NULL}}

  \item{offset}{offset added to kernel centres (logtrans
  only) or \code{NULL}}

  \item{nn}{non-negativity correction method (simple
  boundary correction only)}
}
\description{
  Internal functions not designed to be used directly, but
  are all exported to make them visible to users.
}
\details{
  Internal functions not designed to be used directly. No
  error checking of the inputs is carried out, so user must
  be know what they are doing. They are undocumented, but
  are made visible to the user.

  Mostly, these are used in the kernel density estimation
  functions.
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}. Based on code by
  Anna MacDonald produced for MATLAB.
}
\seealso{
  \code{\link[stats:density]{density}},
  \code{\link[evmix:kden]{kden}} and
  \code{\link[evmix:bckden]{bckden}}.
}

