\name{internal}
\alias{beta1bckdenx}
\alias{beta1pbckdenx}
\alias{beta2bckdenx}
\alias{beta2pbckdenx}
\alias{copulabckdenx}
\alias{copulapbckdenx}
\alias{gamma1bckdenx}
\alias{gamma1pbckdenx}
\alias{gamma2bckdenx}
\alias{gamma2pbckdenx}
\alias{internal}
\alias{kdenx}
\alias{logpbckdenx}
\alias{nnbckdenx}
\alias{nnpbckdenx}
\alias{pkdenx}
\alias{pxb}
\alias{reflectbckdenx}
\alias{reflectpbckdenx}
\alias{renormbckdenx}
\alias{renormpbckdenx}
\alias{simplebckdenx}
\alias{simplepbckdenx}
\title{Internal Functions}
\usage{
  kdenx(x, kerncentres, lambda)

  pkdenx(x, kerncentres, lambda)

  simplebckdenx(x, kerncentres, lambda)

  simplepbckdenx(x, kerncentres, lambda)

  renormbckdenx(x, kerncentres, lambda)

  renormpbckdenx(x, kerncentres, lambda)

  reflectbckdenx(x, kerncentres, lambda)

  reflectpbckdenx(x, kerncentres, lambda)

  pxb(x, lambda)

  beta1bckdenx(x, kerncentres, lambda, xmax)

  beta1pbckdenx(x, kerncentres, lambda, xmax)

  beta2bckdenx(x, kerncentres, lambda, xmax)

  beta2pbckdenx(x, kerncentres, lambda, xmax)

  gamma1bckdenx(x, kerncentres, lambda)

  gamma1pbckdenx(x, kerncentres, lambda)

  gamma2bckdenx(x, kerncentres, lambda)

  gamma2pbckdenx(x, kerncentres, lambda)

  copulabckdenx(x, kerncentres, lambda, xmax)

  copulapbckdenx(x, kerncentres, lambda, xmax)

  logpbckdenx(x, kerncentres, lambda, offset)

  nnbckdenx(x, kerncentres, lambda)

  nnpbckdenx(x, kerncentres, lambda)
}
\arguments{
  \item{x}{quantile}

  \item{kerncentres}{kernel centres (typically sample
  data)}

  \item{lambda}{bandwidth for KDE}

  \item{xmax}{upper bound on support, for copula and beta
  based KDE's only}

  \item{offset}{offset added to kernel centres, log
  transform based KDE}
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
  \code{\link[evmix:kden]{kden}} and
  \code{\link[evmix:bckden]{bckden}}.
}

