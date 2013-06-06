\name{lkden}
\alias{lkden}
\alias{nlkden}
\title{Cross-validation Log-likelihood of Kernel Density Estimator Using normal Kernel}
\usage{
  lkden(x, lambda = NULL, extracentres = NULL, log = TRUE)

  nlkden(lambda, x, extracentres = NULL, finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{extracentres}{extra kernel centres used in KDE, but
  likelihood contribution not evaluated, or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}

  \item{lambda}{bandwidth for normal kernel (standard
  deviation of normal)}

  \item{log}{logical, if TRUE then log density}
}
\value{
  \code{\link[evmix:lkden]{lkden}} gives cross-validation
  (log-)likelihood and \code{\link[evmix:lkden]{nlkden}}
  gives the negative cross-validation log-likelihood.
}
\description{
  Cross-validation log-likelihood and negative
  log-likelihood for the kernel density estimator using a
  normal kernel by treating it as a mixture model.
}
\details{
  The cross-validation likelihood functions for the kernel
  density estimator using a normal density for kernel, as
  used in the maximum likelihood fitting function
  \code{\link[evmix:fkden]{fkden}}.

  They are designed to be used for MLE in
  \code{\link[evmix:fkden]{fkden}} but are available for
  wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fkden]{fkden}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Cross-validation likelihood is used for kernel density
  component, obtained by leaving each point out in turn and
  evaluating the KDE at the point left out:
  \deqn{L(\lambda)\prod_{i=1}^{n} \hat{f}_{-i}(x_i)} where
  \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(n-1)\lambda}
  \sum_{j=1: j\ne i}^{n} K(\frac{x_i - x_j}{\lambda})} is
  the KDE obtained when the \eqn{i}th datapoint is dropped
  out and then evaluated at that dropped datapoint at
  \eqn{x_i}.

  Normally for likelihood estimation of the bandwidth the
  kernel centres and the data where the likelihood is
  evaluated are the same. However, when using KDE for
  extreme value mixture modelling the likelihood only those
  data in the bulk of the distribution should contribute to
  the likelihood, but all the data (including those beyond
  the threshold) should contribute to the density estimate.
  The \code{extracentres} option allows the use to specify
  extra kernel centres used in estimating the density, but
  not evaluated in the likelihood. The default is to just
  use the existing data, so \code{extracentres=NULL}.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:lkden]{lkden}}, which takes bandwidth
  in the same form as distribution functions. The negative
  log-likelihood is a wrapper for
  \code{\link[evmix:lkden]{lkden}}, designed towards making
  it useable for optimisation (e.g. parameters are given a
  vector as first input).

  The function \code{\link[evmix:lkden]{lkden}} carries out
  the calculations for the log-likelihood directly, which
  can be exponentiated to give actual likelihood using
  (\code{log=FALSE}).
}
\note{
  Invalid bandwidth parameter will give \code{0} for
  likelihood, \code{log(0)=-Inf} for cross-validation
  log-likelihood and \code{-log(0)=Inf} for negative
  cross-validation log-likelihood.

  See \code{\link[evmix:gpd]{fgpd}} for explanation of
  \code{finitelik}.

  Error checking of the inputs is carried out and will
  either stop or give warning message as appropriate.
}
\section{Warning}{
  See warning in \code{\link[evmix:fkden]{fkden}}
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}
}
\references{
  \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}

  \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}

  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  Bowman, A.W. (1984). An alternative method of
  cross-validation for the smoothing of density estimates.
  Biometrika 71(2), 353-360.

  Duin, R.P.W. (1976). On the choice of smoothing
  parameters for Parzen estimators of probability density
  functions. IEEE Transactions on Computers C25(11),
  1175-1179.

  MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B.,
  Reale, M. and Russell, G. (2011). A flexible extreme
  value mixture model. Computational Statistics and Data
  Analysis 55(6), 2137-2157.
}
\seealso{
  \code{\link[stats:density]{density}}

  Other kden: \code{\link{dkden}}, \code{\link{fkden}},
  \code{\link{kden}}, \code{\link{pkden}},
  \code{\link{qkden}}, \code{\link{rkden}}
}

