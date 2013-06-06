\name{lkdengpd}
\alias{lkdengpd}
\alias{nlkdengpd}
\title{Cross-validation Log-likelihood of Kernel Density Estimator Using
Normal Kernel and GPD Tail Extreme Value Mixture Model}
\usage{
  lkdengpd(x, lambda = NULL, u = 0, sigmau = 1, xi = 0,
    phiu = TRUE, log = TRUE)

  nlkdengpd(pvector, x, phiu = TRUE, finitelik = FALSE)
}
\arguments{
  \item{x}{vector of sample data}

  \item{phiu}{logical}

  \item{pvector}{vector of initial values of mixture model
  parameters (\code{nmean}, \code{nsd}, \code{u},
  \code{sigmau}, \code{xi}) or \code{NULL}}

  \item{finitelik}{logical, should log-likelihood return
  finite value for invalid parameters}

  \item{lambda}{bandwidth for normal kernel (standard
  deviation of normal)}

  \item{u}{threshold}

  \item{sigmau}{scale parameter (non-negative)}

  \item{xi}{shape parameter}

  \item{log}{logical, if TRUE then log density}
}
\value{
  \code{\link[evmix:lkdengpd]{lkdengpd}} gives
  cross-validation (log-)likelihood and
  \code{\link[evmix:lkdengpd]{nlkdengpd}} gives the
  negative cross-validation log-likelihood.
}
\description{
  Cross-validation log-likelihood and negative
  log-likelihood for the kernel density estimator using a
  normal kernels and GPD Tail Extreme Value Mixture Model.
}
\details{
  The cross-validation likelihood functions for the kernel
  density estimator using normal kernel for the bulk below
  the threshold and GPD for upper tail. As used in the
  maximum likelihood fitting function
  \code{\link[evmix:fkdengpd]{fkdengpd}}.

  They are designed to be used for MLE in
  \code{\link[evmix:fkdengpd]{fkdengpd}} but are available
  for wider usage, e.g. constructing your own extreme value
  mixture models.

  See \code{\link[evmix:fkden]{fkden}} and
  \code{\link[evmix:fgpd]{fgpd}} for full details.

  Cross-validation likelihood is used for kernel density
  component, but standard likelihood is used for GPD
  component. The cross-validation likelihood for the KDE is
  obtained by leaving each point out in turn, evaluating
  the KDE at the point left out:
  \deqn{L(\lambda)\prod_{i=1}^{nb} \hat{f}_{-i}(x_i)} where
  \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(n-1)\lambda}
  \sum_{j=1: j\ne i}^{n} K(\frac{x_i - x_j}{\lambda})} is
  the KDE obtained when the \eqn{i}th datapoint is dropped
  out and then evaluated at that dropped datapoint at
  \eqn{x_i}. Notice that the KDE sum is indexed over all
  datapoints (\eqn{j=1, ..., n}, except datapoint \eqn{i})
  whether they are below the threshold or in the upper
  tail. But the likelihood product is evaluated only for
  those data below the threshold (\eqn{i=1, ..., n_b}). So
  the \eqn{j = n_b+1, ..., n} datapoints are extra kernel
  centres from the data in the upper tails which are used
  in the KDE but the likelihood is not evaluated there.

  Log-likelihood calculations are carried out in
  \code{\link[evmix:lkdengpd]{lkdengpd}}, which takes
  bandwidth in the same form as distribution functions. The
  negative log-likelihood is a wrapper for
  \code{\link[evmix:lkdengpd]{lkdengpd}}, designed towards
  making it useable for optimisation (e.g. parameters are
  given a vector as first input).

  The function \code{\link[evmix:lkdengpd]{lkdengpd}}
  carries out the calculations for the log-likelihood
  directly, which can be exponentiated to give actual
  likelihood using (\code{log=FALSE}).
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
  \code{\link[evmix:kdengpd]{kdengpd}},
  \code{\link[evmix:kden]{kden}},
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[stats:density]{density}}
}

