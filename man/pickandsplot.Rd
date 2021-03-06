% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pickandsplot.r
\name{pickandsplot}
\alias{pickandsplot}
\title{Pickands Plot}
\usage{
pickandsplot(data, orderlim = NULL, tlim = NULL, y.alpha = FALSE,
  alpha = 0.05, ylim = NULL, legend.loc = "topright",
  try.thresh = quantile(data, 0.9, na.rm = TRUE),
  main = "Pickand's Plot", xlab = "order", ylab = ifelse(y.alpha,
  " tail index - alpha", "shape  - xi"), ...)
}
\arguments{
\item{data}{vector of sample data}

\item{orderlim}{vector of (lower, upper) limits of order statistics
to plot estimator, or \code{NULL} to use default values}

\item{tlim}{vector of (lower, upper) limits of range of threshold
to plot estimator, or \code{NULL} to use default values}

\item{y.alpha}{logical, should shape xi (\code{FALSE}) or tail index alpha (\code{TRUE}) be given on y-axis}

\item{alpha}{significance level over range (0, 1), or \code{NULL} for no CI}

\item{ylim}{y-axis limits or \code{NULL}}

\item{legend.loc}{location of legend (see \code{\link[graphics:legend]{legend}}) or \code{NULL} for no legend}

\item{try.thresh}{vector of thresholds to consider}

\item{main}{title of plot}

\item{xlab}{x-axis label}

\item{ylab}{y-axis label}

\item{...}{further arguments to be passed to the plotting functions}
}
\value{
\code{\link[evmix:pickandsplot]{pickandsplot}} gives Pickand's plot. It also 
  returns a dataframe containing columns of the order statistics, order, Pickand's
  estimator, it's standard devation and \eqn{100(1 - \alpha)\%} confidence
  interval (when requested).
}
\description{
Produces the Pickand's plot.
}
\details{
Produces the Pickand's plot including confidence intervals.
  
  For an ordered iid sequence \eqn{X_{(1)}\ge X_{(2)}\ge\cdots\ge X_{(n)}} 
  the Pickand's estimator of the reciprocal of the shape parameter \eqn{\xi} 
  at the \eqn{k}th order statistic is given by 
  \deqn{\hat{\xi}_{k,n}=\frac{1}{\log(2)} \log\left(\frac{X_{(k)}-X_{(2k)}}{X_{(2k)}-X_{(4k)}}\right).}
  Unlike the Hill estimator it does not assume positive data, is valid for any \eqn{\xi} and
  is location and scale invariant.
  The Pickands estimator is defined on orders \eqn{k=1, \ldots, \lfloor n/4\rfloor}. 
  
  Once a sufficiently low order statistic is reached the Pickand's estimator will
  be constant, upto sample uncertainty, for regularly varying tails. Pickand's
  plot is a plot of \deqn{\hat{\xi}_{k,n}} against the \eqn{k}. Symmetric asymptotic
  normal confidence intervals assuming Pareto tails are provided.
  
  The Pickand's estimator is for the GPD shape \eqn{\xi}, or the reciprocal of the
  tail index \eqn{\alpha=1/\xi}. The shape is plotted by default using
  \code{y.alpha=FALSE} and the tail index is plotted when \code{y.alpha=TRUE}.
  
  A pre-chosen threshold (or more than one) can be given in
  \code{try.thresh}. The estimated parameter (\eqn{\xi} or \eqn{\alpha}) at
  each threshold are plot by a horizontal solid line for all higher thresholds. 
  The threshold should be set as low as possible, so a dashed line is shown
  below the pre-chosen threshold. If Pickand's estimator is similar to the
  dashed line then a lower threshold may be chosen.
  
  If no order statistic (or threshold) limits are provided 
  \code{orderlim = tlim = NULL} then the lowest order statistic is set to \eqn{X_{(1)}} and
  highest possible value \eqn{X_{\lfloor n/4\rfloor}}. However, Pickand's estimator is always
  output for all \eqn{k=1, \ldots, \lfloor n/4\rfloor}.
  
  The missing (\code{NA} and \code{NaN}) and non-finite values are ignored.
  
  The lower x-axis is the order \eqn{k}. The upper axis is for the corresponding threshold.
}
\note{
Asymptotic Wald type CI's are estimated for non-\code{NULL} signficance level \code{alpha}
for the shape parameter, assuming exactly GPD tails. When plotting on the tail index scale,
then a simple reciprocal transform of the CI is applied which may well be sub-optimal.

Error checking of the inputs (e.g. invalid probabilities) is carried out and
will either stop or give warning message as appropriate.
}
\section{Acknowledgments}{
 Thanks to Younes Mouatasim, Risk Dynamics, Brussels for reporting various bugs in these functions.
}

\examples{
\dontrun{
par(mfrow = c(2, 1))

# Reproduce graphs from Figure 4.7 of Resnick (2007)
data(danish, package="evir")

# Pickand's plot
pickandsplot(danish, orderlim=c(1, 150), ylim=c(-0.1, 2.2),
 try.thresh=c(), alpha=NULL, legend.loc=NULL)
 
# Using default settings
pickandsplot(danish)
}
}
\references{
Pickands III, J.. (1975). Statistical inference using extreme order statistics. Annal of Statistics 3(1), 119-131.

Dekkers A. and de Haan, S. (1989). On the estimation of the extreme-value index and large quantile estimation.
Annals of Statistics 17(4), 1795-1832.

Resnick, S. (2007). Heavy-Tail Phenomena - Probabilistic and Statistical Modeling. Springer.
}
\seealso{
\code{\link[smoothtail:pickands]{pickands}}
}
\author{
Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}
}
