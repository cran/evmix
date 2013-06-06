\name{mrlplot}
\alias{mrlplot}
\title{Mean Residual Life Plot}
\usage{
  mrlplot(data, tlim = NULL, nt = min(100, length(data)),
    p.or.n = FALSE, alpha = 0.05, ylim = NULL,
    legend.loc = "bottomleft",
    try.thresh = quantile(data, 0.9, na.rm = TRUE),
    main = "Mean Residual Life Plot", xlab = "Threshold u",
    ylab = "Mean Excess", ...)
}
\arguments{
  \item{data}{vector of sample data}

  \item{tlim}{vector of (lower, upper) limits of range of
  threshold to plot MRL, or \code{NULL} to use default
  values}

  \item{nt}{number of thresholds for which to evaluate MRL}

  \item{alpha}{logical, significance level (0, 1)}

  \item{p.or.n}{logical, should tail fraction
  (\code{FALSE}) or number of exceedances (\code{TRUE}) be
  given on upper x-axis}

  \item{ylim}{y-axis limits or \code{NULL}}

  \item{legend.loc}{location of legend (see
  \code{\link[graphics:legend]{legend}})}

  \item{try.thresh}{vector of threshold to fit GPD using
  MLE and show theoretical MRL}

  \item{main}{title of plot}

  \item{xlab}{x-axis label}

  \item{ylab}{y-axis label}

  \item{...}{further arguments to be passed to the plotting
  functions}
}
\value{
  \code{\link[evmix:mrlplot]{mrlplot}} gives the mean
  residual life plot. It also returns a matrix containing
  columns of the threshold, number of exceedances, mean
  excess, standard devation of excesses and \eqn{100(1 -
  \alpha)\%} confidence interval. The standard deviation
  and confidence interval are \code{NA} for less than 5
  exceedances.
}
\description{
  Plots the sample mean residual life (MRL) plot.
}
\details{
  Plots the sample mean residual life plot, which is also
  known as the mean excess plot.

  The mean residual life above a threshold \eqn{u} is given
  by \code{mean(x[x > u]) - u}, i.e. the sample mean of the
  exceedances less the threshold. If the generalised Pareto
  distribution (GPD)is an appropriate model for the
  excesses above \eqn{u}, then for any higher thresholds
  \eqn{v > u} the MRL will be linear with intercept
  \eqn{(\sigma_u - \xi *u)/(1 - \xi)} and gradient
  \eqn{\xi/(1 - \xi)}.

  Symmetric central limit theorem based confidence
  intervals are provided for all mean excesses, provided
  there are at least 5 exceedances. The sampling density
  for the MRL is shown by a greyscale image, where lighter
  greys indicate low density.

  A pre-chosen threshold (or more than one) can be given in
  \code{try.thresh}. The GPD is fitted to the excesses
  using maximum likelihood estimation. The estimated
  parameters are used to plot the linear function for all
  higher thresholds using a solid line. The threshold
  should set as low as possible, so a dashed line is shown
  below the pre-chosen threshold. If the MRL is similar to
  the dashed line then a lower threshold may be chosen.

  If no threshold limits are provided \code{tlim = NULL}
  then the lowest threshold is set to be just below the
  median data point and the maximum threshold is set to the
  6th largest datapoint.

  The range of permitted thresholds is just below the
  minimum datapoint and the second largest value. If there
  are less unique values of data within the threshold range
  than the number of threshold evalations requested, then
  instead of a sequence of thresholds they will be set to
  each unique datapoint, i.e. the MRL will only be
  evaluated where there is data.

  The missing (\code{NA} and \code{NaN}) and non-finite
  values are ignored.

  The lower x-axis is the threshold and an upper axis
  either gives the number of exceedances (\code{p.or.n =
  FALSE}) or proportion of excess (\code{p.or.n = TRUE}).
  Note that unlike the \code{gpd} related functions the
  missing values are ignored, so do not add to the lower
  tail fraction. But ignoring the missing values is
  consistent with all the other mixture model functions.
}
\note{
  If the user specifies the threshold range, the thresholds
  above the second largest are dropped. A warning message
  is given if any thresholds have at most 5 exceedances, in
  which case the confidence interval is not calculated as
  it is unreliable due to small sample. If there are less
  than 10 exceedances of the minimum threshold then the
  function will stop.

  Error checking of the inputs (e.g. invalid probabilities)
  is carried out and will either stop or give warning
  message as appropriate.
}
\examples{
x = rnorm(1000)
mrlplot(x)
mrlplot(x, tlim = c(0, 2.2))
mrlplot(x, tlim = c(0, 2), try.thresh = c(0.5, 1, 1.5))
mrlplot(x, tlim = c(0, 3), try.thresh = c(0.5, 1, 1.5))
}
\author{
  Yang Hu and Carl Scarrott
  \email{carl.scarrott@canterbury.ac.nz}
}
\references{
  Scarrott, C.J. and MacDonald, A. (2012). A review of
  extreme value threshold estimation and uncertainty
  quantification. REVSTAT - Statistical Journal 10(1),
  33-59. Available from
  \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

  Coles S.G. (2004). An Introduction to the Statistical
  Modelling of Extreme Values. Springer-Verlag: London.
}
\seealso{
  \code{\link[evmix:gpd]{gpd}} and
  \code{\link[evd:mrlplot]{mrlplot}} from
  \code{\link[evd:mrlplot]{evd}} library
}

