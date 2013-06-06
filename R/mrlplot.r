#' @export
#' 
#' @title Mean Residual Life Plot
#'
#' @description Plots the sample mean residual life (MRL) plot.
#'
#' @param data       vector of sample data
#' @param tlim       vector of (lower, upper) limits of range of threshold
#'                   to plot MRL, or \code{NULL} to use default values
#' @param nt         number of thresholds for which to evaluate MRL
#' @param alpha      logical, significance level (0, 1)
#' @param p.or.n     logical, should tail fraction (\code{FALSE}) or number of
#'                   exceedances (\code{TRUE}) be given on upper x-axis
#' @param ylim       y-axis limits or \code{NULL}
#' @param legend.loc location of legend (see \code{\link[graphics:legend]{legend}})
#' @param try.thresh vector of threshold to fit GPD using MLE and show theoretical MRL
#' @param main       title of plot
#' @param xlab       x-axis label
#' @param ylab       y-axis label
#' @param ...        further arguments to be passed to the plotting functions
#' 
#' @details Plots the sample mean residual life plot, which is also known as the mean
#' excess plot. 
#' 
#' The mean residual life above a threshold \eqn{u} is given by \code{mean(x[x > u]) - u},
#' i.e. the sample mean of the exceedances less the threshold. If the generalised Pareto
#' distribution (GPD)is an appropriate model for the excesses above \eqn{u}, then for any higher
#' thresholds \eqn{v > u} the MRL will be linear with intercept \eqn{(\sigma_u - \xi *u)/(1 - \xi)}
#' and gradient \eqn{\xi/(1 - \xi)}.
#' 
#' Symmetric central limit theorem based confidence intervals are provided for all mean
#' excesses, provided there are at least 5 exceedances. The sampling density for the MRL
#' is shown by a greyscale image, where lighter greys indicate low density.
#' 
#' A pre-chosen threshold (or more than one) can be given in \code{try.thresh}. The GPD is
#' fitted to the excesses using maximum likelihood estimation. The estimated parameters are
#' used to plot the linear function for all higher thresholds using a solid line. The threshold
#' should set as low as possible, so a dashed line is shown below the pre-chosen threshold.
#' If the MRL is similar to the dashed line then a lower threshold may be chosen.
#' 
#' If no threshold limits are provided \code{tlim = NULL} then the lowest threshold is set
#' to be just below the median data point and the maximum threshold is set to the 6th
#' largest datapoint.
#' 
#' The range of permitted thresholds is just below the minimum datapoint and the
#' second largest value. If there are less unique values of data within the threshold
#' range than the number of threshold evalations requested, then instead of a sequence
#' of thresholds they will be set to each unique datapoint, i.e. the MRL will only be evaluated
#' where there is data.
#' 
#' The missing (\code{NA} and \code{NaN}) and non-finite values are ignored.
#' 
#' The lower x-axis is the threshold and an upper axis either gives the number of 
#' exceedances (\code{p.or.n = FALSE}) or proportion of excess (\code{p.or.n = TRUE}).
#' Note that unlike the \code{gpd} related functions the missing values are ignored, so
#' do not add to the lower tail fraction. But ignoring the missing values is consistent
#' with all the other mixture model functions.
#' 
#' @return \code{\link[evmix:mrlplot]{mrlplot}} gives the mean residual life plot. It also
#' returns a matrix containing columns of the threshold, number of exceedances, mean excess,
#' standard devation of excesses and \eqn{100(1 - \alpha)\%} confidence interval. The standard
#' deviation and confidence interval are \code{NA} for less than 5 exceedances.
#' 
#' @note If the user specifies the threshold range, the thresholds above the second
#' largest are dropped. A warning message is given if any thresholds have at most 5 
#' exceedances, in which case the confidence interval is not calculated as it is
#' unreliable due to small sample. If there are less than 10 exceedances of the minimum
#' threshold then the function will stop.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Coles S.G. (2004). An Introduction to the Statistical Modelling of Extreme Values.
#' Springer-Verlag: London.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[evd:mrlplot]{mrlplot}} from 
#' \code{\link[evd:mrlplot]{evd}} library
#' 
#' @examples
#' x = rnorm(1000)
#' mrlplot(x)
#' mrlplot(x, tlim = c(0, 2.2))
#' mrlplot(x, tlim = c(0, 2), try.thresh = c(0.5, 1, 1.5))
#' mrlplot(x, tlim = c(0, 3), try.thresh = c(0.5, 1, 1.5))

mrlplot <- function(data, tlim = NULL, nt = min(100, length(data)), p.or.n = FALSE,
  alpha = 0.05, ylim = NULL, legend.loc = "bottomleft",
  try.thresh = quantile(data, 0.9, na.rm = TRUE), 
  main = "Mean Residual Life Plot", xlab = "Threshold u", ylab = "Mean Excess", ...) {

  # make sure defaults which result from function evaluations are obtained
  invisible(nt)
  invisible(try.thresh)
  
  # Check properties of inputs
  if (missing(data))
    stop("data must be a non-empty numeric vector")
    
  if (length(data) == 0 | mode(data) != "numeric") 
    stop("data must be a non-empty numeric vector of data")

  if (any(is.infinite(data)))
    stop("infinite cases must be removed")

  if (any(is.na(data)))
    warning("missing values have been removed")

  if (!is.null(tlim)) {
    if (length(tlim) != 2 | mode(tlim) != "numeric") 
      stop("threshold range tlim must be a numeric vector of length 2")

    if (tlim[2] <= tlim[1])
      stop("a range of thresholds must be specified by tlim")
  }
  
  if (length(p.or.n) != 1 | !is.logical(p.or.n))
    stop("p.or.n must be logical")
  
  if (length(nt) != 1 | mode(nt) != "numeric") 
    stop("number of thresholds must be a non-negative integer >= 2")
  
  if (abs(nt - round(nt)) > sqrt(.Machine$double.eps) | nt < 10)
    stop("number of thresholds must be a non-negative integer >= 2")

  if (!is.numeric(alpha) | length(alpha) != 1)
    stop("significance level alpha must be single numeric value")

  if (alpha <= 0 | alpha >= 1)
    stop("significance level alpha must be between (0, 1)")

  if (!is.null(ylim)) {
    if (length(ylim) != 2 | mode(ylim) != "numeric") 
      stop("ylim must be a numeric vector of length 2")

    if (ylim[2] <= ylim[1])
      stop("a range of y axis limits must be specified by ylim")
  }
  
  if (!(legend.loc %in% c("bottomright", "bottom", "bottomleft", "left",
    "topleft", "top", "topright", "right", "center")))
    stop("legend location not correct, see help(legend)")

  if (any(is.na(data)))
    warning("missing values have been removed")

  data = data[which(!is.na(data))]
  if (is.unsorted(data)) {
    data = sort(data)
  }

  if (is.null(tlim)) {
    thresholds = seq(median(data) - 2*.Machine$double.eps, data[length(data) - 6], length.out = nt)
  } else {
    thresholds = seq(tlim[1], tlim[2], length.out = nt)
  }

  n = length(data)
  data = data[data > min(thresholds)]

  # Trick to evaluate MRL at all datapoints if there are not too many
  udata = unique(data)
  if (length(udata) <= nt) {
    warning("less data than number of thresholds requested, so will use unique data as thresholds")
    thresholds = udata[-length(udata)]
  }

  # Check given thresholds
  nminu = sum(data > min(thresholds))
  if (nminu <= 10)
    stop("data must have more than 10 exceedances of lowest threshold")
  
  nmaxu = sum(data > max(thresholds))
  if (nmaxu == 0) {
    warning("thresholds above max of input data are dropped")
    thresholds = thresholds[thresholds <= max(data)]
    nmaxu = sum(data > max(thresholds))
  }
  if (nmaxu <= 5)
    warning("confidence intervals are not shown where there are less than 5 exceedances")

  nt = length(thresholds)
  if (nt < 2)
    stop("must be more than 1 threshold")

  if (!is.null(try.thresh)) {
    if (length(try.thresh) == 0 | mode(try.thresh) != "numeric") 
      stop("threshold to fit GPD to must be numeric scalar or vector")

    if (any((try.thresh < min(thresholds)) | (try.thresh >= max(thresholds))))
      stop("a range of thresholds must be specified by tlim")
  }
    
  me.calc <- function(u, x, alpha) {
    excesses = x[x > u] - u
    nxs = length(excesses)
    meanxs = mean(excesses)
    sdxs = ifelse(nxs <= 5, NA, sd(excesses))
    c(u, nxs, meanxs, sdxs, meanxs + qnorm(c(alpha/2, 1 - alpha/2)) * sdxs/sqrt(nxs))
  }
  me = t(sapply(thresholds, FUN = me.calc, x = data, alpha = 0.05))
  me = as.data.frame(me)
  names(me) = c("u", "nu", "mean.excess", "sd.excess", "cil.excess", "ciu.excess")
  
  # Assume CLT holds for mean excess, calculate density over all thresholds
  mes = range(me[, 5:6], na.rm = TRUE)
  merange = seq(mes[1] - (mes[2] - mes[1])/10, mes[2] + (mes[2] - mes[1])/10, length.out = 200)
  allmat = matrix(merange, nrow = nt, ncol = 200, byrow = TRUE)
  memat = matrix(me[, 3], nrow = nt, ncol = 200, byrow = FALSE)
  sdmat = matrix(me[, 4]/sqrt(me[, 2]), nrow = nt, ncol = 200, byrow = FALSE)
  z = (allmat - memat)/sdmat
  z[abs(z) > 3] = NA

  if (is.null(ylim)) {
    ylim = range(merange, na.rm = TRUE)
    ylim = ylim + c(-1, 1) * diff(ylim)/10
  }
  
  par(mar = c(5, 4, 7, 2) + 0.1)
  image(thresholds, merange, dnorm(z), col = gray(seq(1, 0.3, -0.01)),
    main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
  matplot(matrix(thresholds, nrow = nt, ncol = 3, byrow = FALSE), me[, c(3, 5, 6)], add = TRUE,
    type = "l", lty = c(1, 2, 2), col = "black", lwd = c(2, 1, 1), ...)
  box()

  naxis = rev(ceiling(2^pretty(log2(c(nmaxu, nminu)), 10)))
  naxis = naxis[(naxis > nmaxu) & (naxis < nminu)]
  nxaxis = c(min(thresholds), rev(data)[naxis+1], max(thresholds))
  naxis = c(nminu, naxis, nmaxu) 

  if ((nxaxis[length(nxaxis)] - nxaxis[length(nxaxis) - 1]) < diff(range(thresholds))/10) {
    nxaxis = nxaxis[-(length(nxaxis) - 1)]
    naxis = naxis[-(length(naxis) - 1)]
  }
  if ((nxaxis[2] - nxaxis[1]) < diff(range(thresholds))/20) {
    nxaxis = nxaxis[-2]
    naxis = naxis[-2]
  }

  if (p.or.n) {
    axis(side = 3, at = nxaxis, line = 0, labels = formatC(naxis/n, digits = 2, format = "g"))
    mtext("Tail Fraction phiu", side = 3, line = 2)
  } else {
    axis(side = 3, at = nxaxis, line = 0, labels = naxis)    
    mtext("Number of Excesses", side = 3, line = 2)
  }

  if (!is.null(try.thresh)) {
    ntry = length(try.thresh)
    mleparams = matrix(NA, nrow = 2, ncol = ntry)
    linecols = rep(c("blue", "green", "red"), length.out = ntry)
    for (i in 1:ntry) {
      # fit GPD above threshold using MLE
      fitresults = fgpd(data, try.thresh[i], std.err = FALSE)
      mleparams[, i] = fitresults$mle
      
      # use MLE to estimate intercept and gradient of MRL 
      mrlint = (fitresults$mle[1] - fitresults$mle[2] * try.thresh[i])/(1 - fitresults$mle[2])
      mrlgrad = fitresults$mle[2]/(1 - fitresults$mle[2])
    
      # Suppose to be linear above suitable threshold, different line type before and after
      lines(c(try.thresh[i], max(thresholds)), 
        mrlint + mrlgrad * c(try.thresh[i], max(thresholds)), lwd = 2, lty = 1, col = linecols[i])
      lines(c(min(thresholds), try.thresh[i]), 
        mrlint + mrlgrad * c(min(thresholds), try.thresh[i]), lwd = 2, lty = 2, col = linecols[i])
      abline(v = try.thresh[i], lty = 3, col = linecols[i])
    }
    legend(legend.loc, c("Sample Mean Excess", paste(100*(1 - alpha), "% CI"),
      paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
        "sigmau =", formatC(mleparams[1, 1:min(c(3, ntry))], digits = 2, format = "g"),
        "xi =", formatC(mleparams[2, 1:min(c(3, ntry))], digits = 2, format = "g"))),
      lty = c(1, 2, rep(1, min(c(3, ntry)))),
      lwd = c(2, 1, rep(1, min(c(3, ntry)))),
      col = c("black", "black", linecols), bg = "white")
  } else {
    legend(legend.loc, c("Sample Mean Excess", paste(100*(1 - alpha), "% CI")),
      lty = c(1, 2), lwd = c(2, 1), bg = "white")
  }
  
  invisible(me)
}
