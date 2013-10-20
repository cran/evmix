#' @export
#' 
#' @title Parameter Threshold Stability Plots
#'
#' @description Plots the sample mean residual life (MRL) plot.
#'
#' @inheritParams mrlplot
#' @param ylim.xi     y-axis limits for shape parameter or \code{NULL}
#' @param ylim.sigmau y-axis limits for scale parameter or \code{NULL}
#' 
#' @details The MLE of the (modified) GPD scale and shape (xi) parameters are plotted against 
#' as reange of possible threshold. Known as the threshold stability plots
#' (Coles, 2001). The modified scale parameter is \eqn{\sigma_u - u\xi}. If the
#' GPD is a suitable model for a threshold \eqn{u} then for all higher threshold
#' \eqn{v > u} it will also be suitable, with the shape and modified scale being
#' constant.
#' 
#' In practice there is sample uncertainty in the parameter estimates, which
#' must be taken into account when choosing a threshold.
#' 
#' The usual asymptotic Wald confidence intervals are shown based on the
#' observed information matrix to measure this undertainty. The sampling density
#' of the Wald normal approximation is shown by a greyscale image, where lighter
#' greys indicate low density.
#' 
#' A pre-chosen threshold (or more than one) can be given in \code{try.thresh}.
#' The GPD is fitted to the excesses using maximum likelihood estimation. The
#' estimated parameters are plot as a horizontal line which is solid above this
#' threshold where the parameter from
#' smaller tail fraction should be the same if the GPD is a good model (upto sample uncertainty).
#' The threshold should always be chosen to be as low as possible to reduce sample uncertainty.
#' Therefore, below the pre-chosen threshold, where the GPD should not be a good model, the line
#' is dashed and the parameter estimates should now deviate from the dashed line
#' (otherwise a lower threshold could be used).
# 
#' If no threshold limits are provided \code{tlim = NULL} then the lowest threshold is set
#' to be just below the median data point and the maximum threshold is set to the 11th
#' largest datapoint. This is a slightly lower order statistic compared to that used in the MRL plot 
#' \code{\link[evmix:mrlplot]{mrlplot}} function to account for the fact the maximum likelihood
#' estimation is likely to be very unreliable with 10 or fewer datapoints.
#' 
#' The range of permitted thresholds is just below the minimum datapoint and the
#' second largest value. If there are less unique values of data within the threshold
#' range than the number of threshold evalations requested, then instead of a sequence
#' of thresholds they will be set to each unique datapoint, i.e. MLE will only be applied
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
#' @return \code{\link[evmix:tcplot]{tshapeplot}} and 
#' \code{\link[evmix:tcplot]{tscaleplot}} produces the threshold stability plot for the
#' shape and scale parameter respectively. They also returns a matrix containing columns of
#' the threshold, number of exceedances, MLE shape/scale
#' and their standard devation and \eqn{100(1 - \alpha)\%} Wald confidence interval. Where the
#' observed information matrix is not obtainable the standard deviation and confidence intervals
#' are \code{NA}. For the \code{\link[evmix:tcplot]{tscaleplot}} the modified scale quantities
#' are also provided. \code{\link[evmix:tcplot]{tcplot}} produces both plots on one graph and
#' outputs a merged dataframe of results.
#' 
#' @note If the user specifies the threshold range, the thresholds above the sixth
#' largest are dropped. A warning message is given if any thresholds have at most 10
#' exceedances, in which case the maximum likelihood estimation is unreliable. If there
#' are less than 10 exceedances of the minimum threshold then the function will stop.
#' 
#' By default, no legend is included when using \code{\link[evmix:tcplot]{tcplot}} to get
#' both threshold stability plots.
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
#' @seealso \code{\link[evmix:mrlplot]{mrlplot}} and \code{\link[evd:tcplot]{tcplot}} from 
#' \code{\link[evd:mrlplot]{evd}} library
#' @aliases tcplot tshapeplot tscaleplot
#' @family tcplot
#' 
#' @examples
#' \dontrun{
#' x = rnorm(1000)
#' tcplot(x)
#' tshapeplot(x, tlim = c(0, 2))
#' tscaleplot(x, tlim = c(0, 2), try.thresh = c(0.5, 1, 1.5))
#' tcplot(x, tlim = c(0, 2), try.thresh = c(0.5, 1, 1.5))
#' }

tcplot <- function(data, tlim = NULL, nt = min(100, length(data)), p.or.n = FALSE,
  alpha = 0.05, ylim.xi = NULL, ylim.sigmau = NULL, legend.loc = NULL,
  try.thresh = quantile(data, 0.9, na.rm = TRUE), ...) {
  
  # make sure defaults which result from function evaluations are obtained
  invisible(nt)
  invisible(try.thresh)
  
  # Check properties of inputs
  if (missing(data))
    stop("data must be a non-empty numeric vector")
    
  if (length(data) == 0 | mode(data) != "numeric") 
    stop("data must be a non-empty numeric vector of data")

  if (any(is.infinite(data))) stop("infinite cases must be removed")

  if (any(is.na(data))) warning("missing values have been removed")

  data = data[which(!is.na(data))]
  if (is.unsorted(data)) {
    data = sort(data)
  }

  if (!is.null(tlim)) {
    if (length(tlim) != 2 | mode(tlim) != "numeric") 
      stop("threshold range tlim must be a numeric vector of length 2")

    if (tlim[2] <= tlim[1])
      stop("a range of thresholds must be specified by tlim")
  }
  
  if (length(p.or.n) != 1 | !is.logical(p.or.n)) stop("p.or.n must be logical")
  
  if (length(nt) != 1 | mode(nt) != "numeric") 
    stop("number of thresholds must be a non-negative integer >= 2")
  
  if (abs(nt - round(nt)) > sqrt(.Machine$double.eps) | nt < 2)
    stop("number of thresholds must be a non-negative integer >= 2")

  if (!is.numeric(alpha) | length(alpha) != 1)
    stop("significance level alpha must be single numeric value")

  if (alpha <= 0 | alpha >= 1)
    stop("significance level alpha must be between (0, 1)")

  if (!is.null(ylim.xi)) {
    if (length(ylim.xi) != 2 | mode(ylim.xi) != "numeric") 
      stop("ylim must be a numeric vector of length 2")

    if (ylim.xi[2] <= ylim.xi[1])
      stop("a range of shape y axis limits must be specified by ylim.xi")
  }

  if (!is.null(ylim.sigmau)) {
    if (length(ylim.sigmau) != 2 | mode(ylim.sigmau) != "numeric") 
      stop("ylim must be a numeric vector of length 2")

    if (ylim.sigmau[2] <= ylim.sigmau[1])
      stop("a range of scale y axis limits must be specified by ylim.sigmau")
  }
  
  if (!is.null(legend.loc)) {
    if (!(legend.loc %in% c("bottomright", "bottom", "bottomleft", "left",
      "topleft", "top", "topright", "right", "center")))
      stop("legend location not correct, see help(legend)")
  }

  if (is.null(tlim)) {
    tlim = c(median(data) - 2*.Machine$double.eps, data[length(data) - 11])
  }
      
  par(mfrow = c(2, 1))
  shaperesults = tshapeplot(data, tlim, nt, p.or.n, alpha, ylim.xi, legend.loc, try.thresh, ...)
  scaleresults = tscaleplot(data, tlim, nt, p.or.n, alpha, ylim.sigmau, legend.loc, try.thresh, ...)

  invisible(merge(shaperesults, scaleresults))  
}

#' @export
#' @aliases tcplot tshapeplot tscaleplot
#' @rdname tcplot

tshapeplot <- function(data, tlim = NULL, nt = min(100, length(data)), p.or.n = FALSE,
  alpha = 0.05, ylim = NULL, legend.loc = "bottomleft",
  try.thresh = quantile(data, 0.9, na.rm = TRUE), main = "Shape Threshold Stability Plot", 
  xlab = "Threshold u", ylab = "Shape Parameter", ...) {

  # make sure defaults which result from function evaluations are obtained
  invisible(nt)
  invisible(try.thresh)
  
  # Check properties of inputs
  if (missing(data))
    stop("data must be a non-empty numeric vector")
    
  if (length(data) == 0 | mode(data) != "numeric") 
    stop("data must be a non-empty numeric vector of data")

  if (any(is.infinite(data))) stop("infinite cases must be removed")

  if (any(is.na(data))) warning("missing values have been removed")

  data = data[which(!is.na(data))]
  if (is.unsorted(data)) {
    data = sort(data)
  }

  if (!is.null(tlim)) {
    if (length(tlim) != 2 | mode(tlim) != "numeric") 
      stop("threshold range tlim must be a numeric vector of length 2")

    if (tlim[2] <= tlim[1])
      stop("a range of thresholds must be specified by tlim")
  }
  
  if (length(p.or.n) != 1 | !is.logical(p.or.n)) stop("p.or.n must be logical")
  
  if (length(nt) != 1 | mode(nt) != "numeric") 
    stop("number of thresholds must be a non-negative integer >= 2")
  
  if (abs(nt - round(nt)) > sqrt(.Machine$double.eps) | nt < 2)
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
  
  if (!is.null(legend.loc)) {
    if (!(legend.loc %in% c("bottomright", "bottom", "bottomleft", "left",
      "topleft", "top", "topright", "right", "center")))
      stop("legend location not correct, see help(legend)")
  }
  
  if (is.null(tlim)) {
    tlim = c(median(data) - 2*.Machine$double.eps, data[length(data) - 11])
  }
  
  thresholds = seq(tlim[1], tlim[2], length.out = nt)

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
  if (nmaxu <= 5) {
    warning("thresholds above 6th largest input data are dropped")
    thresholds = thresholds[thresholds < data[length(data) - 5]]
    nmaxu = sum(data > max(thresholds))
  }
  if (nmaxu <= 10) warning("maximum likelihood estimation is unreliable with less than 10 exceedances")

  nt = length(thresholds)
  if (nt < 2)
    stop("must be more than 1 threshold")

  if (!is.null(try.thresh)) {
    if (length(try.thresh) == 0 | mode(try.thresh) != "numeric") 
      stop("threshold to fit GPD to must be numeric scalar or vector")

    if (any((try.thresh < tlim[1]) | (try.thresh >= tlim[2])))
      stop("potential thresholds must be within range specifed by tlim")
  }
    
  mle.calc <- function(x, u, alpha) {
    gpdmle = fgpd(x, u)  
    if (is.null(gpdmle$se)) gpdmle$se = rep(NA, 2)
    
    c(u, sum(x > u), gpdmle$mle, gpdmle$se,
      gpdmle$sigmau + qnorm(c(alpha/2, 1 - alpha/2)) * gpdmle$se[1],
      gpdmle$xi + qnorm(c(alpha/2, 1 - alpha/2)) * gpdmle$se[2])
  }

  mleresults = matrix(NA, nrow = nt, ncol = 10)
  mleresults[1,] = as.vector(mle.calc(data, thresholds[1], alpha))
  for (i in 2:nt) {
    mleresults[i,] = mle.calc(data, thresholds[i], alpha)
  }  
  mleresults = as.data.frame(mleresults)
  names(mleresults) = c("u", "nu", "sigmau", "xi", "se.sigmau", "se.xi", 
    "cil.sigmau", "ciu.sigmau", "cil.xi", "ciu.xi")
  
  xicis = c(mleresults$cil.xi, mleresults$ciu.xi)
  xis = range(xicis[is.finite(xicis)])
  xirange = seq(xis[1] - (xis[2] - xis[1])/10, xis[2] + (xis[2] - xis[1])/10, length.out = 200)
  allmat = matrix(xirange, nrow = nt, ncol = 200, byrow = TRUE)
  ximat = matrix(mleresults$xi, nrow = nt, ncol = 200, byrow = FALSE)
  sdmat = matrix(mleresults$se.xi, nrow = nt, ncol = 200, byrow = FALSE)
  z = (allmat - ximat)/sdmat
  z[abs(z) > 3] = NA

  if (is.null(ylim)) {
    ylim = range(xis, na.rm = TRUE)
    ylim = ylim + c(-1, 1) * diff(ylim)/10
  }

  par(mar = c(5, 4, 7, 2) + 0.1)
  image(thresholds, xirange, dnorm(z), col = gray(seq(1, 0.3, -0.01)),
    main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
  matplot(matrix(thresholds, nrow = nt, ncol = 3, byrow = FALSE), 
    mleresults[, c("xi", "cil.xi", "ciu.xi")],
    add = TRUE, type = "l", lty = c(1, 2, 2), col = "black", lwd = c(2, 1, 1), ...)
  box()

  naxis = rev(ceiling(2^pretty(log2(c(nmaxu, nminu)), 10)))
  naxis = naxis[(naxis > nmaxu) & (naxis < nminu)]
  nxaxis = c(min(thresholds), rev(data)[naxis+1], max(thresholds))
  naxis = c(nminu, naxis, nmaxu) 

  if ((nxaxis[length(nxaxis)] - nxaxis[length(nxaxis) - 1]) < diff(range(thresholds))/20) {
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
      fitresults = fgpd(data, try.thresh[i], std.err = FALSE)
      mleparams[, i] = fitresults$mle

      # Suppose to be constant after suitable threshold, different line type before and after
      lines(c(try.thresh[i], max(thresholds)), rep(fitresults$xi, 2), lwd = 2, lty = 1, col = linecols[i])
      lines(c(min(thresholds), try.thresh[i]), rep(fitresults$xi, 2), lwd = 2, lty = 2, col = linecols[i])
      abline(v = try.thresh[i], lty = 3, col = linecols[i])
    }
    if (!is.null(legend.loc)) {
      legend(legend.loc, c("MLE of Shape", paste(100*(1-alpha), "% CI"),
        paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
          "sigmau =", formatC(mleparams[1, 1:min(c(3, ntry))], digits = 2, format = "g"),
          "xi =", formatC(mleparams[2, 1:min(c(3, ntry))], digits = 2, format = "g"))),
        lty = c(1, 2, rep(1, min(c(3, ntry)))), lwd = c(2, 1, rep(1, min(c(3, ntry)))),
        col = c("black", "black", linecols), bg = "white")
    }
  } else {
    if (!is.null(legend.loc)) {
      legend(legend.loc, c("MLE of Shape", paste(100*(1-alpha), "% CI")),
        lty = c(1, 2), lwd = c(2, 1), bg = "white")
    }
  }
  
  invisible(mleresults)
}

#' @export
#' @aliases tcplot tshapeplot tscaleplot
#' @rdname tcplot

tscaleplot <- function(data, tlim = NULL, nt = min(100, length(data)), p.or.n = FALSE,
  alpha = 0.05, ylim = NULL, legend.loc = "bottomleft",
  try.thresh = quantile(data, 0.9, na.rm = TRUE), main = "Modified Scale Threshold Stability Plot", 
  xlab = "Threshold u", ylab = "Modified Scale Parameter", ...) {

  # make sure defaults which result from function evaluations are obtained
  invisible(nt)
  invisible(try.thresh)
  
  # Check properties of inputs
  if (missing(data))
    stop("data must be a non-empty numeric vector")
    
  if (length(data) == 0 | mode(data) != "numeric") 
    stop("data must be a non-empty numeric vector of data")

  if (any(is.infinite(data))) stop("infinite cases must be removed")

  if (any(is.na(data))) warning("missing values have been removed")

  data = data[which(!is.na(data))]
  if (is.unsorted(data)) {
    data = sort(data)
  }

  if (!is.null(tlim)) {
    if (length(tlim) != 2 | mode(tlim) != "numeric") 
      stop("threshold range tlim must be a numeric vector of length 2")

    if (tlim[2] <= tlim[1])
      stop("a range of thresholds must be specified by tlim")
  }
  
  if (length(p.or.n) != 1 | !is.logical(p.or.n)) stop("p.or.n must be logical")
  
  if (length(nt) != 1 | mode(nt) != "numeric") 
    stop("number of thresholds must be a non-negative integer >= 2")
  
  if (abs(nt - round(nt)) > sqrt(.Machine$double.eps) | nt < 2)
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
  
  if (!is.null(legend.loc)) {
    if (!(legend.loc %in% c("bottomright", "bottom", "bottomleft", "left",
      "topleft", "top", "topright", "right", "center")))
      stop("legend location not correct, see help(legend)")
  }

  if (is.null(tlim)) {
    tlim = c(median(data) - 2*.Machine$double.eps, data[length(data) - 11])
  }
  
  thresholds = seq(tlim[1], tlim[2], length.out = nt)

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
  if (nmaxu <= 5) {
    warning("thresholds above 6th largest input data are dropped")
    thresholds = thresholds[thresholds < data[length(data) - 5]]
    nmaxu = sum(data > max(thresholds))
  }
  if (nmaxu <= 10) warning("maximum likelihood estimation is unreliable with less than 10 exceedances")

  nt = length(thresholds)
  if (nt < 2)
    stop("must be more than 1 threshold")

  if (!is.null(try.thresh)) {
    if (length(try.thresh) == 0 | mode(try.thresh) != "numeric") 
      stop("threshold to fit GPD to must be numeric scalar or vector")

    if (any((try.thresh < tlim[1]) | (try.thresh >= tlim[2])))
      stop("potential thresholds must be within range specifed by tlim")
  }
    
  mle.calc <- function(x, u, alpha) {
    gpdmle = fgpd(x, u)  
    if (is.null(gpdmle$se)) gpdmle$se = rep(NA, 2)
    if (is.null(gpdmle$cov)) {
      gpdmle$cov12 = NA
    } else {
      gpdmle$cov12 = gpdmle$cov[1, 2]
    }
    
    c(u, sum(x > u), gpdmle$mle, gpdmle$se,
      gpdmle$sigmau + qnorm(c(alpha/2, 1 - alpha/2)) * gpdmle$se[1],
      gpdmle$xi + qnorm(c(alpha/2, 1 - alpha/2)) * gpdmle$se[2], gpdmle$cov12)
  }

  mleresults = matrix(NA, nrow = nt, ncol = 11)
  mleresults[1,] = as.vector(mle.calc(data, thresholds[1], alpha))
  for (i in 2:nt) {
    mleresults[i,] = mle.calc(data, thresholds[i], alpha)
  }  
  mleresults = as.data.frame(mleresults)
  names(mleresults) = c("u", "nu", "sigmau", "xi", "se.sigmau", "se.xi", 
    "cil.sigmau", "ciu.sigmau", "cil.xi", "ciu.xi", "cov12")
  
  mleresults$mod.sigmau = mleresults$sigmau - mleresults$xi * mleresults$u
  mleresults$mod.se.sigmau = sqrt(mleresults$se.sigmau^2 - 
      2 * mleresults$u * mleresults$cov12 + (mleresults$u * mleresults$se.xi)^2)
  mleresults$mod.cil.sigmau = mleresults$mod.sigmau + qnorm(alpha/2) * mleresults$mod.se.sigmau
  mleresults$mod.ciu.sigmau = mleresults$mod.sigmau + qnorm(1 - alpha/2) * mleresults$mod.se.sigmau
 
  sigmaucis = c(mleresults$mod.cil.sigmau, mleresults$mod.ciu.sigmau)
  sigmaus = range(sigmaucis[is.finite(sigmaucis)])
  sigmaurange = seq(sigmaus[1] - (sigmaus[2] - sigmaus[1])/10, 
                    sigmaus[2] + (sigmaus[2] - sigmaus[1])/10, length.out = 200)
  allmat = matrix(sigmaurange, nrow = nt, ncol = 200, byrow = TRUE)
  sigmaumat = matrix(mleresults$mod.sigmau, nrow = nt, ncol = 200, byrow = FALSE)
  sdmat = matrix(mleresults$mod.se.sigmau, nrow = nt, ncol = 200, byrow = FALSE)
  z = (allmat - sigmaumat)/sdmat
  z[abs(z) > 3] = NA

  if (is.null(ylim)) {
    ylim = range(sigmaus, na.rm = TRUE)
    ylim = ylim + c(-1, 1) * diff(ylim)/10
  }

  par(mar = c(5, 4, 7, 2) + 0.1)
  image(thresholds, sigmaurange, dnorm(z), col = gray(seq(1, 0.3, -0.01)),
    main = main, xlab = xlab, ylab = ylab, ylim = ylim, ...)
  matplot(matrix(thresholds, nrow = nt, ncol = 3, byrow = FALSE), 
    mleresults[, c("mod.sigmau", "mod.cil.sigmau", "mod.ciu.sigmau")],
    add = TRUE, type = "l", lty = c(1, 2, 2), col = "black", lwd = c(2, 1, 1), ...)
  box()

  naxis = rev(ceiling(2^pretty(log2(c(nmaxu, nminu)), 10)))
  naxis = naxis[(naxis > nmaxu) & (naxis < nminu)]
  nxaxis = c(min(thresholds), rev(data)[naxis+1], max(thresholds))
  naxis = c(nminu, naxis, nmaxu) 

  if ((nxaxis[length(nxaxis)] - nxaxis[length(nxaxis) - 1]) < diff(range(thresholds))/20) {
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
      fitresults = fgpd(data, try.thresh[i], std.err = FALSE)
      mleparams[1, i] = fitresults$sigmau - fitresults$xi * fitresults$u
      mleparams[2, i] = fitresults$xi

      # Suppose to be constant after suitable threshold, different line type before and after
      lines(c(try.thresh[i], max(thresholds)), rep(mleparams[1, i], 2), lwd = 2, lty = 1, col = linecols[i])
      lines(c(min(thresholds), try.thresh[i]), rep(mleparams[1, i], 2), lwd = 2, lty = 2, col = linecols[i])
      abline(v = try.thresh[i], lty = 3, col = linecols[i])
    }
    if (!is.null(legend.loc)) {
      legend(legend.loc, c("MLE of Modified Scale", paste(100*(1 - alpha), "% CI"),
        paste("u =", formatC(try.thresh[1:min(c(3, ntry))], digits = 2, format = "g"),
          "sigmau =", formatC(mleparams[1, 1:min(c(3, ntry))], digits = 2, format = "g"),
          "xi =", formatC(mleparams[2, 1:min(c(3, ntry))], digits = 2, format = "g"))),
        lty = c(1, 2, rep(1, min(c(3, ntry)))),
        lwd = c(2, 1, rep(1, min(c(3, ntry)))),
       col = c("black", "black", linecols), bg = "white")
    }
  } else {
    if (!is.null(legend.loc)) {
      legend(legend.loc, c("MLE of Modified Scale", paste(100*(1 - alpha), "% CI")),
        lty = c(1, 2), lwd = c(2, 1), bg = "white")
    }
  }
  
  invisible(mleresults)
}
