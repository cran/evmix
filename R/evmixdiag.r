#' @export
#' 
#' @title Diagnostic Plots for Extreme Value Mixture Models
#'
#' @description The classic four diagnostic plots for evaluating extreme value mixture models: 
#'   1) return level plot, 2) Q-Q plot, 3) P-P plot and 4) density plot. Each plot is available
#'   individually or as the usual 2x2 collection.
#'
#' @param modelfit   fitted extreme value mixture model object
#' @param upperfocus logical, should plot focus on upper tail?
#' @param ci         logical, should Monte Carlo based CI's be plotted
#' @param alpha      logical, significance level (0, 1)
#' @param N          number of Monte Carlo simulation for CI (N>=10)
#' @param legend     logical, should legend be included
#' @param ...        further arguments to be passed to the plotting functions
#' 
#' @details Model diagnostics are available for all the fitted extreme mixture models in the 
#' \code{\link[evmix:evmix-package]{evmix}} package. These \code{modelfit} is output by all the fitting 
#' functions, e.g. \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:fnormgpd]{fnormgpd}}.
#' 
#' Consistent with \code{\link[evd:plot.uvevd]{plot}} function in the 
#' \code{\link[evd:gpd]{evd}} library the \code{\link[stats:ppoints]{ppoints}} to 
#' estimate the empirical cumulative probabilities. The current default behaviour of this
#' function is to use \deqn{(i-0.5)/n} as the estimate for the \eqn{i}th order statistic of
#' the given sample of size \eqn{n}.
#' 
#' The return level plot quantile (\eqn{x_p} where \eqn{P(X \le x_p)=p} on the
#' \eqn{y}-axis and the (approximate) return period \eqn{1/p} is shown on the
#' \eqn{x}-axis. It is approximate as the tranformation 
#' \eqn{-log(-log(p)) \approx 1/p} is used for the \eqn{x}-axis, which is common in
#' extreme value application as Type I (\eqn{\xi=0}) upper tail behaviour will be
#' linear on this scale. The approximation is better for smaller upper tail probability.
#' 
#' The crosses are the empirical quantiles/return levels (i.e. the ordered sample data)
#' against their corresponding transfomred empirical return period (from 
#' \code{\link[stats:ppoints]{ppoints}}). The solid line is the theoretical return level
#' (quantile) function using the estimated model parameters. The estimated threshold 
#' \code{u} and tail fraction \code{phiu} are shown. For the two tailed models both
#' thresholds \code{ul} and \code{ur} and corresponding tail fractions 
#' \code{phiul} and \code{phiur} are shown. The approximate pointwise confidence intervals
#' for the quantiles are obtained by Monte Carlo simulation using the estimated parameters.
#' Notice that these intervals ignore the parameter estimation uncertainty.
#' 
#' The Q-Q and P-P plots have the empirical values on the \eqn{y}-axis and theoretical values
#' from the fitted model on the \eqn{x}-axis.
#' 
#' The density plot provides a histogram of the sample data overlaid with the fitted density
#' and a standard kernel density estimate using the \code{\link[stats:density]{density}}
#' function. The default settings for the \code{\link[stats:density]{density}} function are used.
#' Note that for distributions with bounded support (e.g. GPD) with high density near the
#' boundary standard kernel density estimators exhibit a negative bias due to leakage past
#' the boundary. So in this case they should not be taken too seriously.
#' 
#' For the kernel density estimates (i.e. \code{kden} and \code{bckden}) there is no threshold, 
#' so no upper tail focus is carried out.
#' 
#' See \code{\link[evd:plot.uvevd]{plot.uvevd}} for more detailed explanations of these
#' types of plots.
#' 
#' @return \code{\link[evmix:evmix.diag]{rlplot}} gives the return level plot, 
#' \code{\link[evmix:evmix.diag]{qplot}} gives the Q-Q plot,
#' \code{\link[evmix:evmix.diag]{pplot}} gives the P-P plot,
#' \code{\link[evmix:evmix.diag]{densplot}} gives density plot and
#' \code{\link[evmix:evmix.diag]{evmix.diag}} gives the collection of all 4.
#' 
#' @note For all mixture models the missing values are removed by the fitting functions 
#' (e.g. \code{\link[evmix:fnormgpd]{fnormgpd}} and \code{\link[evmix:fgng]{fgng}}).
#' However, these are retained in the GPD fitting \code{\link[evmix:fgpd]{fgpd}}, as they 
#' are interpreted as values below the threshold.
#' 
#' By default all the plots focus in on the upper tail, but they can be used to display 
#' the fit over the entire range of support.
#' 
#' You cannot pass \code{xlim} or \code{ylim} to the plotting functions via \code{...}
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' Based on model diagnostic functions in the \code{\link[evd:gpd]{evd}} package.
#' 
#' \url{http://en.wikipedia.org/wiki/Q-Q_plot}
#' 
#' \url{http://en.wikipedia.org/wiki/P-P_plot}
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
#' @seealso \code{\link[stats:ppoints]{ppoints}} and \code{\link[evd:plot.uvevd]{plot.uvevd}}
#' @aliases  evmix.diag rlplot qplot pplot densplot
#' @family   evmix.diag
#' 
#' @examples
#' \dontrun{
#' x = sort(rnorm(1000))
#' fit = fnormgpd(x)
#' evmix.diag(fit)
#' 
#' # repeat without focussing on upper tail
#' par(mfrow=c(2,2))
#' rlplot(fit, upperfocus = FALSE)
#' qplot(fit, upperfocus = FALSE)
#' pplot(fit, upperfocus = FALSE)
#' densplot(fit, upperfocus = FALSE)
#' }

evmix.diag <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000,
  legend = FALSE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "normgpdcon", "gng", "gngcon", 
    "weibullgpd",  "weibullgpdcon", "gammagpd", "gammagpdcon", "mgammagpd",
    "lognormgpd", "lognormgpdcon", "betagpd", "hpd", "hpdcon", "dwm",
    "kden", "kdengpd", "kdengpdcon", "gkg", "bckden", "bckdengpd")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
  if (!is.logical(legend))
    stop("legend must be logical")
  
  if (!is.logical(ci))
    stop("ci must be logical")
  
  if (!is.numeric(alpha) | length(alpha) != 1)
    stop("significance level alpha must be single numeric value")

  if (alpha <= 0 | alpha >= 1)
    stop("significance level alpha must be between (0, 1)")

  if (length(N) != 1 | mode(N) != "numeric") 
    stop("number of simulations must be a non-negative integer")
  
  if (abs(N - round(N)) > sqrt(.Machine$double.eps) | N < 10)
    stop("simulation size must non-negative integer >= 10")

  if (is.unsorted(modelfit$x, na.rm = TRUE))
    modelfit$x=sort(modelfit$x)
  
  par(mfrow = c(2, 2))
  rlplot(modelfit, upperfocus = upperfocus, ci = ci, alpha = alpha, N = N, legend = legend, ...)
  qplot(modelfit, upperfocus = upperfocus, ci = ci, alpha = alpha, N = N, legend = legend, ...)
  pplot(modelfit, upperfocus = upperfocus, ci = ci, alpha = alpha, N = N, legend = legend, ...)
  densplot(modelfit, upperfocus = upperfocus, legend = legend, ...)
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

rlplot <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000,
  legend = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "normgpdcon", "gng", "gngcon", 
    "weibullgpd",  "weibullgpdcon", "gammagpd", "gammagpdcon", "mgammagpd",
    "lognormgpd", "lognormgpdcon", "betagpd", "hpd", "hpdcon", "dwm",
    "kden", "kdengpd", "kdengpdcon", "gkg", "bckden", "bckdengpd")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
  if (!is.logical(legend))
    stop("legend must be logical")
  
  if (!is.logical(ci))
    stop("ci must be logical")
  
  if (!is.numeric(alpha) | length(alpha) != 1)
    stop("significance level alpha must be single numeric value")

  if (alpha <= 0 | alpha >= 1)
    stop("significance level alpha must be between (0, 1)")

  if (length(N) != 1 | mode(N) != "numeric") 
    stop("number of simulations must be a non-negative integer")
  
  if (abs(N - round(N)) > sqrt(.Machine$double.eps) | N < 10)
    stop("simulation size must non-negative integer >= 10")

  if (is.unsorted(modelfit$x, na.rm = TRUE))
    modelfit$x=sort(modelfit$x)
  
  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  if (modelname == "gpd") {
    upperfocus = TRUE  
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf 
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }
  
  n = length(modelfit$x)
  emp.prob = ppoints(n)
  trans.emp.prob = -1/log(emp.prob)

  min.emp.power = -log10(1 - 1/n/100)
  max.emp.power = ceiling(log10(max(trans.emp.prob))) + 1

  the.prob = 1 - 10^(-seq(min.emp.power, max.emp.power, 0.05))
  # Need both tail to be in fine detail
  the.prob = sort(c(the.prob, 1 - the.prob)) 
  trans.the.prob = -1/log(the.prob)

  the.quant = switch(modelname,
    gpd = qgpd(the.prob, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = qnormgpd(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = qnormgpdcon(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gng = qgng(the.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = qgngcon(the.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    weibullgpd = qweibullgpd(the.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = qweibullgpdcon(the.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gammagpd = qgammagpd(the.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = qgammagpdcon(the.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    mgammagpd = qmgammagpd(the.prob, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = qlognormgpd(the.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = qlognormgpd(the.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = qbetagpd(the.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    hpd = qhpd(the.prob, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = qhpdcon(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    dwm = qdwm(the.prob, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    kden = qkden(the.prob, modelfit$kerncentres, modelfit$lambda),
    kdengpd = qkdengpd(the.prob, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    kdengpdcon = qkdengpdcon(the.prob, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = qgkg(the.prob, modelfit$kerncentres, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    bckden = qbckden(the.prob, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = qbckdengpd(the.prob, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))

  # If focus on upper tail then must display at least upper 10%
  # even if tail fraction is smaller, otherwise may not look nice
  if (upperfocus) { 
    xlims = c(-1/log(1 - ifelse(phiu < 0.1, 0.1, phiu)), 10^max.emp.power)
    ylims = c(min(quantile(modelfit$x, 0.9, na.rm = TRUE), u), 
      max(c(modelfit$x, the.quant), na.rm = TRUE))
  } else {
    xlims = c(min(trans.emp.prob), 10^max.emp.power)
    ylims = range(c(modelfit$x, the.quant), na.rm = TRUE)
  }

  if (ci) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.quantiles <- function(i, modelfit, modelname) {
      simdata = switch(modelname,
        gpd = rgpd(n, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = rnormgpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gng = rgng(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = rgngcon(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
        weibullgpd = rweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = rweibullgpdcon(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gammagpd = rgammagpd(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = rgammagpdcon(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        mgammagpd = rmgammagpd(n, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = rbetagpd(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        hpd = rhpd(n, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = rhpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        dwm = rdwm(n, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        kden = rkden(n, modelfit$kerncentres, modelfit$lambda),
        kdengpd = rkdengpd(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        kdengpdcon = rkdengpdcon(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = rgkg(n, modelfit$kerncentres, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        bckden = rbckden(n, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = rbckdengpd(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))  
      sort(simdata, na.last = FALSE)
    }
    sim.q = sapply(1:N, FUN = simulate.new.quantiles, modelfit = modelfit, modelname = modelname)
    ci.q = apply(sim.q, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    # simulation for gpd allow for uncertainty in phiu, so need to ignore CI below threshold
    if (modelname == "gpd") {
      ci.q[1, modelfit$x < u] = NA
      ci.q[2, modelfit$x < u] = NA
    }
    ylims[2] = max(c(ci.q, ylims[2]), na.rm = TRUE)
  }

  plot(trans.emp.prob, modelfit$x, pch = 'x', cex = 0.8, log = "x",
    xlab = "Return Period", ylab = "Return Level",
    xlim = xlims, ylim = ylims, axes = FALSE, ...)
  lines(trans.the.prob, the.quant, col="red", lwd = 2)

  xpower = seq(0, max.emp.power)
  axis(1, at = 10^xpower, labels = formatC(10^xpower, format = "d"))
  axis(2)
  box()

  if (twotail) {
    abline(h = c(modelfit$ul, modelfit$ur), lty = 3)
  } else {
    abline(h = u, lty = 3)
  }

  if (!nothresh)
    abline(v = -1/log(1 - phiu), lty = 3)
  
  if (ci) {
    lines(trans.emp.prob, ci.q[1, ], lty=2)
    lines(trans.emp.prob, ci.q[2, ], lty=2)
    if (legend) {
      legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
        paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
        ifelse(twotail & !upperfocus, 
          "Tail Fractions and Thresholds","Upper Tail Fraction and Threshold")),
        pch = c(4, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white", col = c("black", "red", "black"))
    }
  } else {
    if (legend) {
      legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
        ifelse(twotail & !upperfocus,
          "Tail Fractions and Thresholds","Upper Tail Fraction and Threshold")),
        pch = c(4, rep(-1, 2)), lty = c(0, 1, 3), bg = "white", col = c("black", "red", "black"))
    }
  }
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

qplot <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, 
  legend = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "normgpdcon", "gng", "gngcon", 
    "weibullgpd",  "weibullgpdcon", "gammagpd", "gammagpdcon", "mgammagpd",
    "lognormgpd", "lognormgpdcon", "betagpd", "hpd", "hpdcon", "dwm",
    "kden", "kdengpd", "kdengpdcon", "gkg", "bckden", "bckdengpd")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
  if (!is.logical(legend))
    stop("legend must be logical")
  
  if (!is.logical(ci))
    stop("ci must be logical")
  
  if (!is.numeric(alpha) | length(alpha) != 1)
    stop("significance level alpha must be single numeric value")

  if (alpha <= 0 | alpha >= 1)
    stop("significance level alpha must be between (0, 1)")

  if (length(N) != 1 | mode(N) != "numeric") 
    stop("number of simulations must be a non-negative integer")
  
  if (abs(N - round(N)) > sqrt(.Machine$double.eps) | N < 10)
    stop("simulation size must non-negative integer >= 10")

  if (is.unsorted(modelfit$x, na.rm = TRUE))
    modelfit$x=sort(modelfit$x)
  
  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  if (modelname == "gpd") {
    upperfocus = TRUE  
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf 
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  n = length(modelfit$x)
  emp.prob=ppoints(n)

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)
  
  the.quant = switch(modelname,
    gpd = qgpd(emp.prob, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = qnormgpd(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = qnormgpdcon(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gng = qgng(emp.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = qgngcon(emp.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    weibullgpd = qweibullgpd(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = qweibullgpdcon(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gammagpd = qgammagpd(emp.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = qgammagpdcon(emp.prob, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    mgammagpd = qmgammagpd(emp.prob, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = qlognormgpd(emp.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = qlognormgpd(emp.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = qbetagpd(emp.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    hpd = qhpd(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = qhpdcon(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    dwm = qdwm(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    kden = qkden(emp.prob, modelfit$kerncentres, modelfit$lambda),
    kdengpd = qkdengpd(emp.prob, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    kdengpdcon = qkdengpdcon(emp.prob, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = qgkg(emp.prob, modelfit$kerncentres, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    bckden = qbckden(emp.prob, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = qbckdengpd(emp.prob, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))  

  axislims = range(c(modelfit$x, the.quant), na.rm = TRUE)

  if (ci) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.quantiles <- function(i, modelfit, modelname) {
      simdata = switch(modelname,
        gpd = rgpd(n, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = rnormgpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gng = rgng(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = rgngcon(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
        weibullgpd = rweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = rweibullgpdcon(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gammagpd = rgammagpd(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = rgammagpdcon(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        mgammagpd = rmgammagpd(n, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = rbetagpd(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        hpd = rhpd(n, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = rhpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        dwm = rdwm(n, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        kden = rkden(n, modelfit$kerncentres, modelfit$lambda),
        kdengpd = rkdengpd(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        kdengpdcon = rkdengpdcon(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = rgkg(n, modelfit$kerncentres, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        bckden = rbckden(n, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = rbckdengpd(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))  
      sort(simdata, na.last = FALSE)
    }
    sim.q = sapply(1:N, FUN = simulate.new.quantiles, modelfit = modelfit, modelname = modelname)
    ci.q = apply(sim.q, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    # simulation for gpd allow for uncertainty in phiu, so need to ignore CI below threshold
    if (modelname == "gpd") {
      ci.q[1, modelfit$x < u] = NA
      ci.q[2, modelfit$x < u] = NA
    }
    axislims[2] = max(c(ci.q, axislims[2]), na.rm = TRUE)
    axislims[1] = min(c(ci.q, axislims[1]), na.rm = TRUE)
  }

  if (upperfocus) {
    axislims = c(u, axislims[2])    
  }
  
  plot(the.quant, modelfit$x, pch = 'x',
    xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles", 
    xlim = axislims, ylim = axislims, ...)
  abline(c(0,1))

  if (twotail) {
    abline(h = modelfit$ul, lty = 3)
    abline(h = modelfit$ur, lty = 3)
  } else {
    abline(h = modelfit$u, lty = 3)
  }
  
  if (ci) {
    lines(the.quant, ci.q[1, ], lty=2)
    lines(the.quant, ci.q[2, ], lty=2)
    if (legend) {
      legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
        paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
        ifelse(twotail & !upperfocus, "Thresholds", "Threshold")),
        pch = c(1, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white")
    }
  } else {
    if (legend) {
      legend("bottomright", c("Data", "Fitted Model",
        ifelse(twotail & !upperfocus, "Thresholds" ,"Threshold")),
        pch = c(1, rep(-1, 2)), lty = c(0, 1, 3), bg = "white")
    }
  }
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

pplot <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000,
  legend = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "normgpdcon", "gng", "gngcon", 
    "weibullgpd",  "weibullgpdcon", "gammagpd", "gammagpdcon", "mgammagpd",
    "lognormgpd", "lognormgpdcon", "betagpd", "hpd", "hpdcon", "dwm",
    "kden", "kdengpd", "kdengpdcon", "gkg", "bckden", "bckdengpd")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
  if (!is.logical(legend))
    stop("legend must be logical")
  
  if (!is.logical(ci))
    stop("ci must be logical")
  
  if (!is.numeric(alpha) | length(alpha) != 1)
    stop("significance level alpha must be single numeric value")

  if (alpha <= 0 | alpha >= 1)
    stop("significance level alpha must be between (0, 1)")

  if (length(N) != 1 | mode(N) != "numeric") 
    stop("number of simulations must be a non-negative integer")
  
  if (abs(N - round(N)) > sqrt(.Machine$double.eps) | N < 10)
    stop("simulation size must non-negative integer >= 10")

  if (is.unsorted(modelfit$x, na.rm = TRUE))
    modelfit$x=sort(modelfit$x)
  
  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg")
  
  if (nothresh) {
    upperfocus = FALSE
  }

  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf 
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  n = length(modelfit$x)
  emp.prob = ppoints(n)
  
  the.prob = switch(modelname,
    gpd = pgpd(modelfit$x, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = pnormgpd(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = pnormgpdcon(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gng = pgng(modelfit$x, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = pgngcon(modelfit$x, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    weibullgpd = pweibullgpd(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = pweibullgpdcon(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gammagpd = pgammagpd(modelfit$x, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = pgammagpdcon(modelfit$x, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    mgammagpd = pmgammagpd(modelfit$x, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = plognormgpd(modelfit$x, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = plognormgpd(modelfit$x, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = pbetagpd(modelfit$x, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    hpd = phpd(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = phpdcon(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    dwm = pdwm(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    kden = pkden(modelfit$x, modelfit$kerncentres, modelfit$lambda),
    kdengpd = pkdengpd(modelfit$x, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    kdengpdcon = pkdengpdcon(modelfit$x, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = pgkg(modelfit$x, modelfit$kerncentres, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    bckden = pbckden(modelfit$x, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = pbckdengpd(modelfit$x, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))  

  # If focus on upper tail then must display at least upper 10%
  # even if tail fraction is smaller, otherwise may not look nice
  if (upperfocus) { 
    axislims = c(1 - ifelse(phiu < 0.1, 0.1, phiu), 1)
  } else {
    axislims = c(0, 1)
  }

  if (ci) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.probabilities <- function(i, modelfit, modelname) {
      simdata = switch(modelname,
        gpd = rgpd(n, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = rnormgpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gng = rgng(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = rgngcon(n, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
        weibullgpd = rweibullgpd(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = rweibullgpdcon(n, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gammagpd = rgammagpd(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = rgammagpdcon(n, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        mgammagpd = rmgammagpd(n, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = rlognormgpd(n, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = rbetagpd(n, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        hpd = rhpd(n, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = rhpdcon(n, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        dwm = rdwm(n, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        kden = rkden(n, modelfit$kerncentres, modelfit$lambda),
        kdengpd = rkdengpd(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        kdengpdcon = rkdengpdcon(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = rgkg(n, modelfit$kerncentres, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        bckden = rbckden(n, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = rbckdengpd(n, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))  
      simdata = sort(simdata, na.last = FALSE)
      simprob = switch(modelname,
        gpd = pgpd(simdata, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = pnormgpd(simdata, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
          modelfit$xi, modelfit$phiu),
        normgpdcon = pnormgpdcon(simdata, modelfit$nmean, modelfit$nsd, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gng = pgng(simdata, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        gngcon = pgngcon(simdata, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
        weibullgpd = pweibullgpd(simdata, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpdcon = pweibullgpdcon(simdata, modelfit$wshape, modelfit$wscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gammagpd = pgammagpd(simdata, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpdcon = pgammagpdcon(simdata, modelfit$gshape, modelfit$gscale, modelfit$u,
          modelfit$xi, modelfit$phiu),
        mgammagpd = pmgammagpd(simdata, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
          modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = plognormgpd(simdata, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpdcon = plognormgpd(simdata, modelfit$lnmean, modelfit$lnsd, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = pbetagpd(simdata, modelfit$bshape1, modelfit$bshape2, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        hpd = phpd(simdata, modelfit$nmean, modelfit$nsd, modelfit$xi),
        hpdcon = phpdcon(simdata, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
        dwm = pdwm(simdata, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
          modelfit$sigmau, modelfit$xi),
        kden = pkden(simdata, modelfit$kerncentres, modelfit$lambda),
        kdengpd = pkdengpd(simdata, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu),
        kdengpdcon = pkdengpdcon(simdata, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$xi, modelfit$phiu),
        gkg = pgkg(simdata, modelfit$kerncentres, modelfit$lambda, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
        bckden = pbckden(simdata, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
        bckdengpd = pbckdengpd(simdata, modelfit$kerncentres, modelfit$lambda, modelfit$u,
          modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
          modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))  
      simprob      
    }
    sim.p = sapply(1:N, FUN = simulate.new.probabilities, modelfit = modelfit, modelname = modelname)
    ci.p = apply(sim.p, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  }
  
  plot(the.prob, emp.prob, pch = 'x', cex = 0.8,
    xlab = "Theoretical Probability", ylab = "Empirical Probability",
    xlim = axislims, ylim = axislims, ...)
  abline(c(0,1))

  if (twotail) {
    abline(h = modelfit$phiul, lty = 3)
    abline(h = 1 - modelfit$phiur, lty = 3)
  } else {
    abline(h = 1 - modelfit$phiu, lty = 3)
  }
  
  if (ci) {
    lines(the.prob, ci.p[1, ], lty=2)
    lines(the.prob, ci.p[2, ], lty=2)
    if (legend) {
      legend("bottomright", c("Data", paste("Fitted", modelname, "Model"),
        paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
        ifelse(twotail & !upperfocus, "Tail Fractions", "Tail Fraction")),
        pch = c(1, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white")
    }
  } else {
    if (legend) {
      legend("bottomright", c("Data", "Fitted Model",
        ifelse(twotail & !upperfocus, "Tail Fractions", "Tail Fraction")),
        pch = c(1, rep(-1, 2)), lty = c(0, 1, 3), bg = "white")
    }
  }
}

#' @export
#' @aliases evmix.diag rlplot qplot pplot densplot
#' @rdname evmix.diag

densplot <- function(modelfit, upperfocus = TRUE, legend = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "normgpdcon", "gng", "gngcon", 
    "weibullgpd",  "weibullgpdcon", "gammagpd", "gammagpdcon", "mgammagpd",
    "lognormgpd", "lognormgpdcon", "betagpd", "hpd", "hpdcon", "dwm",
    "kden", "kdengpd", "kdengpdcon", "gkg", "bckden", "bckdengpd")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
  if (!is.logical(legend))
    stop("legend must be logical")
  
  if (is.unsorted(modelfit$x, na.rm = TRUE))
    modelfit$x=sort(modelfit$x)

  nothresh = (modelname == "kden") | (modelname =="bckden") | (modelname =="dwm")
  twotail = (modelname == "gng") | (modelname == "gngcon") | (modelname == "gkg")
  
  if (nothresh) {
    upperfocus = FALSE
  }
  
  # makes following code easier to read
  if (twotail) {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else if (nothresh) {
    phiu = -Inf
    u = Inf
  } else {    
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  xx = seq(min(modelfit$x, na.rm = TRUE), max(modelfit$x, na.rm = TRUE), sd(modelfit$x, na.rm = TRUE)/100)

  the.dens = switch(modelname,
    gpd = dgpd(xx, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = dnormgpd(xx, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau,
      modelfit$xi, modelfit$phiu),
    normgpdcon = dnormgpdcon(xx, modelfit$nmean, modelfit$nsd, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gng = dgng(xx, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    gngcon = dgngcon(xx, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$xil, modelfit$phiul, modelfit$ur, modelfit$xir, modelfit$phiur),
    weibullgpd = dweibullgpd(xx, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpdcon = dweibullgpdcon(xx, modelfit$wshape, modelfit$wscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gammagpd = dgammagpd(xx, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpdcon = dgammagpdcon(xx, modelfit$gshape, modelfit$gscale, modelfit$u,
      modelfit$xi, modelfit$phiu),
    mgammagpd = dmgammagpd(xx, modelfit$mgshape, modelfit$mgscale, modelfit$mgweights,
      modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = dlognormgpd(xx, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpdcon = dlognormgpd(xx, modelfit$lnmean, modelfit$lnsd, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = dbetagpd(xx, modelfit$bshape1, modelfit$bshape2, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    hpd = dhpd(xx, modelfit$nmean, modelfit$nsd, modelfit$xi),
    hpdcon = dhpdcon(xx, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$xi),
    dwm = ddwm(xx, modelfit$wshape, modelfit$wscale, modelfit$cmu, modelfit$ctau,
      modelfit$sigmau, modelfit$xi),
    kden = dkden(xx, modelfit$kerncentres, modelfit$lambda),
    kdengpd = dkdengpd(xx, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu),
    kdengpdcon = dkdengpdcon(xx, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$xi, modelfit$phiu),
    gkg = dgkg(xx, modelfit$kerncentres, modelfit$lambda, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul,
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur),
    bckden = dbckden(xx, modelfit$kerncentres, modelfit$lambda, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax),
    bckdengpd = dbckdengpd(xx, modelfit$kerncentres, modelfit$lambda, modelfit$u,
      modelfit$sigmau, modelfit$xi, modelfit$phiu, modelfit$bcmethod,
      modelfit$proper, modelfit$nn, modelfit$offset, modelfit$xmax))

  kde.est = density(modelfit$x, from = min(modelfit$x, na.rm = TRUE),
    to = max(modelfit$x, na.rm = TRUE), na.rm = TRUE)
  
  # In GPD case NA are treated as below threshold and retained from fitting
  # so need to scale KDE
  if ((modelname == "gpd") & any(is.na(modelfit$x))) {
    kde.est$y = kde.est$y * modelfit$phiu
  }
  
  if (upperfocus) {
    xlims = c(u, max(modelfit$x, na.rm = TRUE))    
    ylims = c(0, 1.4 * max(c(the.dens[xx > u], kde.est$y[kde.est$x > u])))    
  } else {
    xlims = range(modelfit$x, na.rm = TRUE)
    ylims = c(0, 1.1 * max(c(the.dens, kde.est$y)))    
  }
  
  if ((modelname == "gpd") & any(is.na(modelfit$x))) {
    hist(ifelse(is.na(modelfit$x), u - (max(modelfit$x, na.rm = TRUE) - u), modelfit$x),
      breaks = ceiling(sum(!is.na(modelfit$x))/5), freq = FALSE, 
      xlab = "Sample Data", main="", xlim = xlims, ylim = ylims)
  } else {
    hist(modelfit$x, breaks = ceiling(sum(!is.na(modelfit$x))/5), freq = FALSE, 
      xlab = "Sample Data", main="", xlim = xlims, ylim = ylims)    
  }
  rug(modelfit$x, quiet = TRUE)
  box()
    
  lines(xx, the.dens, lwd = 2)
  lines(kde.est,lty = 2, lwd = 1.5)

  if (twotail) {
    abline(v = c(modelfit$ul, modelfit$ur), lty = 3)
  } else {
    abline(v = u, lty = 3)
  }
  if (legend) {
    legend("topright", c("Fitted Density", 
    ifelse(modelname == "gng", "Thresholds", "Threshold"), "KDE"),
    lty = c(1, 3, 2), lwd = c(2, 1, 1.5), bg = "white")
  }
}
