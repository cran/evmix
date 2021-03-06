% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fbckdengpd.r
\name{fbckdengpd}
\alias{fbckdengpd}
\alias{lbckdengpd}
\alias{nlbckdengpd}
\alias{proflubckdengpd}
\alias{nlubckdengpd}
\title{MLE Fitting of Boundary Corrected Kernel Density Estimate for Bulk and GPD Tail Extreme Value Mixture Model}
\usage{
fbckdengpd(x, phiu = TRUE, useq = NULL, fixedu = FALSE,
  pvector = NULL, kernel = "gaussian", bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL,
  add.jitter = FALSE, factor = 0.1, amount = NULL, std.err = TRUE,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE,
  ...)

lbckdengpd(x, lambda = NULL, u = 0, sigmau = 1, xi = 0,
  phiu = TRUE, bw = NULL, kernel = "gaussian", bcmethod = "simple",
  proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL,
  log = TRUE)

nlbckdengpd(pvector, x, phiu = TRUE, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL,
  xmax = NULL, finitelik = FALSE)

proflubckdengpd(u, pvector, x, phiu = TRUE, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL,
  xmax = NULL, method = "BFGS", control = list(maxit = 10000),
  finitelik = TRUE, ...)

nlubckdengpd(pvector, u, x, phiu = TRUE, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL,
  xmax = NULL, finitelik = FALSE)
}
\arguments{
\item{x}{vector of sample data}

\item{phiu}{probability of being above threshold \eqn{(0, 1)} or logical, see Details in 
help for \code{\link[evmix:fnormgpd]{fnormgpd}}}

\item{useq}{vector of thresholds (or scalar) to be considered in profile likelihood or
\code{NULL} for no profile likelihood}

\item{fixedu}{logical, should threshold be fixed (at either scalar value in \code{useq},
or estimated from maximum of profile likelihood evaluated at
sequence of thresholds in \code{useq})}

\item{pvector}{vector of initial values of parameters or \code{NULL} for default
values, see below}

\item{kernel}{kernel name (\code{default = "gaussian"})}

\item{bcmethod}{boundary correction method}

\item{proper}{logical, whether density is renormalised to integrate to unity (where needed)}

\item{nn}{non-negativity correction method (simple boundary correction only)}

\item{offset}{offset added to kernel centres (logtrans only) or \code{NULL}}

\item{xmax}{upper bound on support (copula and beta kernels only) or \code{NULL}}

\item{add.jitter}{logical, whether jitter is needed for rounded kernel centres}

\item{factor}{see \code{\link[base:jitter]{jitter}}}

\item{amount}{see \code{\link[base:jitter]{jitter}}}

\item{std.err}{logical, should standard errors be calculated}

\item{method}{optimisation method (see \code{\link[stats:optim]{optim}})}

\item{control}{optimisation control list (see \code{\link[stats:optim]{optim}})}

\item{finitelik}{logical, should log-likelihood return finite value for invalid parameters}

\item{...}{optional inputs passed to \code{\link[stats:optim]{optim}}}

\item{lambda}{bandwidth for kernel (as half-width of kernel) or \code{NULL}}

\item{u}{scalar threshold value}

\item{sigmau}{scalar scale parameter (positive)}

\item{xi}{scalar shape parameter}

\item{bw}{bandwidth for kernel (as standard deviations of kernel) or \code{NULL}}

\item{log}{logical, if \code{TRUE} then log-likelihood rather than likelihood is output}
}
\value{
\code{\link[evmix:fbckdengpd]{lbckdengpd}}, \code{\link[evmix:fbckdengpd]{nlbckdengpd}},
and \code{\link[evmix:fbckdengpd]{nlubckdengpd}} give the log-likelihood,
negative log-likelihood and profile likelihood for threshold. Profile likelihood
for single threshold is given by \code{\link[evmix:fbckdengpd]{proflubckdengpd}}.
\code{\link[evmix:fbckdengpd]{fbckdengpd}} returns a simple list with the following elements

\tabular{ll}{
 \code{call}:      \tab \code{optim} call\cr
 \code{x}:         \tab data vector \code{x}\cr
 \code{init}:      \tab \code{pvector}\cr
 \code{fixedu}:    \tab fixed threshold, logical\cr
 \code{useq}:      \tab threshold vector for profile likelihood or scalar for fixed threshold\cr
 \code{nllhuseq}:  \tab profile negative log-likelihood at each threshold in useq\cr
 \code{optim}:     \tab complete \code{optim} output\cr
 \code{mle}:       \tab vector of MLE of parameters\cr
 \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
 \code{se}:        \tab vector of standard errors of MLE of parameters\cr
 \code{rate}:      \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
 \code{nllh}:      \tab minimum negative log-likelihood\cr
 \code{n}:         \tab total sample size\cr
 \code{lambda}:    \tab MLE of lambda (kernel half-width)\cr
 \code{u}:         \tab threshold (fixed or MLE)\cr
 \code{sigmau}:    \tab MLE of GPD scale\cr
 \code{xi}:        \tab MLE of GPD shape\cr
 \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
 \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
 \code{bw}:        \tab MLE of bw (kernel standard deviations)\cr
 \code{kernel}:    \tab kernel name\cr
 \code{bcmethod}:  \tab boundary correction method\cr
 \code{proper}:    \tab logical, whether renormalisation is requested\cr
 \code{nn}:        \tab non-negative correction method\cr
 \code{offset}:    \tab offset for log transformation method\cr
 \code{xmax}:      \tab maximum value of scaled beta or copula
}
}
\description{
Maximum likelihood estimation for fitting the extreme value 
mixture model with boundary corrected kernel density estimate for bulk distribution upto the threshold and conditional
GPD above threshold. With options for profile likelihood estimation for threshold and
fixed threshold approach.
}
\details{
The extreme value mixture model with boundary corrected kernel density estimate (BCKDE) for bulk and GPD tail is 
fitted to the entire dataset using maximum likelihood estimation. The estimated
parameters, variance-covariance matrix and their standard errors are automatically
output.

See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. 
Only the different features are outlined below for brevity.

The full parameter vector is
(\code{lambda}, \code{u}, \code{sigmau}, \code{xi}) if threshold is also estimated and
(\code{lambda}, \code{sigmau}, \code{xi}) for profile likelihood or fixed threshold approach.

Negative data are ignored.

Cross-validation likelihood is used for BCKDE, but standard likelihood is used
for GPD component. See help for \code{\link[evmix:fkden]{fkden}} for details,
type \code{help fkden}.

The alternate bandwidth definitions are discussed in the 
\code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default
used in the likelihood fitting. The \code{bw} specification is the same as
used in the \code{\link[stats:density]{density}} function.

The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}}
with the \code{"gaussian"} as the default choice.

Unlike the standard KDE, there is no general rule-of-thumb bandwidth for all these
estimators, with only certain methods having a guideline in the literature, so none
have been implemented. Hence, a bandwidth must always be specified.

The \code{simple}, \code{renorm}, \code{beta1}, \code{beta2} \code{gamma1} and \code{gamma2}
boundary corrected kernel density estimates require renormalisation, achieved
by numerical integration, so are very time consuming.
}
\note{
See notes in \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}.
Only the different features are outlined below for brevity.

No default initial values for parameter vector are provided, so will stop evaluation if
\code{pvector} is left as \code{NULL}. Avoid setting the starting value for the shape parameter to
\code{xi=0} as depending on the optimisation method it may be get stuck.

The data and kernel centres are both vectors. Infinite, missing and negative sample values
(and kernel centres) are dropped.
}
\section{Boundary Correction Methods}{


See \code{\link[evmix:bckden]{dbckden}} for details of BCKDE methods.
}

\section{Warning}{

See important warnings about cross-validation likelihood estimation in 
\code{\link[evmix:fkden]{fkden}}, type \code{help fkden}.

See important warnings about boundary correction approaches in 
\code{\link[evmix:bckden]{dbckden}}, type \code{help bckden}.
}

\section{Acknowledgments}{
 See Acknowledgments in
  \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}. Based on code
by Anna MacDonald produced for MATLAB.
}

\examples{
\dontrun{
set.seed(1)
par(mfrow = c(2, 1))

x = rgamma(500, 2, 1)
xx = seq(-0.1, 10, 0.01)
y = dgamma(xx, 2, 1)

# Bulk model based tail fraction
pinit = c(0.1, quantile(x, 0.9), 1, 0.1) # initial values required for BCKDE
fit = fbckdengpd(x, pvector = pinit, bcmethod = "cutnorm")
hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 10))
lines(xx, y)
with(fit, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bcmethod = "cutnorm"), col="red"))
abline(v = fit$u, col = "red")
  
# Parameterised tail fraction
fit2 = fbckdengpd(x, phiu = FALSE, pvector = pinit, bcmethod = "cutnorm")
with(fit2, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, phiu, bc = "cutnorm"), col="blue"))
abline(v = fit2$u, col = "blue")
legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
  col=c("black", "red", "blue"), lty = 1)
  
# Profile likelihood for initial value of threshold and fixed threshold approach
pinit = c(0.1, 1, 0.1) # notice threshold dropped from initial values
fitu = fbckdengpd(x, useq = seq(1, 6, length = 20), pvector = pinit, bcmethod = "cutnorm")
fitfix = fbckdengpd(x, useq = seq(1, 6, length = 20), fixedu = TRUE, pv = pinit, bc = "cutnorm")

hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 10))
lines(xx, y)
with(fit, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bc = "cutnorm"), col="red"))
abline(v = fit$u, col = "red")
with(fitu, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bc = "cutnorm"), col="purple"))
abline(v = fitu$u, col = "purple")
with(fitfix, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bc = "cutnorm"), col="darkgreen"))
abline(v = fitfix$u, col = "darkgreen")
legend("topright", c("True Density","Default initial value (90\% quantile)",
 "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
 col=c("black", "red", "purple", "darkgreen"), lty = 1)
}

}
\references{
\url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}

\url{http://en.wikipedia.org/wiki/Kernel_density_estimation}

\url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}

\url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}

Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
threshold estimation and uncertainty quantification. REVSTAT - Statistical
Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

Hu, Y. (2013). Extreme value mixture modelling: An R package and simulation study.
MSc (Hons) thesis, University of Canterbury, New Zealand.
\url{http://ir.canterbury.ac.nz/simple-search?query=extreme&submit=Go}

Bowman, A.W. (1984). An alternative method of cross-validation for the smoothing of
density estimates. Biometrika 71(2), 353-360.

Duin, R.P.W. (1976). On the choice of smoothing parameters for Parzen estimators of
probability density functions. IEEE Transactions on Computers C25(11), 1175-1179.

MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B., Reale, M. and Russell, G. (2011).
A flexible extreme value mixture model. Computational Statistics and Data Analysis
55(6), 2137-2157.

MacDonald, A., C. J. Scarrott, and D. S. Lee (2011). Boundary correction, consistency
and robustness of kernel densities using extreme value theory. Submitted.
Available from: \url{http://www.math.canterbury.ac.nz/~c.scarrott}.

Wand, M. and Jones, M.C. (1995). Kernel Smoothing. Chapman && Hall.
}
\seealso{
\code{\link[evmix:kernels]{kernels}}, \code{\link[evmix:kfun]{kfun}},
 \code{\link[stats:density]{density}}, \code{\link[stats:bandwidth]{bw.nrd0}}
and \code{\link[ks:kde]{dkde}} in \code{\link[ks:kde]{ks}} package.
 \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}.

Other kdengpd: \code{\link{bckdengpd}}, \code{\link{fgkg}},
  \code{\link{fkdengpdcon}}, \code{\link{fkdengpd}},
  \code{\link{fkden}}, \code{\link{gkg}},
  \code{\link{kdengpdcon}}, \code{\link{kdengpd}},
  \code{\link{kden}}

Other bckden: \code{\link{bckdengpdcon}},
  \code{\link{bckdengpd}}, \code{\link{bckden}},
  \code{\link{fbckdengpdcon}}, \code{\link{fbckden}},
  \code{\link{fkden}}, \code{\link{kden}}

Other bckdengpd: \code{\link{bckdengpdcon}},
  \code{\link{bckdengpd}}, \code{\link{bckden}},
  \code{\link{fbckdengpdcon}}, \code{\link{fbckden}},
  \code{\link{fkdengpd}}, \code{\link{gkg}},
  \code{\link{kdengpd}}, \code{\link{kden}}

Other bckdengpdcon: \code{\link{bckdengpdcon}},
  \code{\link{bckdengpd}}, \code{\link{bckden}},
  \code{\link{fbckdengpdcon}}, \code{\link{fbckden}},
  \code{\link{fkdengpdcon}}, \code{\link{gkgcon}},
  \code{\link{kdengpdcon}}

Other fbckdengpd: \code{\link{bckdengpd}}
}
\author{
Yang Hu and Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}
}
\concept{bckden}
\concept{bckdengpd}
\concept{bckdengpdcon}
\concept{fbckdengpd}
\concept{kdengpd}
