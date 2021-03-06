% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fkden.r
\name{fkden}
\alias{fkden}
\alias{lkden}
\alias{nlkden}
\title{Cross-validation MLE Fitting of Kernel Density Estimator, With Variety of Kernels}
\usage{
fkden(x, linit = NULL, bwinit = NULL, kernel = "gaussian",
  extracentres = NULL, add.jitter = FALSE, factor = 0.1,
  amount = NULL, std.err = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...)

lkden(x, lambda = NULL, bw = NULL, kernel = "gaussian",
  extracentres = NULL, log = TRUE)

nlkden(lambda, x, bw = NULL, kernel = "gaussian",
  extracentres = NULL, finitelik = FALSE)
}
\arguments{
\item{x}{vector of sample data}

\item{linit}{initial value for bandwidth (as kernel half-width) or \code{NULL}}

\item{bwinit}{initial value for bandwidth (as kernel standard deviations) or \code{NULL}}

\item{kernel}{kernel name (\code{default = "gaussian"})}

\item{extracentres}{extra kernel centres used in KDE, 
but likelihood contribution not evaluated, or \code{NULL}}

\item{add.jitter}{logical, whether jitter is needed for rounded kernel centres}

\item{factor}{see \code{\link[base:jitter]{jitter}}}

\item{amount}{see \code{\link[base:jitter]{jitter}}}

\item{std.err}{logical, should standard errors be calculated}

\item{method}{optimisation method (see \code{\link[stats:optim]{optim}})}

\item{control}{optimisation control list (see \code{\link[stats:optim]{optim}})}

\item{finitelik}{logical, should log-likelihood return finite value for invalid parameters}

\item{...}{optional inputs passed to \code{\link[stats:optim]{optim}}}

\item{lambda}{bandwidth for kernel (as half-width of kernel) or \code{NULL}}

\item{bw}{bandwidth for kernel (as standard deviations of kernel) or \code{NULL}}

\item{log}{logical, if \code{TRUE} then log-likelihood rather than likelihood is output}
}
\value{
Log-likelihood is given by \code{\link[evmix:fkden]{lkden}} and it's
  wrappers for negative log-likelihood from \code{\link[evmix:fkden]{nlkden}}.
  Fitting function \code{\link[evmix:fkden]{fkden}} returns a simple list with the
  following elements

\tabular{ll}{
 \code{call}:        \tab \code{optim} call\cr
 \code{x}:           \tab (jittered) data vector \code{x}\cr
 \code{kerncentres}: \tab actual kernel centres used \code{x}\cr
 \code{init}:        \tab \code{linit} for lambda\cr
 \code{optim}:       \tab complete \code{optim} output\cr
 \code{mle}:         \tab vector of MLE of bandwidth\cr
 \code{cov}:         \tab variance of MLE of bandwidth\cr
 \code{se}:          \tab standard error of MLE of bandwidth\cr
 \code{nllh}:        \tab minimum negative cross-validation log-likelihood\cr
 \code{n}:           \tab total sample size\cr
 \code{lambda}:      \tab MLE of lambda (kernel half-width)\cr
 \code{bw}:          \tab MLE of bw (kernel standard deviations)\cr
 \code{kernel}:      \tab kernel name\cr
}
}
\description{
Maximum (cross-validation) likelihood estimation for fitting kernel density estimator
for a variety of possible kernels, by treating it as a mixture model.
}
\details{
The kernel density estimator (KDE) with one of possible kernels is
fitted to the entire dataset using maximum (cross-validation) likelihood estimation.
The estimated bandwidth, variance and standard error are automatically output. 

The alternate bandwidth definitions are discussed in the
\code{\link[evmix:kernels]{kernels}}, with the \code{lambda} used here but 
\code{bw} also output. The \code{bw} specification is the same as used in the
\code{\link[stats:density]{density}} function.

The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}} help
documentation with the \code{"gaussian"} as the default choice.

Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored.

Cross-validation likelihood is used for kernel density component, obtained by
leaving each point out in turn and evaluating the KDE at the point left out:
   \deqn{L(\lambda)\prod_{i=1}^{n} \hat{f}_{-i}(x_i)}
where 
   \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(n-1)\lambda} \sum_{j=1: j\ne i}^{n} K(\frac{x_i - x_j}{\lambda})}
is the KDE obtained when the \eqn{i}th datapoint is dropped out and then 
evaluated at that dropped datapoint at \eqn{x_i}.

Normally for likelihood estimation of the bandwidth the kernel centres and
the data where the likelihood is evaluated are the same. However, when using
KDE for extreme value mixture modelling the likelihood only those data in the
bulk of the distribution should contribute to the likelihood, but all the
data (including those beyond the threshold) should contribute to the density
estimate. The \code{extracentres} option allows the use to specify extra
kernel centres used in estimating the density, but not evaluated in the
likelihood. Suppose the first \code{nb} data are below the threshold, followed
by \code{nu} exceedances of the threshold, so \eqn{i = 1,\ldots,nb, nb+1, \ldots, nb+nu}.
The cross-validation likelihood using the extra kernel centres is then:
   \deqn{L(\lambda)\prod_{i=1}^{nb} \hat{f}_{-i}(x_i)}
where 
   \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(nb+nu-1)\lambda} \sum_{j=1: j\ne i}^{nb+nu} K(\frac{x_i - x_j}{\lambda})}
which shows that the complete set of data is used in evaluating the KDE, but only those
below the threshold contribute to the cross-validation likelihood. The default is to
use the existing data, so \code{extracentres=NULL}.

The following functions are provided:
\itemize{
 \item \code{\link[evmix:fkden]{fkden}} - maximum (cross-validation) likelihood fitting with all the above options;
 \item \code{\link[evmix:fkden]{lkden}} - cross-validation log-likelihood;
 \item \code{\link[evmix:fkden]{nlkden}} - negative cross-validation log-likelihood;
}
The log-likelihood functions are provided for wider usage, e.g. constructing
profile likelihood functions.

The log-likelihood and negative log-likelihood are also provided for wider
usage, e.g. constructing your own extreme value
mixture models or profile likelihood functions. The parameter
\code{lambda} must be specified in the negative log-likelihood
\code{\link[evmix:fkden]{nlkden}}.

Log-likelihood calculations are carried out in
\code{\link[evmix:fkden]{lkden}}, which takes bandwidths as inputs in
the same form as distribution functions. The negative log-likelihood is a
wrapper for \code{\link[evmix:fkden]{lkden}}, designed towards making
it useable for optimisation (e.g. \code{lambda} given as first input).

Defaults values for the bandwidth \code{linit} and \code{lambda} are given in the fitting 
\code{\link[evmix:fkden]{fkden}} and cross-validation likelihood functions
\code{\link[evmix:fkden]{lkden}}. The bandwidth \code{linit} must be specified in
the negative log-likelihood function \code{\link[evmix:fkden]{nlkden}}. 

Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored,
which is inconsistent with the \code{\link[evd:fpot]{evd}} library which assumes the 
missing values are below the threshold.

The function \code{\link[evmix:fkden]{lkden}} carries out the calculations
for the log-likelihood directly, which can be exponentiated to give actual
likelihood using (\code{log=FALSE}).

The default optimisation algorithm is "BFGS", which requires a finite negative 
log-likelihood function evaluation \code{finitelik=TRUE}. For invalid 
parameters, a zero likelihood is replaced with \code{exp(-1e6)}. The "BFGS" 
optimisation algorithms require finite values for likelihood, so any user 
input for \code{finitelik} will be overridden and set to \code{finitelik=TRUE} 
if either of these optimisation methods is chosen.

It will display a warning for non-zero convergence result comes from 
\code{\link[stats:optim]{optim}} function call or for common indicators of lack
of convergence (e.g. estimated bandwidth equal to initial value).

If the hessian is of reduced rank then the variance covariance (from inverse hessian)
and standard error of parameters cannot be calculated, then by default 
\code{std.err=TRUE} and the function will stop. If you want the parameter estimates
even if the hessian is of reduced rank (e.g. in a simulation study) then
set \code{std.err=FALSE}.
}
\note{
When \code{linit=NULL} then the initial value for the \code{lambda}
bandwidth is calculated 
using \code{\link[stats:bandwidth]{bw.nrd0}} function and transformed using 
\code{\link[evmix:kfun]{klambda}} function.

The extra kernel centres \code{extracentres} can either be a vector of data or \code{NULL}.

Invalid parameter ranges will give \code{0} for likelihood, \code{log(0)=-Inf} for
log-likelihood and \code{-log(0)=Inf} for negative log-likelihood. 

Infinite and missing sample values are dropped.

Error checking of the inputs is carried out and will either stop or give warning message
as appropriate.
}
\section{Warning}{

Two important practical issues arise with MLE for the kernel bandwidth:
1) Cross-validation likelihood is needed for the KDE bandwidth parameter
as the usual likelihood degenerates, so that the MLE \eqn{\hat{\lambda} \rightarrow 0} as
\eqn{n \rightarrow \infty}, thus giving a negative bias towards a small bandwidth.
Leave one out cross-validation essentially ensures that some smoothing between the kernel centres
is required (i.e. a non-zero bandwidth), otherwise the resultant density estimates would always
be zero if the bandwidth was zero.

This problem occassionally rears its ugly head for data which has been heavily rounded,
as even when using cross-validation the density can be non-zero even if the bandwidth is zero.
To overcome this issue an option to add a small jitter should be added to the data
(\code{x} only) has been included in the fitting inputs, using the 
\code{\link[base:jitter]{jitter}} function, to remove the ties. The default options red in the 
\code{\link[base:jitter]{jitter}} are specified above, but the user can override these.
Notice the default scaling \code{factor=0.1}, which is a tenth of the default value in the
\code{\link[base:jitter]{jitter}}
function itself.

A warning message is given if the data appear to be rounded
(i.e. more than 5% of data are tied). If the estimated bandwidth is too small, then
data rounding is the likely culprit. Only use the jittering when the MLE of
the bandwidth is far too small. 

2) For heavy tailed populations the bandwidth is positively biased, giving oversmoothing
(see example). The bias is due to the distance between the upper (or lower) order statistics not
necessarily decaying to zero as the sample size tends to infinity. Essentially, as the distance
between the two largest (or smallest) sample datapoints does not decay to zero, some smoothing between
them is required (i.e. bandwidth cannot be zero). One solution to this problem is to trim
the data at a suitable threshold to remove the problematic tail from the inference for the bandwidth, 
using either the \code{\link[evmix:fkdengpd]{fkdengpd}} function for a single heavy tail
or the \code{\link[evmix:fgkg]{fgkg}} function
if both tails are heavy. See MacDonald et al (2013).
}

\section{Acknowledgments}{
 See Acknowledgments in
  \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}. Based on code
by Anna MacDonald produced for MATLAB.
}

\examples{
\dontrun{
set.seed(1)
par(mfrow = c(1, 1))

nk=50
x = rnorm(nk)
xx = seq(-5, 5, 0.01)
fit = fkden(x)
hist(x, nk/5, freq = FALSE, xlim = c(-5, 5), ylim = c(0,0.6)) 
rug(x)
for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit$lambda)*0.05)
lines(xx,dnorm(xx), col = "black")
lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
lines(density(x, bw = fit$bw), lwd = 2, lty = 2,  col = "blue")
legend("topright", c("True Density", "KDE fitted evmix",
"KDE Using density, default bandwidth", "KDE Using density, c-v likelihood bandwidth"),
lty = c(1, 1, 2, 2), lwd = c(1, 2, 2, 2), col = c("black", "red", "green", "blue"))

par(mfrow = c(2, 1))

# bandwidth is biased towards oversmoothing for heavy tails
nk=100
x = rt(nk, df = 2)
xx = seq(-8, 8, 0.01)
fit = fkden(x)
hist(x, seq(floor(min(x)), ceiling(max(x)), 0.5), freq = FALSE, xlim = c(-8, 10)) 
rug(x)
for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit$lambda)*0.05)
lines(xx,dt(xx , df = 2), col = "black")
lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "red")
legend("topright", c("True Density", "KDE fitted evmix, c-v likelihood bandwidth"),
lty = c(1, 1), lwd = c(1, 2), col = c("black", "red"))

# remove heavy tails from cv-likelihood evaluation, but still include them in KDE within likelihood
# often gives better bandwidth (see MacDonald et al (2011) for justification)
nk=100
x = rt(nk, df = 2)
xx = seq(-8, 8, 0.01)
fit2 = fkden(x[(x > -4) & (x < 4)], extracentres = x[(x <= -4) | (x >= 4)])
hist(x, seq(floor(min(x)), ceiling(max(x)), 0.5), freq = FALSE, xlim = c(-8, 10)) 
rug(x)
for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit2$lambda)*0.05)
lines(xx,dt(xx , df = 2), col = "black")
lines(xx, dkden(xx, x, lambda = fit2$lambda), lwd = 2, col = "red")
lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "blue")
legend("topright", c("True Density", "KDE fitted evmix, tails removed",
"KDE fitted evmix, tails included"),
lty = c(1, 1, 1), lwd = c(1, 2, 2), col = c("black", "red", "blue"))
}

}
\references{
\url{http://en.wikipedia.org/wiki/Kernel_density_estimation}

\url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}

Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
threshold estimation and uncertainty quantification. REVSTAT - Statistical
Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}

Hu Y. and Scarrott, C.J. (2018). evmix: An R Package for Extreme Value Mixture Modeling, 
Threshold Estimation and Boundary Corrected Kernel Density Estimation. Journal of
Statistical Software 84(5), 1-27. doi: 10.18637/jss.v084.i05.

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
\code{\link[base:jitter]{jitter}}, \code{\link[stats:density]{density}} and
\code{\link[stats:bandwidth]{bw.nrd0}}

Other kden: \code{\link{bckden}}, \code{\link{fbckden}},
  \code{\link{fgkgcon}}, \code{\link{fgkg}},
  \code{\link{fkdengpdcon}}, \code{\link{fkdengpd}},
  \code{\link{kdengpdcon}}, \code{\link{kdengpd}},
  \code{\link{kden}}

Other kdengpd: \code{\link{bckdengpd}},
  \code{\link{fbckdengpd}}, \code{\link{fgkg}},
  \code{\link{fkdengpdcon}}, \code{\link{fkdengpd}},
  \code{\link{gkg}}, \code{\link{kdengpdcon}},
  \code{\link{kdengpd}}, \code{\link{kden}}

Other bckden: \code{\link{bckdengpdcon}},
  \code{\link{bckdengpd}}, \code{\link{bckden}},
  \code{\link{fbckdengpdcon}}, \code{\link{fbckdengpd}},
  \code{\link{fbckden}}, \code{\link{kden}}

Other fkden: \code{\link{kden}}
}
\author{
Yang Hu and Carl Scarrott \email{carl.scarrott@canterbury.ac.nz}.
}
\concept{bckden}
\concept{fkden}
\concept{kden}
\concept{kdengpd}
