#shell('R CMD Rd2pdf man')

#' @examples
#' par(mfrow=c(2,1))
#' x = rgamma(1000, shape = 2)
#' xx = seq(-1, 10, 0.01)
#' y = dgamma(xx, shape = 2)
#' 
#' # Bulk model base tail fraction
#' fit = fmgammagpd(x, phiu = TRUE, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, y)
#' lines(xx, dmgammagpd(xx, mgshape = fit$mgshape, mgscale = fit$mgscale, u = fit$u,
#'   sigmau = fit$sigmau, xi = fit$xi, phiu = TRUE), col="red")
#' abline(v = fit$u)
#'   
#' # Parameterised tail fraction
#' fit2 = fmgammagpd(x, phiu = FALSE, std.err = FALSE)
#' plot(xx, y, type = "l")
#' lines(xx, dmgammagpd(xx, mgshape = fit$mgshape, mgscale = fit$mgscale, u = fit$u,
#'   sigmau = fit$sigmau, xi = fit$xi, phiu = TRUE), col="red")
#' lines(xx, dmgammagpd(xx, mgshape = fit2$mgshape, mgscale = fit2$mgscale, u = fit2$u,
#'   sigmau = fit2$sigmau, xi = fit2$xi, phiu = fit2$phiu), col="blue")
#' abline(v = fit$u, col = "red")
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)

#' @examples
#' par(mfrow=c(2,2))
#' x = rmgammagpd(1000, gshape = 2)
#' xx = seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dmgammagpd(xx, gshape = 2))
#' 
#' # three tail behaviours
#' plot(xx, pmgammagpd(xx, gshape = 2), type = "l")
#' lines(xx, pmgammagpd(xx, gshape = 2, xi = 0.3), col = "red")
#' lines(xx, pmgammagpd(xx, gshape = 2, xi = -0.3), col = "blue")
#' legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rmgammagpd(1000, gshape = 2, u = 3, phiu = 0.2)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dmgammagpd(xx, gshape = 2, u = 3, phiu = 0.2))
#' 
#' plot(xx, dmgammagpd(xx, gshape = 2, u = 3, xi=0, phiu = 0.2), type = "l")
#' lines(xx, dmgammagpd(xx, gshape = 2, u = 3, xi=-0.2, phiu = 0.2), col = "red")
#' lines(xx, dmgammagpd(xx, gshape = 2, u = 3, xi=0.2, phiu = 0.2), col = "blue")
#' legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)


if (){
    print(c(gshape,gscale,u,sigmau,xi,phiu,phib, dgamma(u, shape = gshape, scale = gscale),min(syu)))
}

  dbckdencvi <- function(i, isamp, x, lambda, bcmethod, proper, nn, offset, xmax) {
    sum(dbckden(x[which(isamp == i)], kerncentres = x[-which(isamp == i)],
    lambda = lambda, bcmethod = "simple",
    proper = proper, nn = nn, offset = offset, xmax = xmax, log = TRUE))
  }
  
  dbckdeni <- function(i, x, lambda, bcmethod, proper, nn, offset, xmax) {
    dbckden(x[i], kerncentres = x[-i], lambda = lambda, bcmethod = bcmethod,
    proper = proper, nn = nn, offset = offset, xmax = xmax, log = TRUE)
  }
      
  if ((lambda <= 0) | ((bcmethod == "copula") & (lambda >= 1))) {
    l = -Inf
   } else {
     if (bcmethod == "simple") {
       k = 50
       isamp = c(sample(rep(1:k, each = floor(n/k)), k*floor(n/k), replace = FALSE),
         sample(1:k, n - k*floor(n/k), replace = TRUE))

       l = sum(sapply(1:k, FUN = dbckdencvi, isamp = isamp, x = kerncentres, lambda = lambda,
         bcmethod = "simple", proper = proper, nn = nn, offset = offset, xmax = xmax))
     } else {
       l = sum(sapply(1:n, FUN = dbckdeni, x = kerncentres, lambda = lambda,
         bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax))
     }
   }
   

  which0 = which(p < min(pk))
  n0 = length(which0)
  which1 = which(p > max(pk))
  n1 = length(which1)

# extrapolation has numerical solution (slow but not many of these)
  
  # find q to minimise difference between p and pkden(q)
  pkdexmin <- function(q, prob, kerncentres, lambda) {
    abs(mean(pnorm(q, mean = kerncentres, sd = lambda)) - prob)
  }

  findqkdexmin = function(prob, kerncentres, lambda, qinit) {
    q = try(nlm(pkdexmin, qinit, prob = prob, lambda = lambda, kerncentres = kerncentres,
      gradtol = 1e-10, steptol = 1e-10)$estimate)

    if (inherits(q, "try-error")) {
      q = NA
    }
    q
  }
  
  # extrapolate lower tail
  for (i in which0) {
    q[i] = findqkdexmin(p[i], kerncentres = kerncentres, lambda = lambda, qinit = qk[1])
  }

  # extrapolate upper tail
  for (i in which1) {
    q[i] = findqkdexmin(p[i], kerncentres = kerncentres, lambda = lambda, qinit = qk[length(qk)])
  }



n=10
x = rgamma(n, shape = 1, scale = 1)
xx = seq(-0.5, 15, 0.01)
plot(xx, pgamma(xx, shape = 1, scale = 1), type = "l", ylim=c(0,1.1),lwd=2)
rug(x);abline(v=0.6)
temp=pbckden(xx, x, lambda = 0.1, bcmethod = "simple", proper = TRUE, nn="jf96");lines(xx, temp, lwd = 2, col = "red");abline(v=range(x))
lines(qbckden(temp, x, lambda = 0.1, bcmethod = "simple", proper = TRUE, nn="jf96"),temp,lty=2,col="blue")

lambdas=seq(0,8,0.01)
ll=sapply(lambdas, FUN=nlbckden, x=x, bcmethod = "simple", proper = TRUE, nn = "jf96")
plot(lambdas, ll)


  # CV-likelihood bandwidth estimator is biased if upper tail included
  if ((bcmethod == "simple") & is.null(extracentres)) {
    q75 = quantile(x, 0.75)
    indata = x[x <= q75]
    extracentres = x[x > q75]
  } else {
    indata = x
  }
  






nk=5000; bw=0.005; sh=0.5; sc=1;meth="beta2"
x = rbeta(nk, shape1 = sh, shape2 = sc)*5
xx = seq(-0.1, 5.1, 0.001)
plot(xx, dbeta(xx/5, shape1 = sh, shape2 = sc)/5, type="l", ylim = c(0, 5))
rug(x)
for (i in 1:min(nk,20)) lines(xx, dbckden(xx/5, x[i], lambda = bw, bcmethod = meth, xmax = 5, proper = FALSE)*0.2/5, col = "blue")
lines(xx, dbckden(xx, x, lambda = bw, bcmethod = meth, xmax = 5, proper = TRUE), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
legend("topright", c("True Density", "Modified Beta Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))

nk=50; bw=0.1; sh=2; sc=3;meth="beta2"
x = rbeta(nk, shape1 = sh, shape2 = sc)*5
xx = seq(-0.1, 5.1, 0.01)
plot(xx, pbeta(xx/5, shape1 = sh, shape2 = sc), type="l", ylim = c(0, 1))
rug(x)
temp=pbckden(xx, x, lambda = bw, bcmethod = meth, xmax = 5, proper = TRUE)
lines(xx, temp, lwd = 2, col = "red")
lines(qbckden(temp, x, lambda = bw, bcmethod = meth, xmax = 5, proper = TRUE), temp, lwd = 2, lty=2, col = "blue")
legend("topleft", c("True Density", "Modified Beta Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))




nk=5000; bw=0.1; sh=4; sc=2;meth="copula"
nk=5000; bw=0.1; sh=1; sc=1;meth="copula"
nk=5000; bw=0.1; sh=0.5; sc=1;meth="copula"
x = rbeta(nk, shape1 = sh, shape2 = sc)*5
xx = seq(-0.1, 5.1, 0.001)
plot(xx, dbeta(xx/5, shape1 = sh, shape2 = sc)/5, ylim = c(0, 3))
rug(x)
for (i in 1:min(nk,20)) lines(xx, dbckden(xx, x[i], lambda=bw,method=meth,xmax=5)*0.05/5, col = "blue")
lines(xx, dbckden(xx, x, method=meth, xmax=5,proper=TRUE), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
legend("topright", c("True Density", "BC KDE Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))


n=1000;x=rt(n,df=6);hist(x,n/5, freq=FALSE,xlim=c(-6,6))
fit=fkden(x);xx=seq(-6,6,0.01);lines(xx,dkden(xx,x,lambda=fit$lambda),col="blue",lwd=2)
lines(xx,dt(xx,df=6),col="black")
rug(x);for (i in 1:length(x)){lines(xx,dnorm(xx,x[i],sd=fit$lambda)*0.05)}
lines(density(x),col="green",lwd=1.5)
lines(density(x,bw=fit$lambda),col="red",lwd=2,lty=2)

xx=seq(-1,10,0.01);y=0;bw=0.5;x11();
plot(xx,dbckden(xx, y, lambda = bw,proper=TRUE),xlim=range(xx));abline(v=0);
y=0.5;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="red");
lines(xx,dnorm(xx,y,bw)/pnorm(y,sd=bw),col="blue",lwd=2);
y=1;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="purple");
y=2;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="orange")
y=3;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="cyan")

xx=seq(-8,8,0.01);
y=0;bw=2;x11();plot(-(xx-y),dbckden(xx, y, lambda = bw,proper=TRUE),xlim=range(xx));
lines(xx,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="blue",lwd=2);abline(v=0)
y=0.5;lines(-(xx-y),dbckden(xx, y, lambda = bw,proper=TRUE),col="red",lwd=2);
lines(xx,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="blue",lwd=1,lty=2)
y=1;lines(-(xx-y),dbckden(xx, y, lambda = bw,proper=TRUE),col="purple",lwd=2);
lines(xx,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="blue",lwd=1,lty=2)
y=2;lines(-(xx-y),dbckden(xx, y, lambda = bw,proper=TRUE),col="orange",lwd=2);
lines(xx,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="orange",lwd=1,lty=2)


nk=5000; bw=0.2; sh=3; sc=2;meth="simple"
x = rchisq(nk,df=5)
xx = seq(-1, 12, 0.01)
plot(xx, dchisq(xx,df=5), ylim = c(0, 1))
rug(x)
for (i in 1:min(nk,50)) lines(xx, dbckden(xx, x[i], lambda = bw, method=meth)*0.05, col = "blue")
lines(xx, dbckden(xx, x, lambda = bw, method=meth), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
legend("topright", c("True Density", "BC KDE Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))

xx=seq(-8,8,0.01);
y=0;bw=2;x11();plot(xx,dbckden(xx, y, lambda = bw,proper=TRUE),xlim=range(xx));
lines(xx+y,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="blue",lwd=2);abline(v=0)
y=0.5;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="red",lwd=2);
lines(xx+y,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="blue",lwd=1,lty=2)
y=1;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="purple",lwd=2);
lines(xx+y,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="blue",lwd=1,lty=2)
y=2;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="orange",lwd=2);
lines(xx+y,(pnorm(y/bw)+(xx/bw-y/bw)*dnorm(y/bw))*dnorm(xx/bw)/bw/(pnorm(y/bw)*(pnorm(y/bw)-y*dnorm(y/bw)/bw)-dnorm(y/bw)^2),col="orange",lwd=1,lty=2)

y=0;bw=2;x11();plot(xx,dbckden(xx, y, lambda = bw,proper=TRUE),xlim=range(xx));
y=0.5;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="red",lwd=2);
y=1;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="purple",lwd=2);
y=2;lines(xx,dbckden(xx, y, lambda = bw,proper=TRUE),col="orange",lwd=2);

nk=5000
bw=1
x = rgamma(nk, shape = 1,scale=2)
xx = seq(-1, 12, 0.01)
plot(xx, dgamma(xx, shape = 1,scale=2), ylim = c(0, 1))
rug(x)
for (i in 1:50) lines(xx, dbckden(xx, x[i], lambda = bw,proper=FALSE)*0.05, col = "blue")
lines(xx, dbckden(xx, x, lambda = bw), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
legend("topright", c("True Density", "BC KDE Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))
lines(xx, dbckden(xx, x, lambda = 0.2, method="renorm"), lwd = 2, col = "orange")

nk=50; bw=0.5; sh=3; sc=2;meth="simple"
x = rchisq(nk,df=2)
xx = seq(-1, 12, 0.01)
plot(xx, dchisq(xx,df=2), ylim = c(0, 1))
rug(x)
for (i in 1:min(nk,50)) lines(xx, dbckden(xx, x[i], lambda = bw, method=meth)*0.05, col = "blue")
lines(xx, dbckden(xx, x, lambda = bw, method=meth), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
legend("topright", c("True Density", "BC KDE Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))

nk=100
bw=0.5
x = rgamma(nk, shape = 1,scale=2)
xx = seq(-1, 12, 0.01)
plot(xx, dgamma(xx, shape = 1,scale=2), ylim = c(0, 1))
rug(x)
for (i in 1:nk) lines(xx, dbckden(xx, x[i], lambda = bw,proper=FALSE)*0.05, col = "blue")
lines(xx, dbckden(xx, x, lambda = bw), lwd = 2, col = "red")
lines(density(x), lty = 2, lwd = 2, col = "green")
legend("topright", c("True Density", "BC KDE Using evmix", "KDE Using density function",
"Boundary Corrected Kernels"),
lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))

n=500
x = rgamma(n, shape = 1, scale = 2)
xx = seq(-0.1, 15, 0.01)
plot(xx, pgamma(xx, shape = 1, scale = 2), type = "l")
rug(x);abline(v=0.6)
temp=pbckden(xx, x, lambda = 0.3, method = "simple", proper = TRUE,nn="jf96");lines(xx, temp, lwd = 2, col = "red")
lines(qbckden(temp, x, lambda = 0.3, method = "simple", proper = TRUE,nn="jf96"),temp,lty=2,col="blue")

n=500
x = rgamma(n, shape = 3, scale = 2)
xx = seq(-0.1, 15, 0.1)
plot(xx, pgamma(xx, shape = 3, scale = 2), type = "l")
rug(x);abline(v=0.6)
temp=pbckden(xx, x, lambda = 0.3, method = "simple", proper = TRUE,nn="jf96");lines(xx, temp, lwd = 2, col = "red")
lines(qbckden(temp, x, lambda = 0.3, method = "simple", proper = TRUE,nn="jf96"),temp,lty=2,col="blue")

n=500
x = rgamma(n, shape = 2, scale = 2)
hist(rbckden(1000, x, lambda = 0.3, method = "simple"),100)

#     # correction only applies to point near boundary, so integrate density near there
#     # for those more than 10*lambda from boundary use standard kde
#     if (x > 10*lambda) {
#       p = mean(pnorm(x, kerncentres, lambda))
#     } else {
#       # Re-use density function to do numerical integration
#       bckdenint = try(integrate(dbckden, lower = 0, upper = x, kerncentres = kerncentres,
#         lambda = lambda, method = "simple", proper = proper, nn = nn,
#         subdivisions = 10000, rel.tol = 1.e-9, stop.on.error = FALSE))
# 
#       if (inherits(bckdenint, "try-error")) {
#         bckdenint$value = NA
#         warning("failed to renormalise boundary corrected KDE to integrate to unity")
#       }
#       p = bckdenint$value
#     }
#     p

# distance of kerncentres to lower boundary (p in paper)
    truncpoint = x/lambda

    # Use notation of Jones (1993) to make easier to check
    a0 = pnorm(truncpoint)
    a1 = -dnorm(truncpoint)  
    a2 = a0 + truncpoint * a1

    # weights in local linear fitting
    denom = (a2*a0 - a1^2)
    lx = a2/denom
    mx = a1/denom
    
    u = (x - kerncentres)/lambda

    mean(lx*(pnorm(u) - pnorm(0, kerncentres, lambda)) + mx*(dnorm(0, kerncentres, lambda) - dnorm(u)))













evmix.diag <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "weibullgpd", "gammagpd", "betagpd", "lognormgpd", "gng")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
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
  rlplot(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...)
  qplot(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...)
  pplot(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...)
  densplot(modelfit, upperfocus = TRUE, ...)
}
  
  

rlplot <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "weibullgpd", "gammagpd", "betagpd", "lognormgpd", "gng")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
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
    normgpd = qnormgpd(the.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpd = qweibullgpd(the.prob, modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpd = qnormgpd(the.prob, modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = qnormgpd(the.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = qnormgpd(the.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gng = qgng(the.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  

  # makes following code easier to read
  if (modelname == "gng") {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  # If focus on upper tail then must display at least upper 10%
  # even if tail fraction is smaller, otherwise may not look nice
  if (upperfocus) { 
    xlims = c(-1/log(1 - ifelse(phiu < 0.1, 0.1, phiu)), 10^max.emp.power)
    ylims = c(min(quantile(modelfit$x, 0.9, na.rm = TRUE), u), 
      max(c(modelfit$x, the.quant), na.rm = TRUE))
  } else {
    xlims = c(min(trans.emp.prob), 10^max.emp.power)
    ylims = range(c(modelfit$x, the.quant), na.rm = T)
  }
  if (ci) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.quantiles <- function(i, modelfit, modelname){
      simdata = switch(modelname,
        gpd = rgpd(length(modelfit$x), modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(length(modelfit$x), modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpd = rweibullgpd(length(modelfit$x), modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpd = rnormgpd(length(modelfit$x), modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = rnormgpd(length(modelfit$x), modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = rnormgpd(length(modelfit$x), modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gng = rgng(length(modelfit$x), modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  
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

  if (modelname == "gng") {
    abline(h = c(modelfit$ul, modelfit$ur), lty = 3)
    mtext(paste("ul =", formatC(modelfit$ul, digits = 3, format = "g")), side = 2, line = 2, at = modelfit$ul)
    mtext(paste("ur =", formatC(modelfit$ur, digits = 3, format = "g")), side = 2, line = 2, at = modelfit$ur)
  } else {
    abline(h = u, lty = 3)
    mtext(paste("u =", formatC(u, digits = 3, format = "g")), side = 2, line = 2, at = u)    
  }

  abline(v = -1/log(1 - phiu), lty = 3)
  mtext(paste("1/phiu =", formatC(1/phiu, format="d")), side = 1, line = 2, at = 1/phiu, las = 1)
  
  if (ci) {
    lines(trans.emp.prob, ci.q[1, ], lty=2)
    lines(trans.emp.prob, ci.q[2, ], lty=2)
    legend("bottomright", c("Data",paste("Fitted",modelname,"Model"),
      paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
      ifelse((modelname == "gng") & !upperfocus,
        "Tail Fractions and Thresholds","Upper Tail Fraction and Threshold")),
      pch = c(4, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white", col = c("black", "red", "black"))
  } else {
    legend("bottomright", c("Data",paste("Fitted",modelname,"Model"),
      ifelse((modelname == "gng") & !upperfocus,
        "Tail Fractions and Thresholds","Upper Tail Fraction and Threshold")),
      pch = c(4, rep(-1, 2)), lty = c(0, 1, 3), bg = "white", col = c("black", "red", "black"))
  }
}

qplot <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "weibullgpd", "gammagpd", "betagpd", "lognormgpd", "gng")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
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
  
  n = length(modelfit$x)
  emp.prob=ppoints(n)

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  if (modelname == "gpd"){
    upperfocus = TRUE  
  }
  
  the.quant = switch(modelname,
    gpd = qgpd(emp.prob, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = qnormgpd(emp.prob, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpd = qweibullgpd(emp.prob, modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpd = qnormgpd(emp.prob, modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = qnormgpd(emp.prob, modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = qnormgpd(emp.prob, modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gng = qgng(emp.prob, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  

  # makes following code easier to read
  if (modelname == "gng") {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  axislims = range(c(modelfit$x, the.quant), na.rm = T)
  
  if (ci) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.quantiles <- function(i, fit, modelname){
      simdata = switch(modelname,
        gpd = rgpd(length(modelfit$x), modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(length(modelfit$x), modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpd = rweibullgpd(length(modelfit$x), modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpd = rnormgpd(length(modelfit$x), modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = rnormgpd(length(modelfit$x), modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = rnormgpd(length(modelfit$x), modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gng = rgng(length(modelfit$x), modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  
      sort(simdata, na.last = FALSE)
    }
    sim.q = sapply(1:N, FUN = simulate.new.quantiles, fit = modelfit, modelname = modelname)
    ci.q = apply(sim.q, FUN=quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
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

  if (modelname == "gng") {
    abline(h = modelfit$ul, lty = 3)
    mtext(paste("ul =", formatC(modelfit$ul, digits = 3, format = "g")), side = 2, line = 2, at = modelfit$ul)
    abline(h = modelfit$ur, lty = 3)
    mtext(paste("ur =", formatC(modelfit$ur, digits = 3, format = "g")), side = 2, line = 2, at = modelfit$ur)
  } else {
    abline(h = modelfit$u, lty = 3)
    mtext(paste("u =", formatC(modelfit$u, digits = 3, format = "g")), side = 2, line = 2, at = modelfit$u)
  }
  
  if (ci) {
    lines(the.quant, ci.q[1, ], lty=2)
    lines(the.quant, ci.q[2, ], lty=2)
    legend("bottomright", c("Data",paste("Fitted", modelname, "Model"),
      paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
      ifelse((modelname == "gng") & !upperfocus, "Thresholds","Threshold")),
      pch = c(1, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white")
  } else {
    legend("bottomright", c("Data","Fitted Model",
      ifelse((modelname == "gng") & !upperfocus, "Thresholds","Threshold")),
      pch = c(1, rep(-1, 2)), lty = c(0, 1, 3), bg = "white")
  }
}

pplot <- function(modelfit, upperfocus = TRUE, ci = TRUE, alpha = 0.05, N = 1000, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call),"\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "weibullgpd", "gammagpd", "betagpd", "lognormgpd", "gng")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
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
  
  n = length(modelfit$x)
  emp.prob = ppoints(n)
  
  the.prob = switch(modelname,
    gpd = pgpd(modelfit$x, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = pnormgpd(modelfit$x, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpd = pweibullgpd(modelfit$x, modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpd = pnormgpd(modelfit$x, modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = pnormgpd(modelfit$x, modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = pnormgpd(modelfit$x, modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gng = pgng(modelfit$x, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  

  # makes following code easier to read
  if (modelname == "gng") {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  # If focus on upper tail then must display at least upper 10%
  # even if tail fraction is smaller, otherwise may not look nice
  if (upperfocus) { 
    axislims = c(1 - ifelse(phiu < 0.1, 0.1, phiu), 1)
  } else {
    axislims = c(0, 1)
  }

  if (ci) {
    # obtain CI by Monte Carlo simulation (ignores estimation uncertainty)
    simulate.new.probabilities <- function(i, modelfit, modelname){
      simdata = switch(modelname,
        gpd = rgpd(length(modelfit$x), modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = rnormgpd(length(modelfit$x), modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpd = rweibullgpd(length(modelfit$x), modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpd = rnormgpd(length(modelfit$x), modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = rnormgpd(length(modelfit$x), modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = rnormgpd(length(modelfit$x), modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gng = rgng(length(modelfit$x), modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  
      simdata = sort(simdata, na.last = FALSE)
      simprob = switch(modelname,
        gpd = pgpd(simdata, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        normgpd = pnormgpd(simdata, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        weibullgpd = pweibullgpd(simdata, modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gammagpd = pnormgpd(simdata, modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        betagpd = pnormgpd(simdata, modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        lognormgpd = pnormgpd(simdata, modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
        gng = pgng(simdata, modelfit$nmean, modelfit$nsd, 
          modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
          modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))
      simprob      
    }
    sim.p = sapply(1:N, FUN = simulate.new.probabilities, modelfit = modelfit, modelname = modelname)
    ci.p = apply(sim.p, FUN = quantile, MARGIN = 1, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  }
  
  plot(the.prob, emp.prob, pch = 'x', cex = 0.8,
    xlab = "Theoretical Probability", ylab = "Empirical Probability",
    xlim = axislims, ylim = axislims, ...)
  abline(c(0,1))

  if (modelname == "gng") {
    abline(h = modelfit$phiul, lty = 3)
    mtext(paste("phiul =", formatC(modelfit$phiul, digits = 3, format = "g")), side = 2, line = 2, at = modelfit$phiul)
    abline(h = 1 - modelfit$phiur, lty = 3)
    mtext(paste("phiur =", formatC(modelfit$phiur, digits = 3, format = "g")), side = 2, line = 2, at = 1 - modelfit$phiur)
  } else {
    abline(h = 1 - modelfit$phiu, lty = 3)
    mtext(paste("phiu =", formatC(modelfit$phiu, digits = 3, format = "g")), side = 2, line = 2, at = 1 - modelfit$phiu)
  }
  
  if (ci) {
    lines(the.prob, ci.p[1, ], lty=2)
    lines(the.prob, ci.p[2, ], lty=2)
    legend("bottomright", c("Data",paste("Fitted", modelname, "Model"),
      paste("Simulated Pointwise ", formatC(100 * (1 - alpha), format = "f", digits=1), "% CI", sep=""),
      ifelse((modelname == "gng") & !upperfocus, "Tail Fractions","Tail Fraction")),
      pch = c(1, rep(-1, 3)), lty = c(0, 1, 2, 3), bg = "white")
  } else {
    legend("bottomright", c("Data","Fitted Model",
      ifelse((modelname == "gng") & !upperfocus, "Tail Fractions","Tail Fraction")),
      pch = c(1, rep(-1, 2)), lty = c(0, 1, 3), bg = "white")
  }
}


densplot <- function(modelfit, upperfocus = TRUE, ...) {
  
  # Check properties of inputs
  if (missing(modelfit))
    stop("modelfit from extreme value mixture model must be given")
   
  if (!exists("call", modelfit))
    stop("modelfit from extreme value mixture model must be given")

  modelname = substring(strsplit(deparse(modelfit$call), "\\(")[[c(1, 1)]], 2)

  allmodelnames = c("gpd", "normgpd", "weibullgpd", "gammagpd", "betagpd", "lognormgpd", "gng")

  if (!(modelname %in% allmodelnames))
    stop("invalid extreme value mixture model given")

  if (!is.logical(upperfocus))
    stop("upperfocus must be logical")
  
  if (is.unsorted(modelfit$x, na.rm = TRUE))
    modelfit$x=sort(modelfit$x)

  # makes following code easier to read
  if (modelname == "gng") {
    phiu = modelfit$phiur
    u = modelfit$ur
  } else {
    phiu = modelfit$phiu
    u = modelfit$u    
  }

  xx = seq(min(modelfit$x, na.rm = T), max(modelfit$x, na.rm = TRUE), sd(modelfit$x, na.rm = T)/100)

  the.dens = switch(modelname,
    gpd = dgpd(xx, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    normgpd = dnormgpd(xx, modelfit$nmean, modelfit$nsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    weibullgpd = dweibullgpd(xx, modelfit$wshape, modelfit$wscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gammagpd = dnormgpd(xx, modelfit$gshape, modelfit$gscale, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    betagpd = dnormgpd(xx, modelfit$bshape1, modelfit$bshape2, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    lognormgpd = dnormgpd(xx, modelfit$lnmean, modelfit$lnsd, modelfit$u, modelfit$sigmau, modelfit$xi, modelfit$phiu),
    gng = dgng(xx, modelfit$nmean, modelfit$nsd, 
      modelfit$ul, modelfit$sigmaul, modelfit$xil, modelfit$phiul, 
      modelfit$ur, modelfit$sigmaur, modelfit$xir, modelfit$phiur))  

  kde.est = density(modelfit$x, from = min(modelfit$x, na.rm = T),
    to = max(modelfit$x, na.rm = T), na.rm = TRUE)
  
  # In GPD case NA are treated as below threshold and retained from fitting
  # so need to scale KDE
  if ((modelname == "gpd") & any(is.na(modelfit$x))) {
    kde.est$y = kde.est$y * modelfit$phiu
  }
  
  if (upperfocus) {
    xlims = c(u, max(modelfit$x, na.rm = TRUE))    
  } else {
    xlims = range(modelfit$x, na.rm = TRUE)
  }
  ylims = c(0, 1.1 * max(c(the.dens, kde.est$y)))    
  
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

  if (modelname == "gng") {
    abline(v = c(modelfit$ul, modelfit$ur), lty = 3)
  } else {
    abline(v = modelfit$u, lty = 3)
  }
  legend("topright", c("Fitted Density", 
    ifelse(modelname=="gng", "Thresholds", "Threshold"), "KDE"),
    lty = c(1, 3, 2), lwd = c(2, 1, 1.5), bg = "white")
}



x = sort(rnorm(1000))
#x=ifelse(x<=0.5,NA,x)
#x = ifelse(runif(1000)<0.01, NA,sort(rnorm(1000)))
#fit = fnormgpd(x, phiu = TRUE)
fit = fgng(x, phiul = TRUE, phiur = TRUE, std.err=FALSE)
#fit = fgpd(x, u = 0.5,std.err=FALSE)
#x11();
par(mfrow=c(2,2))
rlplot(fit,upperfocus=FALSE,ci=TRUE)
qplot(fit,upperfocus=FALSE,ci=TRUE)
pplot(fit,upperfocus=FALSE,ci=TRUE)
densplot(fit, upperfocus=TRUE)




# 
#   if (modelname %in% c("gammagpd", "weibullgpd", "betagpd")){
#     if (any(modelfit$x == 0, na.rm = TRUE)) {
#       y = log(modelfit$x + sd(modelfit$x)/1000)
#     } else {
#       y = log(modelfit$x)
#     }
#     kde.est = density(y, na.rm = TRUE)
#     kde.est$x = exp(kde.est$x)
#     kde.est$y = kde.est$y/kde.est$x
#   } else if ((modelname == "gpd") & (any(is.na(modelfit$x)))) {
#     xu = modelfit$x[which(modelfit$x > modelfit$u)] - modelfit$u
#     y = log(xu)
#     kde.est = density(y)
#     kde.est$x = exp(kde.est$x) + modelfit$u
#     kde.est$y = kde.est$y/kde.est$x
#   } else {
#     kde.est = 
#   }
#   
