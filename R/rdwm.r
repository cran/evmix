#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# random number generation for dynamically weighted mixture model
rdwm = function(n = 1, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2), xi = 0) {
  
  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(q), length(wshape), length(wscale), 
              length(cmu), length(ctau), length(sigmau), length(xi))
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")
  
  if (mode(wshape) != "numeric" | mode(wscale) != "numeric" | mode(cmu) != "numeric" |
    mode(ctau) != "numeric" | mode(sigmau) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
 
  if (any(!is.finite(c(wshape, wscale, cmu, ctau, sigmau, xi))))
    stop("parameters must be numeric")
  
  if (min(wshape) < 0)
    stop("weibull shape must be non-negative")
  
  if (min(wscale) < 0)
    stop("weibull scale must be non-negative")
  
  if (min(ctau) < 0)
    stop("Cauchy scale must be non-negative")
  
  if (min(sigmau) <= 0)
    stop("scale must be non-negative")  
  
  if (any(xi == 1))
    stop("shape cannot be 1")
  
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # Simulation scheme proposed by authors
  r = rep(NA, length.out = n)
  i = 1
  while (i <= n) {
    u = runif(1)
    if (u < 0.5) {
      rw = rweibull(1, shape = wshape[i], scale = wscale[i])
      pw = pcauchy(rw, location = cmu[i], scale = ctau[i])
      # accept or reject
      v = runif(1)
      if (v <= (1 - pw)) {
        r[i] = rw
        i = i + 1
      }
    } else {
      rg = rgpd(1, sigmau = sigmau[i], xi = xi[i])
      pg = pcauchy(rg, location = cmu[i], scale = ctau[i])

      # accept or reject
      v = runif(1)
      if (v <= pg){
        r[i] = rg
        i = i + 1
      }
    }
  } 
  
  r
}
