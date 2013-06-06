#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# inverse cumulative distribution function for dynamically weighted mixture model
qdwm = function(p, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, lower.tail = TRUE, qinit = NULL) {
  
  # Check properties of inputs
  if (missing(p))
    stop("p must be a non-empty numeric vector")
  
  if (length(p) == 0 | mode(p) != "numeric") 
    stop("p must be a non-empty numeric vector")
  
  if (min(p, na.rm = TRUE) < 0 | max(p, na.rm = TRUE) > 1)
    stop("probability must be between 0 and 1 (inclusive)")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(p), length(wshape), length(wscale), 
              length(cmu), length(ctau), length(sigmau), length(xi),
              ifelse(length(qinit) == 0, 1, length(qinit)))
  n = max(linputs)
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector, with vectors all same length")
  
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
  
  if (!is.logical(lower.tail))
    stop("lower.tail must be logical")
  
  if (length(lower.tail) != 1)
    stop("lower.tail must be of length 1")
  
  if (!lower.tail) p = 1 - p
  
  if (!is.null(qinit)) {
    if (mode(qinit) != "numeric")
      stop("qinit must be numeric")
    if (any(!is.finite(qinit)))
      stop("qinit must be numeric")
    if (any(qinit < 0))
      stop("qinit must be non-negative")    
    qinit = rep(qinit, length.out = n)
  } else {
    qinit = rep(NA, length.out = n)    
  }
  
  p = rep(p, length.out = n)
  wshape = rep(wshape, length.out = n)
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  pdmmmin = function(q, cprob, wshape, wscale, cmu, ctau, sigmau, xi) {
     cdfmm = pdwm(q, wshape, wscale, cmu, ctau, sigmau, xi)
     if (is.na(cdfmm)) {
        qdiff = 1e6
     } else {
        qdiff = abs(cdfmm - cprob)
     }
    qdiff
  }
  
  findqdmm = function(cprob, wshape, wscale, cmu, ctau, sigmau, xi, qinit) {
    if (is.na(qinit)) {
      qwbl = qweibull(cprob, shape = wshape, scale = wscale)
      qgp = qgpd(cprob, sigmau = sigmau, xi = xi)
      qinit = mean(c(qwbl, qgp))
    }
    
    gt = try(nlm(pdmmmin, qinit, cprob, wshape, wscale, cmu, ctau, sigmau, xi,
      gradtol = 1e-10, steptol = 1e-10)$estimate)

    if (inherits(gt, "try-error")) {
      gt = try(nlm(pdmmmin, qgpd(cprob, sigmau = sigmau, xi = xi),
        cprob, wshape, wscale, cmu, ctau, sigmau, xi,
        gradtol = 1e-10, steptol = 1e-10)$estimate)
      
      if (inherits(gt, "try-error")) {
        gt = NA
      }
    }
    return(gt)
  }
 
  q = rep(NA, length.out = n)
  for (i in 1:n) {
    q[i] = findqdmm(p[i], wshape[i], wscale[i], cmu[i], ctau[i], sigmau[i], xi[i], qinit[i])
  }     

  q                     
}

