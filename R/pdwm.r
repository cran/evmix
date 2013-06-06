#' @export
#' @aliases dwm ddwm pdwm qdwm rdwm
#' @rdname  dwm

# cumulative distribution function for dynamically weighted mixture model
pdwm = function(q, wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  if (missing(q))
    stop("q must be a non-empty numeric vector")
  
  if (length(q) == 0 | mode(q) != "numeric") 
    stop("q must be a non-empty numeric vector")
  
  if (any(is.infinite(q)))
    warning("infinite quantiles are set to NaN")
  
  q[is.infinite(q)]=NA # user will have to deal with infinite cases
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be same length or scalar
  linputs = c(length(q), length(wshape), length(wscale), 
            length(cmu), length(ctau), length(sigmau), length(xi))
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

  q = rep(q, length.out = n)
  wshape = rep(wshape, length.out = n)  
  wscale = rep(wscale, length.out = n)
  cmu = rep(cmu, length.out = n)
  ctau = rep(ctau, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
    
  rx <- function(x, wshape, wscale, cmu, ctau, sigmau, xi) {
    (dgpd(x, sigmau = sigmau, xi = xi) - dweibull(x, shape = wshape, scale = wscale))*atan((x - cmu)/ctau)
  }

  rxw <- function(x, wshape, wscale, cmu, ctau) {
     (1 - pcauchy(x, location = cmu, scale = ctau))*dweibull(x, shape = wshape, scale = wscale)
  }

  rxg <- function(x, cmu, ctau, sigmau, xi){
     pcauchy(x, location = cmu, scale = ctau)*dgpd(x, sigmau = sigmau, xi = xi)
  }

  p = q
  
  whichnonmiss = which(!is.na(q))

  z1 = z2 = z = rep(NA, length.out = n)
  for (i in 1:n) {
    if ((max(linputs[-1]) == 1) & (i == 1)) {
      r = try(integrate(rx, wshape = wshape[i], wscale = wscale[i],
        cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
        lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         
  
      if (inherits(r, "try-error")) {
        z[i] = NA
      } else {              
        z[i] = 1 + r/pi
      }
    } else if ((max(linputs[-1]) == 1) & (i > 1)) {
      z[i] = z[1]      
    } else {
      r = try(integrate(rx, wshape = wshape[i], wscale = wscale[i],
        cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
        lower = 0, upper = Inf, subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         
  
      if (inherits(r, "try-error")) {
        z[i] = NA
      } else {              
        z[i] = 1 + r/pi
      }
    }
    
    r1 = try(integrate(rxw, wshape = wshape[i], wscale = wscale[i], cmu = cmu[i], ctau = ctau[i],
      lower = 0, upper = q[i], subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         

    if (inherits(r1, "try-error")) {
      z1[i] = NA
    } else {              
      z1[i] = r1
    }

    r2 = try(integrate(rxg, cmu = cmu[i], ctau = ctau[i], sigmau = sigmau[i], xi = xi[i],
      lower = 0, upper = q[i], subdivisions = 10000, rel.tol = 1e-10, stop.on.error = FALSE)$value)         

    if (inherits(r2, "try-error")) {
      z2[i] = NA
    } else {              
      z2[i] = r2
    }
  }

  p[whichnonmiss] = (z1[whichnonmiss] + z2[whichnonmiss])/z[whichnonmiss]

  if (!lower.tail) p = 1 - p

  p
}

