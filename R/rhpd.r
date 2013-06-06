#' @export
#' @aliases  hpd dhpd phpd qhpd rhpd
#' @rdname hpd

# random number generation for hybrid Pareto model
rhpd <- function(n = 1, nmean = 0, nsd = 1, xi = 0) {
  
  # Check properties of inputs
  if (length(n) != 1 | mode(n) != "numeric") 
    stop("sample size must be positive integer")
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps) | n < 1)
    stop("sample size must be positive integer")
  
  # parameter inputs could be single values or vectors to allow for nonstationary modelling
  # all input vectors must be length n or scalar
  linputs = c(length(nmean), length(nsd), length(xi))
  
  if (sum(linputs[linputs != 1] != n) > 0)
    stop("Data and parameters must be either scalar or vector of length n")
  
  if (mode(nmean) != "numeric" | mode(nsd) != "numeric" | mode(xi) != "numeric")
    stop("parameters must be numeric")
  
  if (any(!is.finite(c(nmean, nsd, xi))))
    stop("parameters must be numeric")
  
  if (min(nsd) <= 0)
    stop("normal standard deviation must be non-negative")

  if (any(xi == 1))
    stop("shape cannot be 1")
    
  qhpd(runif(n), nmean, nsd, xi)
}

