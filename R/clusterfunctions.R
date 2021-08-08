## clusterfunctions.R
##
## Contains the generic functions:
##  - clusterkernel
##  - clusterfield
##  - clusterradius.
##
##   $Revision: 1.7 $  $Date: 2021/08/08 08:54:55 $
##

clusterkernel <- function(model, ...) {
  UseMethod("clusterkernel")
}

clusterkernel.kppm <- function(model, ...) {
  kernelR <- Kpcf.kppm(model, what = "kernel")
  f <- function(x, y = 0, ...){
    kernelR(sqrt(x^2+y^2))
  }
  return(f)
}

clusterkernel.character <- function(model, ...){
  info <- spatstatClusterModelInfo(model, onlyPCP = TRUE)
  internalkernel <- info$kernel
  dots <- list(...)
  par <- c(kappa = 1, scale = dots$scale)
  par <- info$checkpar(par, old = TRUE)
  nam <- info$clustargsnames
  margs <- NULL
  if(!is.null(nam))
    margs <- dots[nam]
  f <- function(x, y = 0, ...){
    internalkernel(par = par, rvals = sqrt(x^2+y^2), margs = margs)
  }
  return(f)
}

clusterfield <- function(model, locations = NULL, ...) {
    UseMethod("clusterfield")
}

clusterfield.kppm <- function(model, locations = NULL, ...) {
    f <- clusterkernel(model)
    if(is.null(locations)){
        if(!is.stationary(model))
            stop("The model is non-stationary. The argument ",
                 sQuote("locations"), " must be given.")
        locations <- centroid.owin(Window(model), as.ppp = TRUE)
    }
    clusterfield.function(f, locations, ..., mu = model$mu)
}

clusterfield.character <- function(model, locations = NULL, ...){
    f <- clusterkernel(model, ...)
    clusterfield.function(f, locations, ...)
}

clusterfield.function <- function(model, locations = NULL, ..., mu = NULL) {
  if(is.null(locations)){
    locations <- ppp(.5, .5, window=square(1))
  } else if(!is.ppp(locations))
    stop("Argument ", sQuote("locations"), " must be a point pattern (ppp).")

  if("sigma" %in% names(list(...)) && "sigma" %in% names(formals(model)))
    warning("Currently ", sQuote("sigma"),
            "cannot be passed as an extra argument to the kernel function. ",
            "Please redefine the kernel function to use another argument name.")

  rslt <- density(locations, kernel=model, ..., edge=FALSE)
  if(is.null(mu))
    return(rslt)
  
  mu <- as.im(mu, W=rslt)
  if(min(mu)<0)
    stop("Cluster reference intensity ", sQuote("mu"), " is negative.")

  return(rslt*mu)
}

clusterradius <- function(model, ...){
    UseMethod("clusterradius")
}

clusterradius.character <- function(model, ..., thresh = NULL, precision = FALSE){
  info <- spatstatClusterModelInfo(model, onlyPCP=FALSE)
  if(!isTRUE(info$isPCP)) {
    warning("cluster radius is only defined for cluster processes", call.=FALSE)
    return(NA)
  }
  rmax <- info$range(..., thresh = thresh)
  if(precision && is.function(info$ddist)){
    ddist <- function(r) info$ddist(r, ...)
    prec <- integrate(ddist, 0, rmax)
    attr(rmax, "prec") <- prec
  }
  return(rmax)
}

clusterradius.kppm <- function(model, ..., thresh = NULL, precision = FALSE){
    a <- list(model = model$clusters,
              thresh = thresh,
              precision = precision)
    a <- append(a, as.list(c(model$clustpar, model$clustargs)))
    do.call(clusterradius.character, a)
}

reach.kppm <- function(x, ..., epsilon) {
  thresh <- if(missing(epsilon)) NULL else epsilon
  if(x$isPCP) 
    return(2 * clusterradius.kppm(x, ..., thresh=thresh))
  ## use pair correlation
  g <- pcfmodel(x)
  ## find upper bound
  if(is.null(thresh)) thresh <- 0.01
  f <- function(r) { g(r) - 1 - thresh }
  scal <- as.list(x$par)$scale %orifnull% 1
  for(a in scal * 2^(0:10)) { if(f(a) < 0) break; }
  if(f(a) > 0) return(Inf)
  ## solve g(r) = 1 + epsilon
  b <- uniroot(f, c(0, a))$root
  return(b)
}
