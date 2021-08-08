#'
#'    varcount.R
#'
#'   Variance of N(B)
#'
#'  $Revision: 1.11 $  $Date: 2021/08/08 08:56:54 $
#'

varcount <- function(model, B, ..., dimyx=NULL) {
  stopifnot(is.owin(B) || is.im(B) || is.function(B))
  g <- pcfmodel(model)
  if(!is.function(g))
    stop("Pair correlation function cannot be computed")
  R <- reach(model, epsilon=0.001) 
  if(!isTRUE(is.finite(R))) R <- NULL
  if(is.owin(B)) {
    lambdaB <- predict(model, locations=B, ngrid=dimyx, type="intensity")
    v <- varcountEngine(g, B, lambdaB, R=R)
  } else {
    f <- if(is.im(B)) B else as.im(B, W=as.owin(model), ..., dimyx=dimyx)
    B <- as.owin(f)
    lambdaB <- predict(model, locations=B, type="intensity")
    v <- varcountEngine(g, B, lambdaB, f, R=R)
  } 
  return(v)
}

varcountEngine <- local({

  varcountEngine <- function(g, B, lambdaB, f=1, R=NULL) {
    if(missing(f) || identical(f, 1)) {
      v <- integral(lambdaB) + covterm(g, B, lambdaB, R=R)
    } else if(min(f) >= 0) {
      ## nonnegative integrand
      v <- integral(lambdaB * f^2) + covterm(g, B, lambdaB * f, R=R)
    } else if(max(f) <= 0) {
      ## nonpositive integrand
      v <- integral(lambdaB * f^2) + covterm(g, B, lambdaB * (-f), R=R)
    } else {
      ## integrand has both positive and negative parts
      lamfplus <- eval.im(lambdaB * pmax(0, f))
      lamfminus <- eval.im(lambdaB * pmax(0, -f))
      v <- integral(lambdaB * f^2) +
        (covterm(g, B, lamfplus, R=R) +
         covterm(g, B, lamfminus, R=R)
          - covterm(g, B, lamfplus, lamfminus, R=R)
         - covterm(g, B, lamfminus, lamfplus, R=R))
    }
    return(v)
  }

  covterm <- function(g, B, f, f2, R=NULL) {
    ## Evaluate \int_B \int_B (g(u-v) - 1) f(u) f(v) du dv    
    ## or       \int_B \int_B (g(u-v) - 1) f(u) f2(v) du dv
    g1 <- function(r) { g(r) - 1 }
    if(missing(f2)) {
      ## \int_B \int_B (g(u-v) - 1) f(u) f(v) du dv
      H <- distcdf(B, dW=f, nr=NULL)
      ## ensure H is fine enough
      H <- refine(H, R)
      ## integrate 
      a <- integral(f)^2 * as.numeric(stieltjes(g1, H))
    } else {
      ## \int_B \int_B (g(u-v) - 1) f(u) f2(v) du dv
      H <- distcdf(B, dW=f, dV=f2, nr=NULL)
      ## ensure H is fine enough
      H <- refine(H, R)
      ## integrate 
      a <- integral(f) * integral(f2) * as.numeric(stieltjes(g1, H))
    }
    return(a)
  }


  refine <- function(H, R=NULL, nrequire=100, verbose=FALSE, force=FALSE) {
    ## Ensure H has enough quadrature points to compute an integral on [0,R]
    if(!is.null(R) && is.finite(R)) {
      rstep <- mean(diff(H$r))
      nuseful <- R/rstep
      if(verbose) {
        splat("R=", R, "rstep=", rstep)
        splat("H(R)=", as.function(H)(R))
      }
    } else {
      nuseful <- Inf
    }
    if(verbose) splat("nuseful = ", nuseful)
    if(nuseful < nrequire) {
      ## interpolate H
      if(verbose) {
        plot(H, xlim=c(0, R/2))
      }
      H <- interpCDF(H, n=ceiling(nrequire/nuseful))
      if(verbose) {
        plot(H, add=TRUE, xlim=c(0,R), col=2)
        splat("New rstep=", mean(diff(H$r)))
      }
      if(force) {
        ## force CDF to be nondecreasing and to start from 0
        Fr <- H[["f"]]
        Fr[1] <- 0
        Fr <- cummax(Fr)
        H[["f"]] <- Fr
      }
    }
    return(H)
  }
  
  interpCDF <- function(H, ..., method=c("smooth.spline", "loess"),
                        delta=NULL, n=NULL) {
    method <- match.arg(method)
    rname <- fvnames(H, ".x")
    rold <- H[[rname]]
    rpos <- (rold > 0)
    if(is.null(delta) == is.null(n))
      stop("Exactly one of the arguments 'delta' or 'n' should be given")
    if(!is.null(n)) {
      delta <- mean(diff(rold))/n
    } else {
      check.1.real(delta)
      stopifnot(delta > 0)
    }
    rnew <- seq(min(rold), max(rold), by=delta)
    ## initialise result
    newvalues <- vector(mode="list", length=ncol(H))
    names(newvalues) <- colnames(H)
    newvalues[[rname]] <- rnew
    ## process each column of function values
    nama <- fvnames(H, ".a")
    for(ynam in nama) {
      yy <- H[[ynam]]
      ok <- is.finite(yy) & rpos
      yok <- yy[ok]
      rok <- rold[ok]
      switch(method,
             smooth.spline = {
               ss <- smooth.spline(x=rok, y=yok/rok^2, ...)
               yhat <- predict(ss, rnew)$y * rnew^2
             },
             loess = {
               df <- data.frame(x=rok, y=yok/rok^2)
               lo <- loess(y ~ x, df, ...)
               yhat <- predict(lo, data.frame(x=rnew)) * rnew^2
             })
      newvalues[[ynam]] <- yhat
    }
    newH <- as.data.frame(newvalues)
    ## copy attributes
    anames <- setdiff(names(attributes(H)),
                      c("row.names", "dim", "dimnames", "names", "tsp"))
    for(e in anames)
      attr(newH, e) <- attr(H, e)
    return(newH)
  }

  varcountEngine
})


