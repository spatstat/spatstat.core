#'
#'    varcount.R
#'
#'   Variance of N(B)
#'
#'  $Revision: 1.14 $  $Date: 2021/08/12 06:35:51 $
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

  varcountEngine <- function(g, B, lambdaB, f=1, R=NULL,
                             what=c("variance","excess","pairs","squared")) {
    ## variance = var(N)
    ## excess   = var(N) - E(N)
    ## pairs    = E[ N(N-1) ]
    ## squared  = E[N^2]
    what <- match.arg(what)
    g1 <- function(r) { g(r) - 1 }
    dr <- if(is.null(R)) NULL else R/100
    if(missing(f) || identical(f, 1)) {
      v <- switch(what,
                  variance = integral(lambdaB) + dublin(g1, B, lambdaB, dr=dr),
                  excess   =                     dublin(g1, B, lambdaB, dr=dr),
                  pairs    =                     dublin(g,  B, lambdaB, dr=dr),
                  squared  = integral(lambdaB) + dublin(g,  B, lambdaB, dr=dr))
    } else if(min(f) >= 0) {
      ## nonnegative integrand
      v <- switch(what,
        variance = integral(lambdaB * f^2) + dublin(g1, B, lambdaB * f, dr=dr),
        excess   =                           dublin(g1, B, lambdaB * f, dr=dr),
        pairs    =                           dublin(g,  B, lambdaB * f, dr=dr),
        squared  = integral(lambdaB * f^2) + dublin(g,  B, lambdaB * f, dr=dr))
    } else if(max(f) <= 0) {
      ## nonpositive integrand
      v <- switch(what,
        variance=integral(lambdaB * f^2) + dublin(g1, B, lambdaB * (-f), dr=dr),
        excess  =                          dublin(g1, B, lambdaB * (-f), dr=dr),
        pairs   =                          dublin(g,  B, lambdaB * (-f), dr=dr),
        squared =integral(lambdaB * f^2) + dublin(g,  B, lambdaB * (-f), dr=dr))
    } else {
      ## integrand has both positive and negative parts
      lamfplus <- eval.im(lambdaB * pmax(0, f))
      lamfminus <- eval.im(lambdaB * pmax(0, -f))
      h <- switch(what,
                  variance = g1,
                  excess   = g1,
                  pairs    = g,
                  squared  = g)
      co <- (dublin(h, B, lamfplus, dr=dr) 
            + dublin(h, B, lamfminus, dr=dr)
            - dublin(h, B, lamfplus, lamfminus, dr=dr)
            - dublin(h, B, lamfminus, lamfplus, dr=dr))
      v <- switch(what,
                  variance = integral(lambdaB * f^2) + co,
                  excess   =                           co,
                  pairs    =                           co,
                  squared  = integral(lambdaB * f^2) + co)
    }
    return(v)
  }

  dublin <- function(h, B, f, f2, dr=NULL) {
    ## Double integral
    ##          \int_B \int_B h(|u-v|) f(u) f(v) du dv    
    ## or       \int_B \int_B h(|u-v|) f(u) f2(v) du dv
    ## Assume h, f, f2 are nonnegative
    ## dr = maximum spacing of argument
    if(missing(f2)) {
      ## \int_B \int_B h(|u-v|) f(u) f(v) du dv
      M <- distcdf(B, dW=f, nr=NULL, delta=dr)
      ## integrate 
      a <- integral(f)^2 * as.numeric(stieltjes(h, M))
    } else {
      ## \int_B \int_B h(|u-v|) f(u) f2(v) du dv
      M <- distcdf(B, dW=f, dV=f2, nr=NULL, delta=dr)
      ## integrate 
      a <- integral(f) * integral(f2) * as.numeric(stieltjes(h, M))
    }
    return(a)
  }


  varcountEngine
})


