#' estimate covariance function of a random field
#' assuming stationary (and optionally - isotropic)
#'
#'   Naive moment estimator
#' 
#'     Originally written for Max Chatfield
#'     original: Adrian Baddeley 15-19 may 2020
#'     $Revision: 1.11 $ $Date: 2021/05/03 02:40:27 $

spatcov <- function(X, Y=X, ...,
                    correlation=FALSE, isotropic=TRUE,
                    clip=TRUE, pooling=TRUE) {
  stopifnot(is.im(X))
  eX <- X - mean(X)
  if(correlation) eX <- eX/sqrt(mean(eX^2))
  if(missing(Y) || is.null(Y)) {
    #' spatial covariance of X
    A <- imcov(eX)
  } else {
    #' spatial cross-covariance of X and Y
    stopifnot(is.im(Y))
    eY <- Y - mean(Y)
    if(correlation) eY <- eY/sqrt(mean(eY^2))
    A <- imcov(eX, eY)
  }
  B <- setcov(Window(X))
  if(!(isotropic && pooling)) {
    #' first estimate covariance as function of vector argument
    Z <- A/B
    #' deal with numerical errors at extremes
    pixelarea <- with(X, xstep * ystep)
    Z[B < pixelarea] <- 0
  }
  if(isotropic) {
    #' result is a function of lag distance
    if(pooling) {
      mA <- rotmean(A)
      mB <- rotmean(B)
      f <- eval.fv(mA/mB)
    } else {
      f <- rotmean(Z)
    }
    #' give it more meaningful labels
    f <- rebadge.fv(f,
                    new.ylab=quote(C(r)),
                    new.fname="C",
                    tags=fvnames(f, ".y"),
                    new.tags="est",
                    new.desc="estimate of %s",
                    new.labl="hat(%s)(r)")
    if(clip) 
      attr(f, "alim") <- c(0, shortside(Frame(X))/2)
    result <- f
  } else {
    #' result is an image representing a function of lag vector
    Z <- A/B
    #' return an image representing a function of lag vector
    if(clip) {
      Box <- Frame(Z)
      b <- sidelengths(Box)
      Bclip <- trim.rectangle(Box, b[1]/4, b[2]/4)
      Z <- Z[Bclip, drop=FALSE, tight=TRUE]
    }
    result <- Z
  } 
  return(result)
}

pairMean <- function(fun, W, V=NULL, ..., normalise=TRUE) {
  #' fun is a function of pairwise distance
  if(!is.function(fun))
    stop("fun should be a function in the R language")
  #' W is the domain over which to integrate
  W <- as.owin(W)
  FD <- distcdf(W, V, ..., savedenom=!normalise)
  result <- as.numeric(stieltjes(fun, FD, ...))
  if(!normalise) result <- result * attr(FD, "denom")
  return(result)
}

  
