#'
#'   Selection of threshold 
#'
#'   Copyright (c) 2020 Adrian Baddeley, Warick Brown, Robin K. Milne,
#'   Gopalan Nair, Suman Rakshit, Tom Lawrence, Aloke Phatak, Shih Ching Fu
#'
#'   GNU Public Licence >= 2
#'
#'   $Revision: 1.4 $ $Date: 2022/05/20 07:39:55 $
#'
#' 
#' threshold selection
#' inputs:
#'     X   deposit locations
#'     Z   covariate

thresholdSelect <- function(X, Z, method=c("Y", "LL", "AR", "t", "C"), Zname) {
  if(!is.ppp(X)) stop("X should be a point pattern (class ppp)")
  if(missing(Zname)) Zname <- short.deparse(substitute(Z))
  method <- match.arg(method)
  a <- spatialCovariateEvidence(X, Z, jitter=FALSE)$values
  FF <- ecdf(a$ZX)
  GG <- ecdf(a$Zvalues)
  n <- npoints(X)
  A <- area(Window(a$Zimage))
  zrange <- range(range(a$ZX), range(a$Zimage))
  zz <- seq(zrange[1], zrange[2], length.out=1028)
  nz <- n * (pz <- FF(zz))
  Az <- A * (sz <- GG(zz))
  Cz <- log((nz/Az)/((n-nz)/(A-Az)))
  yy <- switch(method,
               C = Cz,
               t = Cz/sqrt(1/nz + 1/(n-nz)),
               LL = {  n * log(nz/Az) - (n-nz) * Cz - n },
               AR = { sqrt(sz * (1-sz)) * (nz/sz - (n-nz)/(1-sz)) },
               Y = { pz - sz })
  yy[!is.finite(yy)] <- -Inf
  critname <- switch(method,
                     C = "WofE contrast",
                     t = "studentised contrast",
                     LL = "profile log likelihood",
                     AR = "Akman-Raftery criterion",
                     Y = "Youden criterion")
  bw.optim(yy, zz, optimum="max",
           cvname=method, hname=Zname, 
           criterion=critname,
           unitname=if(inherits(Z, "distfun")) unitname(X) else NULL)
}

#' confidence interval for threshold

thresholdCI <- local({
  
  thresholdCI <- function(X, Z, confidence=0.95, nsim=1000, parametric=FALSE) {
    #' bootstrap confidence interval for Youden estimate only.
    if(!is.ppp(X)) stop("X should be a point pattern (class ppp)")
    a <- spatialCovariateEvidence(X, Z, jitter=FALSE)$values
    FF <- ecdf(a$ZX)
    GG <- ecdf(a$Zvalues)
    est <- Youden(FF,GG)
    b <- simthresh(FF, GG, npoints(X), nsim, parametric)
    zCI <- quantCI(b$z, est[["z"]], confidence=confidence)
    sCI <- quantCI(b$s, est[["s"]], confidence=confidence)
    rbind(z=zCI, s=sCI)
  }

  #' Underlying code based on cumulative distribution functions
  #' inputs:
  #'    F = ecdf of covariate values for data points
  #'    G = ecdf of covariate values for study region

  Youden <- function(F, G) {
    zz <- get("x", envir=environment(F))
    iopt <- which.max(F(zz) - G(zz))
    zopt <- zz[iopt]
    sopt <- G(zopt)
    return(c(z=zopt, s=sopt))
  }

  Fpredicted <- function(F, G, zest) {
    if(missing(zest)) zest <- Youden(F,G)[["z"]]
    plow <- F(zest)
    glow <- G(zest)
    #' mixture of unif[0, glow] and unif[glow, 1] with weights plow, 1-plow
    zz <- get("x", envir=environment(G))
    pp <- get("y", envir=environment(G))
    qq <- ifelse(pp < glow, plow*(pp/glow), plow + (1-plow)*(pp-glow)/(1-glow))
    FF <- approxfun(zz, qq, rule=2)
    return(FF)
  }

  inversefunction <- function(F) {
    zz <- get("x", envir=environment(F))
    pz <- get("y", envir=environment(F))
    Finv <- approxfun(pz, zz, rule=2)
    return(Finv)
  }

  simthresh <- function(F, G, ndata, nsim=100, parametric) {
    check.1.integer(nsim)
    stopifnot(nsim > 1)
    if(parametric) F <- Fpredicted(F, G)
    Finv <- inversefunction(F)
    zout <- sout <- numeric(nsim)
    zz <- get("x", envir=environment(G))
    for(isim in 1:nsim) {
      zsim <- Finv(runif(ndata))
      Fhat <- ecdf(zsim)
      iopt <- which.max(Fhat(zz) - G(zz))
      zopt <- zz[iopt]
      sopt <- G(zopt)
      zout[isim] <- zopt
      sout[isim] <- sopt
    }
    return(data.frame(z=zout, s=sout))
  }

  quantCI <- function(x, xest, confidence=0.95) {
    xleft <- quantile(x[x<=xest], 1-confidence)
    xright <- quantile(x[x>=xest], confidence)
    achieved <- mean(x >= xleft & x <= xright)
    return(c(lo=unname(xleft), hi=unname(xright), conf=achieved))
  }

  thresholdCI
})


