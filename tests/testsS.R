#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.core
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.core)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#'   tests/sdr.R
#'
#'   $Revision: 1.2 $ $Date: 2020/05/01 09:59:59 $

if(FULLTEST) {
local({
  AN <- sdr(bei, bei.extra, method="NNIR")
  AV <- sdr(bei, bei.extra, method="SAVE")
  AI <- sdr(bei, bei.extra, method="SIR")
  AT <- sdr(bei, bei.extra, method="TSE")
  subspaceDistance(AN$B, AV$B)
  dimhat(AN$M)
})
}
##
##  tests/segments.R
##   Tests of psp class and related code
##                      [SEE ALSO: tests/xysegment.R]
##
##  $Revision: 1.30 $  $Date: 2020/11/09 09:01:13 $


local({
  if(ALWAYS) { # C code
    #' tests of density.psp
    Y <- edges(letterR)
    Window(Y) <- grow.rectangle(Frame(Y), 0.4)
    YC <- density(Y, 0.2, method="C", edge=FALSE, dimyx=64)
    YI <- density(Y, 0.2, method="interpreted", edge=FALSE, dimyx=64)
    YF <- density(Y, 0.2, method="FFT", edge=FALSE, dimyx=64)
    xCI <- max(abs(YC/YI - 1))
    xFI <- max(abs(YF/YI - 1))
    cat(paste("xCI =", xCI, "\txFI =", signif(xFI, 5)), fill=TRUE)
    if(xCI > 0.01) stop(paste("density.psp C algorithm relative error =", xCI))
    if(xFI > 0.1) stop(paste("density.psp FFT algorithm relative error =", xFI))

    B <- square(0.3)
    density(Y, 0.2, at=B)
    density(Y, 0.2, at=B, edge=TRUE, method="C")
    Z <- runifpoint(3, B)
    density(Y, 0.2, at=Z)
    density(Y, 0.2, at=Z, edge=TRUE, method="C")
  }
})



#
## tests/sigtraceprogress.R
#
## Tests of *.sigtrace and *.progress
#
## $Revision: 1.5 $ $Date: 2020/05/01 09:59:59 $

if(FULLTEST) {
local({
  plot(dclf.sigtrace(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dclf.progress(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.sigtrace(redwood, nsim=5, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.progress(redwood, nsim=5, alternative="greater", rmin=0.02,
                   verbose=FALSE))
  ## test 'leave-two-out' algorithm
  a <- dclf.sigtrace(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                     verbose=FALSE)
  aa <- dclf.progress(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                      verbose=FALSE)
  b <- dg.sigtrace(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2)
  bb <- dg.progress(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2,
                    verbose=FALSE)
  ## other code blocks
  e <- mad.progress(redwood, nsim=5)
  e <- mad.progress(redwood, nsim=19, alpha=0.05)
  f <- dclf.progress(redwood, nsim=5, scale=function(x) x^2)
  f <- dclf.progress(redwood, nsim=5, normalize=TRUE, deflate=TRUE)
  g <- dg.progress(redwood, nsim=5, scale=function(x) x^2)
  g <- dg.progress(redwood, nsim=5, normalize=TRUE, deflate=TRUE)
})
}
#
# tests/slrm.R
#
# $Revision: 1.3 $ $Date: 2020/05/01 09:59:59 $
#
# Test slrm fitting and prediction when there are NA's
#

if(ALWAYS) {
local({
  X <- copper$SouthPoints
  W <- owin(poly=list(x=c(0,35,35,1),y=c(1,1,150,150)))
  Y <- X[W]
  fit <- slrm(Y ~ x+y)
  pred <- predict(fit)
  extractAIC(fit)
  fitx <- update(fit, . ~ x)
  simulate(fitx, seed=42)
  if(FULLTEST) {
    unitname(fitx)
    unitname(fitx) <- "km"

    mur <- solapply(murchison,rescale, 1000, "km")
    mur$dfault <- distfun(mur$faults)
    fut <- slrm(gold ~ dfault, data=mur, splitby="greenstone")
    A <- model.images(fut)
  }
})
}


#'
#'    tests/ssf.R
#'
#'   Tests of 'ssf' class
#'
#'   $Revision: 1.3 $ $Date: 2020/05/01 09:59:59 $
#'

if(FULLTEST) {
local({
  Y <- cells[1:5]
  X <- rsyst(Window(Y), 5)
  Z <- runifpoint(3, Window(Y))
  f1 <- ssf(X, nncross(X,Y,what="dist"))
  f2 <- ssf(X, nncross(X,Y,what="dist", k=1:2))
  image(f1)
  g1 <- as.function(f1)
  g1(Z)
  g2 <- as.function(f2)
  g2(Z)
  plot(f1, style="contour")
  plot(f1, style="imagecontour")
  contour(f1)
  apply.ssf(f2, 1, sum)
  range(f1)
  min(f1)
  max(f1)
  integral(f1, weights=tile.areas(dirichlet(X)))
})
}
#
#   tests/step.R
#
#   $Revision: 1.5 $  $Date: 2020/05/01 09:59:59 $
#
# test for step() operation
#
if(FULLTEST) {
local({
  Z <- as.im(function(x,y){ x^3 - y^2 }, nztrees$window)
  fitP <- ppm(nztrees ~x+y+Z, covariates=list(Z=Z))
  step(fitP)
  fitS <- update(fitP, Strauss(7))
  step(fitS)
  fitM <- ppm(amacrine ~ marks*(x+y),
              MultiStrauss(types=levels(marks(amacrine)), radii=matrix(0.04, 2, 2)))
  step(fitM)
})
}



#'
#'   tests/sumfun.R
#'
#'   Tests of code for summary functions
#'       including score residual functions etc
#'
#'   $Revision: 1.6 $ $Date: 2020/05/01 09:59:59 $

if(ALWAYS) { # involves C code 
local({
  W <- owin(c(0,1), c(-1/2, 0))
  Gr <- Gest(redwood, correction="all",domain=W)
  Fr <- Fest(redwood, correction="all",domain=W)
  Jr <- Jest(redwood, correction="all",domain=W)
  
  F0 <- Fest(redwood[FALSE], correction="all")
  Fh <- Fest(humberside, domain=erosion(Window(humberside), 100))

  FIr <- Finhom(redwood, savelambda=TRUE)
  JIr <- Jinhom(redwood, savelambda=TRUE)
  
  Ga <- Gcross(amacrine, correction="all")
  Ia <- Iest(amacrine, correction="all")
  lam <- intensity(amacrine)
  lmin <- 0.9 * min(lam)
  nJ <- sum(marks(amacrine) == "off")
  FM <- FmultiInhom(amacrine, marks(amacrine) == "off",
                    lambdaJ=rep(lam["off"], nJ),
                    lambdamin = lmin)
  GM <- GmultiInhom(amacrine, marks(amacrine) == "on",
                    marks(amacrine) == "off",
                    lambda=lam[marks(amacrine)],
                    lambdamin=lmin,
                    ReferenceMeasureMarkSetI=42)

  pt <- psst(cells, interaction=Strauss(0.1), fun=nndcumfun)

  a <- compileCDF(D=nndist(redwood),
                  B=bdist.points(redwood),
                  r=seq(0, 1, length=256))

  ## distance argument spacing and breakpoints
  e <- check.finespacing(c(0,1,2), eps=0.1, action="silent")
  b <- as.breakpts(pi, 20)
  b <- as.breakpts(42, max=pi, npos=20)
  b <- even.breaks.owin(letterR)
})
}
