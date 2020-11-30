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
#'  tests/aucroc.R
#'
#'  AUC and ROC code
#'
#'  $Revision: 1.6 $ $Date: 2020/11/02 06:26:45 $

local({
  if(FULLTEST) {
    fit <- kppm(redwood ~ I(y-x))
    a <- roc(fit)
    b <- auc(fit)
    fet <- ppm(amacrine~x+y+marks)
    d <- roc(fet)
    e <- auc(fet)
  }
})
## tests/cdf.test.R

local({
  AC <- split(ants, un=FALSE)$Cataglyphis
  AM <- split(ants, un=FALSE)$Messor
  DM <- distmap(AM)
  if(ALWAYS) {
    ## (1) check cdf.test with strange data
    ## Marked point patterns with some marks not represented
    ## should produce a warning, rather than a crash:
    cdf.test(AC, DM)
  }
  if(FULLTEST) {
    ## should be OK:
    cdf.test(unmark(AC), DM)
    cdf.test(unmark(AC), DM, "cvm")
    cdf.test(unmark(AC), DM, "ad")
    ## other code blocks
    cdf.test(finpines, "x")
  }
  if(FULLTEST) {
    ## (2) Monte Carlo test for Gibbs model
    fit <- ppm(cells ~ 1, Strauss(0.07))
    cdf.test(fit, "x", nsim=9)

    ## cdf.test.slrm
    fut <- slrm(japanesepines ~ x + y)
    Z <- distmap(japanesepines)
    cdf.test(fut, Z)
  }
})
#'    tests/circular.R
#'
#'    Circular data and periodic distributions
#'
#'    $Revision: 1.4 $  $Date: 2020/04/28 12:58:26 $


local({
  if(ALWAYS) {
    a <- pairorient(redwood, 0.05, 0.15, correction="none")
    rose(a)
  }
  if(FULLTEST) {
    b <- pairorient(redwood, 0.05, 0.15, correction="best")
    rose(b, start="N", clockwise=TRUE)
  }
  if(ALWAYS) {
    #' arcs on the circle 
    #'       (depends on numerical behaviour)
    set.seed(19171025)
    aa <- replicate(7, runif(1, 0, 2*pi) + c(0, runif(1, 0, pi)),
                    simplify=FALSE)
    bb <- circunion(aa)

    assertsingle <- function(x, a, id) {
      y <- circunion(x)
      if(length(y) != 1 || max(abs(y[[1]] - a)) > .Machine$double.eps)
        stop(paste("Incorrect result from circunion in case", id),
             call.=FALSE)
      invisible(NULL)
    }

    assertsingle(list(c(pi/3, pi), c(pi/2, 3*pi/2)),
                 c(pi/3, 3*pi/2),
                 1)
    assertsingle(list(c(0, pi/2), c(pi/4, pi)),
                 c(0,pi),
                 2)
    assertsingle(list(c(-pi/4, pi/2), c(pi/4, pi)),
                 c((2-1/4)*pi, pi),
                 3)
  }
})

  
#'
#'     contact.R
#'
#'   Check machinery for first contact distributions
#'
#'   $Revision: 1.6 $  $Date: 2020/04/28 12:58:26 $

local({
  if(ALWAYS) {
    #' reduce complexity
    Y <- as.mask(heather$coarse, dimyx=c(100, 50))
    
    X <- runifpoint(100, win = complement.owin(Y))
    G <- Gfox(X, Y)
    J <- Jfox(X, Y)

    Y <- as.polygonal(Y)
    X <- runifpoint(100, win = complement.owin(Y))
    G <- Gfox(X, Y)
    J <- Jfox(X, Y)

    op <- spatstat.options(exactdt.checks.data=TRUE)
    U <- exactdt(X)
    spatstat.options(op)
  }
})

reset.spatstat.options()
#'
#'   tests/contrib.R
#'
#'   Tests for user-contributed code in spatstat
#'
#'   $Revision: 1.2 $  $Date: 2020/04/28 12:58:26 $

local({
  #' Jinhom
  #' Marie-Colette van Lieshout and Ottmar Cronie
  X <- redwood3
  fit <- ppm(X ~ polynom(x,y,2))
  lam <- predict(fit)
  lamX <- fitted(fit, dataonly=TRUE)
  lmin <- 0.9 * min(lam)
  g1 <- Ginhom(X, lambda=fit, update=TRUE)
  if(FULLTEST) {
    g2 <- Ginhom(X, lambda=fit, update=FALSE, lmin = lmin)
    g3 <- Ginhom(X, lambda=lam,  lmin=lmin)
    g4 <- Ginhom(X, lambda=lamX, lmin=lmin)
  }
  if(ALWAYS) {
    f2 <- Finhom(X, lambda=fit, update=FALSE)
  }
  if(FULLTEST) {
    f1 <- Finhom(X, lambda=fit, update=TRUE)
    f3 <- Finhom(X, lambda=lam,  lmin=lmin)
  }
})
# tests/correctC.R
# check for agreement between C and interpreted code
# for interpoint distances etc.
# $Revision: 1.7 $ $Date: 2020/04/28 12:58:26 $

if(ALWAYS) { # depends on hardware
local({
  eps <- .Machine$double.eps * 4

  checkagree <- function(A, B, blurb) {
    maxerr <- max(abs(A-B))
    cat("Discrepancy", maxerr, "for", blurb, fill=TRUE)
    if(maxerr > eps) 
      stop(paste("Algorithms for", blurb, "disagree"))
    return(TRUE)
  }

  ## pairdist.ppp
  set.seed(190901)
  X <- rpoispp(42)
  dC <- pairdist(X, method="C")
  dR <- pairdist(X, method="interpreted")
  checkagree(dC, dR, "pairdist()")

  dCp <- pairdist(X, periodic=TRUE, method="C")
  dRp <- pairdist(X, periodic=TRUE, method="interpreted")
  checkagree(dCp, dRp, "pairdist(periodic=TRUE)")

  dCp2 <- pairdist(X, periodic=TRUE, squared=TRUE, method="C")
  dRp2 <- pairdist(X, periodic=TRUE, squared=TRUE, method="interpreted")
  checkagree(dCp2, dRp2, "pairdist(periodic=TRUE, squared=TRUE)")

  ## crossdist.ppp
  Y <- rpoispp(42)
  dC <- crossdist(X, Y, method="C")
  dR <- crossdist(X, Y, method="interpreted")
  checkagree(dC, dR, "crossdist()")

  dC <- crossdist(X, Y, periodic=TRUE, method="C")
  dR <- crossdist(X, Y, periodic=TRUE, method="interpreted")
  checkagree(dC, dR, "crossdist(periodic=TRUE)")

  dC2 <- crossdist(X, Y, periodic=TRUE, squared=TRUE, method="C")
  dR2 <- crossdist(X, Y, periodic=TRUE, squared=TRUE, method="interpreted")
  checkagree(dC2, dR2, "crossdist(periodic=TRUE, squared=TRUE)")

  # nndist.ppp
  nnC <- nndist(X, method="C")
  nnI <- nndist(X, method="interpreted")
  checkagree(nnC, nnI, "nndist()")

  nn3C <- nndist(X, k=3, method="C")
  nn3I <- nndist(X, k=3, method="interpreted")
  checkagree(nn3C, nn3I, "nndist(k=3)")

  # nnwhich.ppp
  nwC <- nnwhich(X, method="C")
  nwI <- nnwhich(X, method="interpreted")
  checkagree(nwC, nwI, "nnwhich()")

  nw3C <- nnwhich(X, k=3, method="C")
  nw3I <- nnwhich(X, k=3, method="interpreted")
  checkagree(nw3C, nw3I, "nnwhich(k=3)")

  # whist
  set.seed(98123)
  x <- runif(1000)
  w <- sample(1:5, 1000, replace=TRUE)
  b <- seq(0,1,length=101)
  op <- spatstat.options(Cwhist=TRUE)
  aT <- whist(x,b,w)
  spatstat.options(Cwhist=FALSE)
  aF <- whist(x,b,w)
  if(!all(aT == aF))
    stop("Algorithms for whist disagree")
  spatstat.options(op)
})

reset.spatstat.options()
}
