#'
#'  distcdf.R
#'
#' cdf of |X1-X2| when X1,X2 are iid uniform in W, etc
#'
#'  $Revision: 1.19 $  $Date: 2022/02/12 09:07:38 $
#'

distcdf <- local({
  
  distcdf <- function(W, V=W, ..., dW=1, dV=dW, nr=1024,
                      regularise=TRUE, savedenom=FALSE, delta=NULL) {
    reflexive <- (missing(V) || is.null(V)) && (missing(dV) || is.null(dV))
    diffuse <- is.owin(W) && is.owin(V)
    uniformW <- is.null(dW) || identical(dW, 1)
    uniformV <- is.null(dV) || identical(dV, 1)
    uniform <- uniformW && uniformV

    if(is.owin(W)) {
      W <- as.mask(as.owin(W), ...)
      dW <- as.im(dW, W=W)
    } else if(is.ppp(W)) {
      if(uniformW) {
        #' discrete uniform distribution on W
        dW <- pixellate(W, ...)
      } else {
        #' dW should be a weight or vector of weights
        if(!is.vector(dW) || !is.numeric(dW))
          stop("If W is a point pattern, dW should be a vector of weights")
        if(length(dW) == 1L) {
          dW <- rep(dW, npoints(W))
        } else stopifnot(length(dW) == npoints(W))
        dW <- pixellate(W, weights=dW, ...)
      }
    } else stop("W should be a point pattern or a window")

    if(!reflexive) {
      if(is.owin(V)) {
        V <- as.mask(as.owin(V), ...)
        dV <- as.im(dV, W=V)
      } else if(is.ppp(V)) {
        if(uniformV) {
          #' discrete uniform distribution on V
          dV <- pixellate(V, ...)
        } else {
          #' dV should be a weight or vector of weights
          if(!is.vector(dV) || !is.numeric(dV))
            stop("If V is a point pattern, dV should be a vector of weights")
          if(length(dV) == 1L) {
            dV <- rep(dV, npoints(V))
          } else stopifnot(length(dV) == npoints(V))
          dV <- pixellate(V, weights=dV, ...)
        }
      } else stop("V should be a point pattern or a window")
      
      if(!uniformV && min(dV) < 0) 
        stop("Negative values encountered in dV")
    }

    #' compute
    if(diffuse && uniform) {
      #' uniform distributions on windows 
      g <- if(reflexive) setcov(W, ...) else setcov(W, V, ...)
    } else {
      g <- if(reflexive) imcov(dW) else imcov(dW, dV)
    }
    r <- as.im(function(x,y) { sqrt(x^2 + y^2) }, g)
    pix <- with(r, max(xstep, ystep))
    #' extract
    rvals <- as.vector(as.matrix(r))
    gvals <- as.vector(as.matrix(g))
    rmax <- max(rvals)
    #' histogram
    if(is.null(nr)) 
      nr <- max(1024, 512 * ceiling(rmax/(pix*512)))
    rgrid <- seq(0, rmax, length=nr)
    ## dr <- rmax/(nr-1)
    h <- whist(rvals, breaks=rgrid, weights=gvals/sum(gvals))
    ch <- c(0,cumsum(h))
    #' regularise at very short distances
    if(regularise) {
      pix <- with(r, max(xstep, ystep))
      suspect <- which(rgrid <= 10 * pix)
      reference <- which(rgrid <= 20 * pix)
      weigh <- pmin(seq_along(ch), min(reference))^2
      fit <- lm(ch ~ I(rgrid^2) + I(rgrid^3) - 1,
                subset=reference, weights=weigh)
      ch[suspect] <- predict(fit)[suspect]
      ## enforce cdf properties
      ch[1] <- 0
      ch <- cummax(ch)
    }
    #' ok
    result <- fv(data.frame(r=rgrid, f=ch),
                 "r", quote(CDF(r)),
                 "f", , range(rvals), c("r","%s(r)"),
                 c("Interpoint distance","Cumulative probability"),
                 fname="CDF")
    #' refine spacing, if required
    if(!is.null(delta))
      result <- refine(result, delta)
    #'
    if(savedenom) {
      denomW <- integral(dW)
      denomV <- if(reflexive) denomW else integral(dV)
      attr(result, "denom") <- denomW * denomV
    }
    return(result)
  }

  refine <- function(H, delta=NULL, verbose=FALSE, force=FALSE) {
    ## H is CDF of pairwise distances
    ## Ensure H has spacing at most 'delta'
    check.1.real(delta)
    stopifnot(is.finite(delta) && (delta > 0))
    rstep <- mean(diff(H$r))
    inflate <- rstep/delta
    if(verbose) 
      splat("delta=", delta, "rstep=", rstep, "inflate=", inflate)
    if(inflate > 1) {
      ## interpolate H
      if(verbose) {
        plot(H, xlim=c(0, R/2))
      }
      H <- interpCDF(H, n=ceiling(inflate))
      if(verbose) {
        plot(H, add=TRUE, xlim=c(0,R), col=2)
        splat("New rstep=", mean(diff(H$r)))
      }
      if(force) {
        ## force CDF to be nondecreasing and to start from 0
        Hr <- H[["f"]]
        Hr[1] <- 0
        Hr <- cummax(Hr)
        H[["f"]] <- Hr
      }
    }
    return(H)
  }
  
  interpCDF <- function(H, ..., method=c("smooth.spline", "loess"),
                        delta=NULL, n=NULL) {
    ## H is CDF of pairwise distance
    ## Interpolate H by smoothing H(r)/r^2
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

  distcdf
  
})
