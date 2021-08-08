#'
#'  distcdf.R
#'
#' cdf of |X1-X2| when X1,X2 are iid uniform in W, etc
#'
#'  $Revision: 1.15 $  $Date: 2021/08/08 05:25:40 $
#'

distcdf <- function(W, V=W, ..., dW=1, dV=dW, nr=1024,
                    regularise=TRUE, savedenom=FALSE) {
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
  dr <- rmax/(nr-1)
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
  #'
  if(savedenom) {
    denomW <- integral(dW)
    denomV <- if(reflexive) denomW else integral(dV)
    attr(result, "denom") <- denomW * denomV
  }
  return(result)
}

