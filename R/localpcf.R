#
#   localpcf.R
#
#  $Revision: 1.27 $  $Date: 2022/05/21 08:53:38 $
#
#

localpcf <- function(X, ..., delta=NULL, rmax=NULL, nr=512, stoyan=0.15,
                     rvalue=NULL) {
  if(length(list(...)) > 0)
    warning("Additional arguments ignored")
  stopifnot(is.ppp(X))
  localpcfengine(X, delta=delta, rmax=rmax, nr=nr, stoyan=stoyan, rvalue=rvalue)
}

localpcfinhom <- function(X, ..., delta=NULL, rmax=NULL, nr=512, stoyan=0.15,
                          lambda=NULL, sigma=NULL, varcov=NULL,
                          update=TRUE, leaveoneout=TRUE,
                          rvalue=NULL) {
  stopifnot(is.ppp(X))
  a <- resolve.lambda(X, lambda, ...,
                      sigma=sigma, varcov=varcov,
                      update=update, leaveoneout=leaveoneout)
  result <- localpcfengine(X, ...,
                           delta=delta, rmax=rmax, nr=nr, stoyan=stoyan,
                           lambda=a$lambda, rvalue=rvalue)
  if(is.null(rvalue) && a$danger)
    attr(result, "dangerous") <- a$dangerous
  return(result)
}
 
localpcfengine <- function(X, ...,
                           delta=NULL, rmax=NULL, nr=512, stoyan=0.15,
                           lambda=NULL, rvalue=NULL) {
  if(!is.null(rvalue)) rmax <- rvalue
  m <- localpcfmatrix(X, delta=delta, rmax=rmax, nr=nr, stoyan=stoyan,
                      lambda=lambda)
  r <- attr(m, "r")
  if(!is.null(rvalue)) {
    nr <- length(r)
    if(r[nr] != rvalue)
      stop("Internal error - rvalue not attained")
    return(as.numeric(m[nr,]))
  }
  delta <- attr(m, "delta")
  nX <- npoints(X)
  if(nX == 0) {
    df <- data.frame(r=r, theo=rep.int(1, length(r)))
    nama <- desc <- labl <- NULL
  } else {
    # border correction
    dbord <- bdist.points(X)
    m[r[row(m)] > dbord[col(m)]] <- NA
    #
    df <- data.frame(m, r=r, theo=rep.int(1, length(r)))
    icode <- unlist(lapply(seq_len(nX), numalign, nmax=nX))
    nama <- paste("est", icode, sep="")
    desc <- paste("estimate of %s for point", icode)
    labl <- paste("%s[", icode, "](r)", sep="")
  }
  names(df) <- c(nama, "r", "theo")
  desc <- c(desc, "distance argument r", "theoretical Poisson %s")
  labl <- c(labl, "r", "%s[pois](r)")
  # create fv object
  g <- fv(df, "r", quote(localg(r)),
          "theo", , c(0, max(r)), labl, desc, fname="localg")
  # default is to display them all
  formula(g) <- . ~ r
  fvnames(g, ".") <- names(df)[names(df) != "r"]
  unitname(g) <- unitname(X)
  attr(g, "delta") <- delta  
  attr(g, "correction") <- "border"
  return(g)
}

localpcfmatrix <- function(X, i=seq_len(npoints(X)), ...,
                           lambda = NULL,
                           delta=NULL, rmax=NULL,
                           nr=512, stoyan=0.15) {
  missi <- missing(i)
  weighted <- !is.null(lambda)
  nX <- npoints(X)
  nY <- if(missi) nX else length(seq_len(nX)[i])
  W <- as.owin(X)
  lambda.ave <- nX/area(W)
  if(is.null(delta)) 
    delta <- stoyan/sqrt(lambda.ave)
  if(is.null(rmax)) 
    rmax <- rmax.rule("K", W, lambda.ave)
  #
  if(nX == 0 || nY == 0) {
    out <- matrix(0, nr, 0)
  } else {
    # sort points in increasing order of x coordinate
    oX <- fave.order(X$x)
    Xsort <- X[oX]
    idXsort <- (1:nX)[oX]
    if(weighted) {
      lambdaXsort <- lambda[oX]
      weightXsort <- 1/lambdaXsort
    }
    if(missi) {
      Y <- X
      oY <- oX
      Ysort   <- Xsort
      idYsort <- idXsort
    } else {
      # i is some kind of index
      Y <- X[i]
      idY <- (1:nX)[i]
      oY <- fave.order(Y$x)
      Ysort <- Y[oY]
      idYsort <- idY[oY]
    }
    nY <- npoints(Y)
    force(nr)
    # call C
    if(!weighted) {
      zz <- .C(SC_locpcfx,
               nn1 = as.integer(nY),
               x1  = as.double(Ysort$x),
               y1  = as.double(Ysort$y),
               id1 = as.integer(idYsort),
               nn2 = as.integer(nX),
               x2  = as.double(Xsort$x),
               y2  = as.double(Xsort$y),
               id2 = as.integer(idXsort),
               nnr = as.integer(nr),
               rmaxi=as.double(rmax),
               del=as.double(delta),
               pcf=as.double(double(nr * nY)),
               PACKAGE="spatstat.core")
    } else {
      zz <- .C(SC_locWpcfx,
               nn1 = as.integer(nY),
               x1  = as.double(Ysort$x),
               y1  = as.double(Ysort$y),
               id1 = as.integer(idYsort),
               nn2 = as.integer(nX),
               x2  = as.double(Xsort$x),
               y2  = as.double(Xsort$y),
               id2 = as.integer(idXsort),
               w2  = as.double(weightXsort),
               nnr = as.integer(nr),
               rmaxi=as.double(rmax),
               del=as.double(delta),
               pcf=as.double(double(nr * nY)),
               PACKAGE="spatstat.core")
    }
    out <- matrix(zz$pcf, nr, nY)
    # reorder columns to match original
    out[, oY] <- out
    # rescale
    out <- out/(2 * pi * if(!weighted) lambda.ave else 1)
  }
  # dress up
  attr(out, "r") <- seq(from=0, to=rmax, length.out=nr)
  attr(out, "delta") <- delta
  class(out) <- c("localpcfmatrix", class(out))
  return(out)
}

print.localpcfmatrix <- function(x, ...) {
  cat("Matrix of local pair correlation estimates\n")
  nc <- ncol(x)
  nr <- nrow(x)
  cat(paste("pcf estimates for", nc, ngettext(nc, "point", "points"), "\n"))
  rval <- attr(x, "r")
  cat(paste("r values from 0 to", max(rval), "in", nr, "steps\n"))
  return(invisible(NULL))
}

plot.localpcfmatrix <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  force(x)
  rval <- attr(x, "r")
  do.call(matplot,
          resolve.defaults(list(rval, quote(x)),
                           list(...),
                           list(type="l", main=xname,
                                xlab="r", ylab="pair correlation")))
}

"[.localpcfmatrix" <-
  function(x, i, ...) {
    r     <- attr(x, "r")
    delta <- attr(x, "delta")
    class(x) <- "matrix"
    if(missing(i)) {
      x <- x[ , ...]
    } else {
      x <- x[i, ...]
      if(is.matrix(i))
        return(x)
      r <- r[i]
    }
    if(!is.matrix(x))
      x <- matrix(x, nrow=length(r))
    attr(x, "r") <- r
    attr(x, "delta") <- delta
    class(x) <- c("localpcfmatrix", class(x))
    return(x)
}

