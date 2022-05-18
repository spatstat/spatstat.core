#'
#'     localKcross.R
#'
#'     original by Ege Rubak
#' 
#'     $Revision: 1.16 $ $Date: 2022/05/18 05:41:49 $

"localLcross" <- function(X, from, to, ..., rmax = NULL, correction = "Ripley") {
  localKcross(X, from, to, ..., rmax = rmax, correction = correction, wantL = TRUE)
}

"localLdot" <- function(X, from, ..., rmax = NULL, correction = "Ripley") {
  localKdot(X, from, ..., rmax = rmax, correction = correction, wantL = TRUE)
}

"localKcross" <-
  function(X, from, to, ..., rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL)
  {
    verifyclass(X, "ppp")
    if(!is.multitype(X, dfok=FALSE)) 
	    stop("Point pattern must be multitype")
    marx <- marks(X)
    if(missing(from))
      from <- levels(marx)[1]
    if(missing(to))
      to <- levels(marx)[2]
    I <- (marx == from)
    if(!any(I))
      stop(paste("No points have mark =", from))
    Iexplain <- paste("points having mark =", from)
    Ikey <- make.parseable(paste(from))
    if(from == to) {
      ## use Kest
      XI <- X[I]
      dont.complain.about(XI)
      result <- do.call(localK,
                        resolve.defaults(list(X=quote(XI),
                                              rmax=rmax,
                                              correction=correction,
                                              verbose=verbose,
                                              rvalue=rvalue),
                                         list(...)))
    } else {
      J <- (marx == to)
      if(!any(J))
        stop(paste("No points have mark =", to))
      Jexplain <- paste("points having mark =", to)
      Jkey <- make.parseable(paste(to))
      result <-localKmultiEngine(X, I, J, ...,
                                 Ikey=Ikey, Jkey=Jkey,
                                 Iexplain=Iexplain, Jexplain=Jexplain,
                                 rmax = rmax, correction=correction,
                                 verbose=verbose, rvalue=rvalue)
    }
    return(result)
  }

"localKdot" <- 
function(X, from, ..., rmax = NULL, correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
  	stop("Point pattern must be multitype")
  marx <- marks(X)
  if(missing(from)) from <- levels(marx)[1]
  
  I <- (marx == from)
  J <- rep.int(TRUE, X$n)  # i.e. all points
  Iexplain <- paste("points having mark =", from)
  Jexplain <- "points of any type"
  Ikey <- make.parseable(paste(from))
  Jkey <- "."
  
  if(!any(I)) stop(paste("No points have mark =", from))
	
  result <- localKmultiEngine(X, I, J, ...,
                              Iexplain=Iexplain, Jexplain=Jexplain,
                              Ikey=Ikey, Jkey=Jkey,
                              rmax = rmax, correction=correction,
                              verbose=verbose, rvalue=rvalue)
  attr(result, "indices") <- list(from=from)
  return(result)
}

"localKcross.inhom" <-
  function(X, from, to, lambdaFrom=NULL, lambdaTo=NULL, ..., rmax = NULL,
           correction = "Ripley",
           sigma=NULL, varcov=NULL,
           lambdaX=NULL, update=TRUE, leaveoneout=TRUE)
  {
    verifyclass(X, "ppp")
    if(!is.multitype(X, dfok=FALSE))
      stop("Point pattern must be multitype")
    marx <- marks(X)
    if(missing(from))
      from <- levels(marx)[1]
    if(missing(to))
      to <- levels(marx)[2]
    I <- (marx == from)
    J <- (marx == to)
    Iexplain <- paste("points having mark =", from)
    Jexplain <- paste("points having mark =", to)
    Ikey <- make.parseable(paste(from))
    Jkey <- make.parseable(paste(to))
    K <- localKmultiEngine(X, I, J, lambdaFrom, lambdaTo, ..., rmax = rmax,
                           Iexplain=Iexplain, Jexplain=Jexplain,
                           Ikey=Ikey, Jkey=Jkey,
                           correction=correction,
                           sigma=sigma, varcov=varcov,
                           lambdaX=lambdaX, update=update,
                           leaveoneout=leaveoneout)
    attr(K, "indices") <- list(from=from, to=to)
    return(K)
  }

localLcross.inhom <- function(X, from, to, lambdaFrom = NULL, lambdaTo = NULL, ..., rmax = NULL) {
  localKcross.inhom(X, from, to, lambdaFrom, lambdaTo, ..., rmax = rmax, wantL = TRUE)
}

"localKmultiEngine" <-
  function(X, from, to, lambdaFrom=NULL, lambdaTo=NULL, ...,
           rmax = NULL, wantL=FALSE,
           correction="Ripley", verbose=TRUE, rvalue=NULL,
           sigma=NULL, varcov=NULL,
           lambdaX=NULL, update=TRUE, leaveoneout=TRUE,
           Iexplain="points satisfying condition I",
           Jexplain="points satisfying condition J",
           Ikey="I",
           Jkey="J")
  {
    npts <- npoints(X)
    W <- Window(X)
    areaW <- area(W)
    lambda.ave <- npts/areaW
    
    from <- ppsubset(X, from)
    to <- ppsubset(X, to)
    if(is.null(from) || is.null(to))
      stop("from and to must be valid subset indices")
    
    if(!any(from)) stop("no points belong to subset from")
    if(!any(to)) stop("no points belong to subset to")
    
    X_from <- X[from]
    X_to <- X[to]
    
    n_from <- sum(from)
    n_to <- sum(to)
    
    lambdaFrom.ave <- n_from/areaW
    lambdaTo.ave <- n_to/areaW

    weighted <- !is.null(lambdaFrom) || !is.null(lambdaTo) || !is.null(lambdaX)
    if(weighted){
      lambdas <- resolve.lambda.cross(X, from, to, lambdaFrom, lambdaTo, ...,
                                      lambdaX = lambdaX,
                                      sigma = sigma, varcov = varcov,
                                      leaveoneout = leaveoneout,
                                      update = update,
                                      Iexplain=Iexplain,
                                      Jexplain=Jexplain)
      lambdaFrom <- lambdas$lambdaI
      lambdaTo <- lambdas$lambdaJ
    }
    
    if(is.null(rvalue)) 
      rmaxdefault <- rmax %orifnull% rmax.rule("K", W, lambda.ave)
    else {
      stopifnot(is.numeric(rvalue))
      stopifnot(length(rvalue) == 1)
      stopifnot(rvalue >= 0)
      rmaxdefault <- rvalue
    }
    breaks <- handle.r.b.args(NULL, NULL, W, rmaxdefault=rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
    
    correction.given <- !missing(correction)
    correction <- pickoption("correction", correction,
                             c(none="none",
                               isotropic="isotropic",
                               Ripley="isotropic",
                               trans="translate",
                               translate="translate",
                               translation="translate",
                               best="best"),
                             multi=FALSE)
    
    correction <- implemented.for.K(correction, W$type, correction.given)
    
    # recommended range of r values
    alim <- c(0, min(rmax, rmaxdefault))
    
    # identify all close pairs
    rmax <- max(r)
    close <- crosspairs(X_from, X_to, rmax)
    # close$i and close$j are serial numbers in X_from and X_to respectively;        
    # map them to original serial numbers in X
    orig <- seq_len(npts)
    imap <- orig[from]
    jmap <- orig[to]
    iX <- imap[close$i]
    jX <- jmap[close$j]
    # eliminate any identical pairs
    if(any(from & to)) {
      ok <- (iX != jX)
      if(!all(ok)) {
        close$i  <- close$i[ok]
        close$j  <- close$j[ok]
        close$d  <- close$d[ok]
        close$xi  <- close$xi[ok]
        close$xj  <- close$xj[ok]
        close$yi  <- close$yi[ok]
        close$yj  <- close$yj[ok]
      }
    }
    # extract information for these pairs (relative to orderings of X_from, X_to)
    DIJ <- close$d
    XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
    I <- close$i
    J <- close$j
    if(weighted) {
      ## lambdaI <- lambdaFrom[I] ## not used
      lambdaJ <- lambdaTo[J]
      ## weightI <- 1/lambdaI  ## not used
      weightJ <- 1/lambdaJ
    } 
    
    # initialise
    df <- as.data.frame(matrix(NA, length(r), n_from))
    labl <- desc <- character(n_from)
    
    if(verbose) state <- list()
    
    switch(correction,
           none={
             # uncorrected! For demonstration purposes only!
             for(i in 1:n_from) {
               ii <- (I == i)
               ## Below
               wh <- whist(DIJ[ii], breaks$val,
                           if(weighted) weightJ[ii] else NULL)  # no edge weights
               Knone <- cumsum(wh)
               ## Tweaking factor to express Kcross.inhom as unweighted average of local contrib.
               if(weighted) Knone <- Knone * lambdaFrom.ave/lambdaFrom[i]
               df[,i] <- Knone
               icode <- numalign(i, n_from)
               names(df)[i] <- paste("un", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("uncorrected estimate of %s",
                                "for point", icode)
               if(verbose) state <- progressreport(i, n_from, state=state)
               
             }
             if(!weighted) df <- df/lambdaTo.ave
           },
           translate={
             # Translation correction
             XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
             edgewt <- edge.Trans(XI, XJ, paired=TRUE)
             if(weighted)
               edgewt <- edgewt * weightJ
             for(i in 1:n_from) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
               Ktrans <- cumsum(wh)
               ## Tweaking factor to express Kcross.inhom as unweighted average of local contrib.
               if(weighted) Ktrans <- Ktrans * lambdaFrom.ave/lambdaFrom[i]
               df[,i] <- Ktrans
               icode <- numalign(i, n_from)
               names(df)[i] <- paste("trans", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("translation-corrected estimate of %s",
                                "for point", icode)
               if(verbose) state <- progressreport(i, n_from, state=state)
             }
             if(!weighted) df <- df/lambdaTo.ave
             h <- diameter(W)/2
             df[r >= h, ] <- NA
           },
           isotropic={
             # Ripley isotropic correction
             edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
             if(weighted)
               edgewt <- edgewt * weightJ
             for(i in 1:n_from) {
               ii <- (I == i)
               wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
               Kiso <- cumsum(wh)
               ## Tweaking factor to express Kcross.inhom as unweighted average of local contrib.
               if(weighted) Kiso <- Kiso * lambdaFrom.ave/lambdaFrom[i]
               df[,i] <- Kiso
               icode <- numalign(i, n_from)
               names(df)[i] <- paste("iso", icode, sep="")
               labl[i] <- makefvlabel(NULL, "hat", character(2), icode)
               desc[i] <- paste("Ripley isotropic correction estimate of %s", 
                                "for point", icode)
               if(verbose) state <- progressreport(i, n_from, state=state)
             }
             if(!weighted) df <- df/lambdaTo.ave
             h <- diameter(W)/2
             df[r >= h, ] <- NA
           })
    # transform values if L required
    if(wantL)
      df <- sqrt(df/pi)
    
    # return vector of values at r=rvalue, if desired
    if(!is.null(rvalue)) {
      nr <- length(r)
      if(r[nr] != rvalue)
        stop("Internal error - rvalue not attained")
      return(as.numeric(df[nr,]))
    }
    ## function value table required
    ## add r and theo
    df <- cbind(df,
                data.frame(r=r,
                           theo=if(wantL) r else (pi * r^2)))
    desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
    labl <- c(labl, c("r", "{%s[%s]^{pois}}(r)"))
    ## Handle 'dot' symbol
    if(identical(Jkey, ".")) {
      Jkeyname <- "symbol(\"\\267\")"
      Jkeylab  <- quote(dot)
      Jkeyexpr <- quote(symbol("\267"))
    } else Jkeyname <- Jkeylab <- Jkeyexpr <- Jkey
    ## Determine fv labels
    if(!wantL) {
      if(!weighted) {
        fnam <- c("K", paste0("list(loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(K[loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(K[list(loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      } else {
        fnam <- c("K", paste0("list(inhom,loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(K[inhom,loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(K[list(inhom,loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      }
    } else {
      if(!weighted) {
        fnam <- c("L", paste0("list(loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(L[loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(L[list(loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      } else {
        fnam <- c("L", paste0("list(inhom,loc,", Ikey, ",", Jkeyname, ")"))
        ylab <- substitute(L[inhom,loc,I,J](r), list(I=Ikey, J=Jkeylab))
        yexp <- substitute(L[list(inhom,loc,I,J)](r), list(I=Ikey, J=Jkeyexpr))
      }
    }
    # create fv object
    K <- fv(df, "r", ylab, "theo", , alim, labl, desc,
            fname=fnam, yexp=yexp)
    # default is to display them all
    formula(K) <- . ~ r
    unitname(K) <- unitname(X)
    attr(K, "correction") <- correction
    if(weighted && lambdas$danger)
      attr(K, "dangerous") <- lambdas$dangerous
    ### TEMPORARY HACK to save info about the "from" points
    attr(K, "Xfrom") <- X_from
    return(K)
  }

