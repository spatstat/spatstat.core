#
#  transect.R
#
# Line transects of pixel images
#
#  $Revision: 1.8 $  $Date: 2021/06/22 05:33:50 $
#

transect.im <- local({

  specify.location <- function(loc, rect) {
    lname <- short.deparse(substitute(loc))
    if(is.numeric(loc) && length(loc) == 2)
      return(list(x=loc[1], y=loc[2]))
    if(is.list(loc))
      return(xy.coords(loc))
    if(!(is.character(loc) && length(loc) == 1))
      stop(paste("Unrecognised format for", sQuote(lname)), call.=FALSE)
    xr <- rect$xrange
    yr <- rect$yrange
    switch(loc,
           bottomleft  = list(x=xr[1],    y=yr[1]),
           bottom      = list(x=mean(xr), y=yr[1]),
           bottomright = list(x=xr[2],    y=yr[1]),
           right       = list(x=xr[2],    y=mean(yr)),
           topright    = list(x=xr[2],    y=yr[2]),
           top         = list(x=mean(xr), y=yr[2]),
           topleft     = list(x=xr[1],    y=yr[2]),
           left        = list(x=xr[1],    y=mean(yr)),
           centre=,
           center      = list(x=mean(xr), y=mean(yr)),
           stop(paste("Unrecognised location",
                      sQuote(lname), "=", dQuote(loc)),
                call.=FALSE)
           )
  }

  transect.im <- 
    function(X, ..., from="bottomleft", to="topright",
             nsample=512, click=FALSE, add=FALSE, curve=NULL) {
      Xname <- short.deparse(substitute(X))
      Xname <- sensiblevarname(Xname, "X")
      stopifnot(is.im(X))
      check.1.integer(nsample)
      if(length(curve)) {
        ## parametric curve
        ## validate specification of curve
        check.named.list(curve, c("f", "tlim"), namopt=c("tname", "tdescrip"))
        stopifnot(is.function(curve$f))
        check.range(curve$tlim)
        ## extract info
        tlim <- curve$tlim
        tname <- curve$tname %orifnull% "t"
        tdescrip <- curve$tdescrip %orifnull% "curve parameter"
        tunits <- NULL
        ## create sample points along curve
        t <- seq(tlim[1L], tlim[2L], length.out=nsample)
        xy <- (curve$f)(t)
        if(is.null(dim(xy)))
          stop("curve$f() should return a matrix or data frame")
        if(ncol(xy) != 2L)
          stop("curve$f() should return a matrix or data frame with 2 columns")
        hasnames <- all(c("x", "y") %in% colnames(xy))
        x <- xy[, if(hasnames) "x" else 1L]
        y <- xy[, if(hasnames) "y" else 2L]
      } else {
        ## straight line transect
        if(click) {
          ## interactive
          if(!add) plot(X)
          from <- spatstatLocator(1)
          points(from)
          to <- spatstatLocator(1)
          points(to)
          segments(from$x, from$y, to$x, to$y)
        } else {
          ## data defining a line segment
          R <- as.rectangle(X)
          from <- specify.location(from, R)
          to   <- specify.location(to,   R)
        }
        ## create sample points along transect
        if(identical(from,to))
          stop(paste(sQuote("from"), "and", sQuote("to"),
                     "must be distinct points"), call.=FALSE)
        u <- seq(0,1,length.out=nsample)
        x <- from$x + u * (to$x - from$x)
        y <- from$y + u * (to$y - from$y)
        leng <- sqrt( (to$x - from$x)^2 +  (to$y - from$y)^2)
        t <- u * leng
        tname <- "t"
        tdescrip <- "distance along transect"
        tunits <- unitname(X)
      }
      ## look up pixel values (may be NA)
      v <- X[list(x=x, y=y), drop=FALSE]
      ## package into fv object
      df <- data.frame(t=t, v=v)
      colnames(df) <- c(tname, Xname)
      fv(df, argu = tname,
         ylab = substitute(Xname(tname),
                           list(Xname=as.name(Xname), tname=as.name(tname))),
         valu=Xname,
         labl = c(tname, paste0("%s", paren(tname))),
         desc = c(tdescrip, "pixel value of %s"),
         unitname = tunits, fname = Xname)
    }

  transect.im
})
