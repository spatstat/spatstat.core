#
#
#  density.psp.R
#
#  $Revision: 1.21 $    $Date: 2022/05/21 08:55:03 $
#
#

density.psp <- function(x, sigma, ..., weights=NULL, edge=TRUE, 
                        method=c("FFT", "C", "interpreted"),
                        at=NULL) {
  verifyclass(x, "psp")
  method <- match.arg(method)
  w <- x$window
  n <- x$n
  if(length(weights)) {
    check.nvector(weights, n, things="segments", oneok=TRUE, vname="weights")
    if(length(weights) == 1) weights <- rep(weights, n)
  } else weights <- NULL
  len <- lengths_psp(x)
  ang <- angles.psp(x, directed=TRUE)
  ux <- unitname(x)
  if(missing(sigma))
    sigma <- 0.1 * diameter(w)
  #' determine locations for evaluation of density
  if(is.null(at)) {
    atype <- "window"
    w <- do.call.matched(as.mask, resolve.defaults(list(w=quote(w), ...)))
  } else if(is.owin(at)) {
    atype <- "window"
    w <- do.call.matched(as.mask, resolve.defaults(list(w=quote(at), ...)))
  } else {  
    atype <- "points"
    atY <- try(as.ppp(at, W=w))
    if(inherits(atY, "try-error"))
      stop("Argument 'at' should be a window or a point pattern", call.=FALSE)
  }
  #' detect empty pattern
  if(n == 0 || all(len == 0)) 
    switch(atype,
           window = return(as.im(0, w)),
           points = return(rep(0, npoints(atY))))
  #' determine prediction coordinates
  switch(atype,
         window = {
           xy <- rasterxy.mask(w)
           xx <- xy$x
           yy <- xy$y
         },
         points = {
           xx <- atY$x
           yy <- atY$y
         })
  #' c o m p u t e
  switch(method,
         interpreted = {
           #' compute matrix contribution from each segment 
           coz <- cos(ang)
           zin <- sin(ang)
           if(is.null(weights)) {
             #' unweighted
             for(i in seq_len(n)) {
               en <- x$ends[i,]
               dx <- xx - en$x0
               dy <- yy - en$y0
               u1 <- dx * coz[i] + dy * zin[i]
               u2 <- - dx * zin[i] + dy * coz[i]
               value <- dnorm(u2, sd=sigma) *
                 (pnorm(u1, sd=sigma) - pnorm(u1-len[i], sd=sigma))
               totvalue <- if(i == 1L) value else (value + totvalue)
             }
           } else {
             #' weighted
             for(i in seq_len(n)) {
               en <- x$ends[i,]
               dx <- xx - en$x0
               dy <- yy - en$y0
               u1 <- dx * coz[i] + dy * zin[i]
               u2 <- - dx * zin[i] + dy * coz[i]
               value <- weights[i] * dnorm(u2, sd=sigma) *
                 (pnorm(u1, sd=sigma) - pnorm(u1-len[i], sd=sigma))
               totvalue <- if(i == 1L) value else (value + totvalue)
             }
           }
           dens <- switch(atype,
                          window = im(totvalue, w$xcol, w$yrow, unitname=ux),
                          points = totvalue)
         },
         C = {
           #' C implementation of the above
           xs <- x$ends$x0
           ys <- x$ends$y0
           xp <- as.numeric(as.vector(xx))
           yp <- as.numeric(as.vector(yy))
           np <- length(xp)
           if(is.null(weights)) {
             #' unweighted
             z <- .C(SC_segdens,
                     sigma = as.double(sigma),
                     ns    = as.integer(n),
                     xs    = as.double(xs),
                     ys    = as.double(ys),
                     alps  = as.double(ang),
                     lens  = as.double(len),
                     np    = as.integer(np),
                     xp    = as.double(xp), 
                     yp    = as.double(yp),
                     z     = as.double(numeric(np)),
                     PACKAGE="spatstat.core")
           } else {
             #' weighted
             z <- .C(SC_segwdens,
                     sigma = as.double(sigma),
                     ns    = as.integer(n),
                     xs    = as.double(xs),
                     ys    = as.double(ys),
                     alps  = as.double(ang),
                     lens  = as.double(len),
                     ws    = as.double(weights),
                     np    = as.integer(np),
                     xp    = as.double(xp), 
                     yp    = as.double(yp),
                     z     = as.double(numeric(np)),
                     PACKAGE="spatstat.core")
           }
           dens <- switch(atype,
                          window = im(z$z, w$xcol, w$yrow, unitname=ux),
                          points = z$z)
         },
         FFT = {
           Y <- pixellate(x, ..., weights=weights, DivideByPixelArea=TRUE)
           dens <- blur(Y, sigma, normalise=edge, bleed=FALSE, ...)
           if(atype == "points") dens <- dens[atY, drop=FALSE]
         })
  if(edge && method != "FFT") {
    edg <- second.moment.calc(midpoints.psp(x), sigma, what="edge", ...)
    switch(atype,
           window = {
             dens <- eval.im(dens/edg)
           },
           points = {
             edgY <- edg[atY, drop=FALSE]
             dens <- dens/edgY
           })
  }
  if(atype == "window")
    dens <- dens[x$window, drop=FALSE]
  attr(dens, "sigma") <- sigma
  return(dens)
}
