#'
#' evalcovar.R
#'
#'   evaluate covariate values at data points and at pixels
#'
#' $Revision: 1.35 $ $Date: 2020/06/14 10:37:28 $
#'

evalCovar <- function(model, covariate, ...) {
  UseMethod("evalCovar")
}

evalCovar.ppm <- local({

  evalCovar.ppm <- function(model, covariate, ...,
                            lambdatype=c("cif", "trend", "intensity"),
                            dimyx=NULL, eps=NULL,
                            interpolate=TRUE, jitter=TRUE, 
                            modelname=NULL, covname=NULL,
                            dataname=NULL, subset=NULL) {
    lambdatype <- match.arg(lambdatype)
    #' evaluate covariate values at data points and at pixels
    ispois <- is.poisson(model)
    csr <- ispois && is.stationary(model)
    #' determine names
    if(is.null(modelname))
      modelname <- if(csr) "CSR" else short.deparse(substitute(model))
    if(is.null(covname)) {
      covname <- singlestring(short.deparse(substitute(covariate)))
      if(is.character(covariate)) covname <- covariate
    }
    if(is.null(dataname))
      dataname <- model$Qname
    
    info <-  list(modelname=modelname, covname=covname,
                  dataname=dataname, csr=csr, ispois=ispois,
                  spacename="two dimensions")
  
    X <- data.ppm(model)
    W <- as.owin(model)

    #' explicit control of pixel resolution
    if(!is.null(dimyx) || !is.null(eps))
      W <- as.mask(W, dimyx=dimyx, eps=eps)

    if(!is.null(subset)) {
      #' restrict to subset if required
      X <- X[subset]
      W <- W[subset, drop=FALSE]
    }
    
    #' evaluate covariate 
    if(is.character(covariate)) {
      #' One of the characters 'x' or 'y'
      #' Turn it into a function.
      ns <- length(covariate)
      if(ns == 0) stop("covariate is empty")
      if(ns > 1) stop("more than one covariate specified")
      covname <- covariate
      covariate <- switch(covariate,
                          x=xcoordfun,
                          y=ycoordfun,
                          stop(paste("Unrecognised covariate",
                                     dQuote(covariate))))
    } 
  
    if(!is.marked(model)) {
      #' ...................  unmarked .......................
      if(is.im(covariate)) {
        type <- "im"
        if(!interpolate) {
          #' look up covariate values 
          ZX <- safelookup(covariate, X)
        } else {
          #' evaluate at data points by interpolation
          ZX <- interp.im(covariate, X$x, X$y)
          #' fix boundary glitches
          if(any(uhoh <- is.na(ZX)))
            ZX[uhoh] <- safelookup(covariate, X[uhoh])
        }
        #' covariate values for pixels inside window
        Z <- covariate[W, drop=FALSE]
        #' corresponding mask
        W <- as.owin(Z)
      } else if(is.function(covariate)) {
        type <- "function"
        #' evaluate exactly at data points
        ZX <- covariate(X$x, X$y)
        if(!all(is.finite(ZX)))
          warning("covariate function returned NA or Inf values")
        #' window
        W <- as.mask(W)
        #' covariate in window
        Z <- as.im(covariate, W=W)
        #' collapse function body to single string
        covname <- singlestring(covname)
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("The covariate should be",
                        "an image, a function(x,y)",
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
      #' values of covariate in window
      Zvalues <- as.vector(Z[W, drop=TRUE])
      #' corresponding fitted [conditional] intensity values
      lambda <- as.vector(predict(model, locations=W,
                                  type=lambdatype)[W, drop=TRUE])
      #' pixel area (constant)
      pixelarea <- with(Z, xstep * ystep)
    } else {
      #' ...................  marked .......................
      if(!is.multitype(model))
        stop("Only implemented for multitype models (factor marks)")
      marx <- marks(X, dfok=FALSE)
      possmarks <- levels(marx)
      npts <- npoints(X)
      #' single image: replicate 
      if(is.im(covariate)) {
        covariate <- rep(list(covariate), times=length(possmarks))
        names(covariate) <- as.character(possmarks)
      }
      #'
      if(is.list(covariate) && all(sapply(covariate, is.im))) {
        #' list of images
        type <- "im"
        if(length(covariate) != length(possmarks))
          stop("Number of images does not match number of possible marks")
        #' evaluate covariate at each data point 
        ZX <- numeric(npts)
        for(k in seq_along(possmarks)) {
          ii <- (marx == possmarks[k])
          covariate.k <- covariate[[k]]
          if(!interpolate) {
            #' look up covariate values 
            values <- safelookup(covariate, X)
          } else {
            #' interpolate
            values <- interp.im(covariate.k, x=X$x[ii], y=X$y[ii])
            #' fix boundary glitches
            if(any(uhoh <- is.na(values)))
              values[uhoh] <- safelookup(covariate.k, X[ii][uhoh])
          }
          ZX[ii] <- values
        }
        #' restrict covariate images to window 
        Z <- solapply(covariate, "[", i=W, drop=FALSE)
        #' extract pixel locations and pixel values
        Zframes <- lapply(Z, as.data.frame)
        #' covariate values at each pixel inside window
        Zvalues <- unlist(lapply(Zframes, getElement, name="value"))
        #' pixel locations 
        locn <- lapply(Zframes, getxy)
        #' tack on mark values
        for(k in seq_along(possmarks))
          locn[[k]] <- cbind(locn[[k]], data.frame(marks=possmarks[k]))
        loc <- do.call(rbind, locn)
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=loc, type=lambdatype)
        #' pixel areas
        pixelarea <- rep(sapply(Z, pixarea), sapply(Z, npixdefined))
      } else if(is.function(covariate)) {
        type <- "function"
        #' evaluate exactly at data points
        ZX <- functioncaller(x=X$x, y=X$y, m=marx, f=covariate, ...)
        #' functioncaller: function(x,y,m,f,...) { f(x,y,m,...) }
        #' same window
        W <- as.mask(W)
        #' covariate in window
        Z <- list()
        for(k in seq_along(possmarks))
          Z[[k]] <- as.im(functioncaller, m=possmarks[k], f=covariate, W=W, ...)
        Zvalues <- unlist(lapply(Z, pixelvalues))
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=W, type=lambdatype)
        lambda <- unlist(lapply(lambda, pixelvalues))
        if(length(lambda) != length(Zvalues))
          stop("Internal error: length(lambda) != length(Zvalues)")
        #' collapse function body to single string
        covname <- singlestring(covname)
        #' pixel areas
        pixelarea <- rep(sapply(Z, pixarea), sapply(Z, npixdefined))
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("For a multitype point process model,",
                        "the covariate should be an image, a list of images,",
                        "a function(x,y,m)", 
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
    }    
    #' ..........................................................

    #' apply jittering to avoid ties
    if(jitter) {
      nX <- length(ZX)
      dZ <- 0.3 * quantile(diff(sortunique(c(ZX, Zvalues))), 1/min(20, nX))
      ZX <- ZX + rnorm(nX, sd=dZ)
      Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
    }

    lambdaname <- if(is.poisson(model)) "intensity" else lambdatype
    lambdaname <- paste("the fitted", lambdaname)
    check.finite(lambda, xname=lambdaname, usergiven=FALSE)
    check.finite(Zvalues, xname="the covariate", usergiven=TRUE)

    #' lambda values at data points
    lambdaX <- predict(model, locations=X, type=lambdatype)

    #' lambda image(s)
    lambdaimage <- predict(model, locations=W, type=lambdatype)
    
    #' wrap up 
    values <- list(Zimage      = Z,
                   lambdaimage = lambdaimage,
                   Zvalues     = Zvalues,
                   lambda      = lambda,
                   lambdaX     = lambdaX,
                   weights     = pixelarea,
                   ZX          = ZX,
                   type        = type)
    return(list(values=values, info=info))
  }

  xcoordfun <- function(x,y,m){x}
  ycoordfun <- function(x,y,m){y}

  pixarea <- function(z) { z$xstep * z$ystep }
  npixdefined <- function(z) { sum(!is.na(z$v)) }
  pixelvalues <- function(z) { as.data.frame(z)[,3L] }
  getxy <- function(z) { z[,c("x","y")] }

  functioncaller <- function(x,y,m,f,...) {
    nf <- length(names(formals(f)))
    if(nf < 2) stop("Covariate function must have at least 2 arguments")
    value <- if(nf == 2) f(x,y) else if(nf == 3) f(x,y,m) else f(x,y,m,...)
    return(value)
  }
            
  
  evalCovar.ppm
})

