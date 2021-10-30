#'
#'  lurkslrm.R
#'
#'  Lurking variable plot for spatial logistic regression model
#'
#'  $Revision: 1.3 $ $Date: 2021/10/30 06:21:03 $
#'

lurking.slrm <- function(object, covariate,
                         type="raw",
                         cumulative=TRUE,
                         ..., 
                         plot.it=TRUE,
                         plot.sd=TRUE,
                         clipwindow=NULL,
                         rv=NULL,
                         envelope=FALSE, nsim=39, nrank=1,
                         typename,
                         covname, 
                         oldstyle=FALSE, check=TRUE,
                         verbose=TRUE,
                         nx=128,
                         splineargs=list(spar=0.5),
                         internal=NULL) {
  cl <- match.call()
  clenv <- parent.frame()

  verifyclass(object, "slrm")

  ## default name for covariate
  if(missing(covname) || is.null(covname)) {
    co <- cl$covariate
    covname <- if(is.name(co)) as.character(co) else
               if(is.expression(co)) format(co[[1]]) else NULL
  }

  Xsim <- NULL
  if(!isFALSE(envelope)) {
    ## compute simulation envelope
    Xsim <- NULL
    if(!isTRUE(envelope)) {
      ## some kind of object
      Y <- envelope
      if(is.list(Y) && all(sapply(Y, is.ppp))) {
        Xsim <- Y
        envelope <- TRUE
      } else if(inherits(Y, "envelope")) {
        Xsim <- attr(Y, "simpatterns")
        if(is.null(Xsim))
          stop("envelope does not contain simulated point patterns")
        envelope <- TRUE
      } else stop("Unrecognised format of argument: envelope")
      nXsim <- length(Xsim)
      if(missing(nsim) && (nXsim < nsim)) {
        warning(paste("Only", nXsim, "simulated patterns available"))
        nsim <- nXsim
      }
    }
  }

  ## may need to refit the model
  if(is.expression(covariate)) {
    ## expression could involve variables that are not stored in object
    neednames <- all.vars(covariate)
    if(!all(neednames %in% colnames(object$Data$df))) 
      object <- update(object, save.all.vars=TRUE)
  }

  ## match type argument
  type <- pickoption("type", type,
                     c(raw="raw",
                       inverse="inverse",
                       pearson="pearson",
                       Pearson="pearson"))
  if(missing(typename))
    typename <- switch(type,
                       raw="raw residuals",
                       inverse="inverse-lambda residuals",
                       pearson="Pearson residuals")

  
  #################################################################
  ## extract data from fitted model
  Data <- object$Data
  ## original data pattern
  X <- Data$response
  ## spatial locations and weights used in fit
  df <- Data$df
  quadpoints <- ppp(df$x, df$y, window=Window(X))
  Z <- as.logical(df[,1])
  datapoints <- quadpoints[Z]
  wts <- exp(df[,"logpixelarea"])
  
  #################################################################
  ## compute the covariate
  if(is.im(covariate)) {
    covvalues <- covariate[quadpoints, drop=FALSE]
    covrange <- internal$covrange %orifnull% range(covariate, finite=TRUE)
  } else if(is.vector(covariate) && is.numeric(covariate)) {
    covvalues <- covariate
    covrange <- internal$covrange %orifnull% range(covariate, finite=TRUE)
    if(length(covvalues) != npoints(quadpoints))
      stop("Length of covariate vector,", length(covvalues), "!=",
           npoints(quadpoints), ", number of quadrature points")
  } else if(is.expression(covariate)) {
    ## Expression involving covariates in the fitted object
    if(!is.null(object$Data$covariates)) {
      ## Expression may involve an external variable
      neednames <- all.vars(covariate)
      missingnames <- setdiff(neednames, colnames(df))
      if(length(missingnames)) {
        ## missing variables should be 'external'
        foundvars <- mget(missingnames,
                          parent.frame(),
                          ifnotfound=rep(list(NULL), length(missingnames)))
        bad <- sapply(foundvars, is.null)
        if(any(bad)) {
          nbad <- sum(bad)
          stop(paste(ngettext(nbad, "Variable", "Variables"),
                     commasep(sQuote(missingnames[bad])),
                     "not found"), call.=FALSE)
        }
        founddata <- mpl.get.covariates(foundvars, quadpoints)
        df <- cbind(df, founddata)
      }
    }
    ## Evaluate expression
    sp <- parent.frame()
    covvalues <- eval(covariate, envir=df, enclos=sp)
    covrange <- internal$covrange %orifnull% range(covvalues, finite=TRUE)
    if(!is.numeric(covvalues))
      stop("The evaluated covariate is not numeric")
  } else 
    stop(paste("The", sQuote("covariate"), "should be either",
               "a pixel image, an expression or a numeric vector"))

  #################################################################
  ## Secret exit
  if(identical(internal$getrange, TRUE))
    return(covrange)
    
  ################################################################
  ## Residuals/marks attached to appropriate locations.
  ## Stoyan-Grabarnik weights are attached to the data points only.
  ## Others (residuals) are attached to all quadrature points.
  resvalues <- 
    if(!is.null(rv)) rv
    else if(type=="eem") eem(object, check=check)
    else residuals(object, type=type, check=check)

  if(inherits(resvalues, "imlist")) {
    if(length(resvalues) > 1) 
      stop("Not implemented for vector-valued residuals")
    resvalues <- resvalues[[1]]
  }

  if(is.im(resvalues))
     resvalues <- resvalues[quadpoints]

  ## NAMES OF THINGS
  ## name of the covariate
  if(is.null(covname)) 
    covname <- if(is.expression(covariate)) covariate else "covariate"
  ## type of residual/mark
  if(missing(typename)) 
    typename <- if(!is.null(rv)) "rv" else ""

  clip <- !is.null(clipwindow)

  ## CALCULATE
  stuff <- LurkEngine(object=object,
                      type=type, cumulative=cumulative, plot.sd=plot.sd,
                      quadpoints=quadpoints,
                      wts=wts,
                      Z=Z,
                      subQset=TRUE,
                      covvalues=covvalues,
                      resvalues=resvalues,
                      clip=clip,
                      clipwindow=clipwindow,
                      cov.is.im=is.im(covariate),
                      covrange=covrange,
                      typename=typename,
                      covname=covname,
                      cl=cl, clenv=clenv,
                      oldstyle=oldstyle, check=check, verbose=verbose,
                      nx=nx, splineargs=splineargs,
                      envelope=envelope, nsim=nsim, nrank=nrank, Xsim=Xsim,
                      internal=internal)
    
  ## ---------------  PLOT ----------------------------------
  if(plot.it && inherits(stuff, "lurk")) {
    plot(stuff, ...)
    return(invisible(stuff))
  } else {
    return(stuff)
  }
}
