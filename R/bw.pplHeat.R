#'
#'  bw.pplHeat.R
#' 
#'  Bandwidth selection for densityHeat.ppp
#'  by point process likelihood cross-validation
#'
#'  Copyright (c) 2020 Adrian Baddeley, Tilman Davies and Suman Rakshit
#'  GNU Public Licence >= 2.0

bw.pplHeat <- function(X, ..., srange=NULL, ns=16, sigma=NULL,
                     leaveoneout=TRUE, verbose=TRUE) {
  
  #' compute intensity estimates
  b <- HeatEstimates.ppp(X, ..., srange=srange, ns=ns, sigma=sigma,
                         leaveoneout=leaveoneout, verbose=verbose)
  lambda <- b$lambda
  h      <- b$h
  hname  <- b$hname
  #' compute likelihood cross-validation criterion
  CV <- rowSums(log(lambda))
  iopt <- which.max(CV)
  result <- bw.optim(CV, h, iopt,
                     criterion="Likelihood cross-validation",
                     hname=hname,
                     unitname=unitname(X))
  return(result)
}

HeatEstimates.ppp <- function(X, ..., srange=NULL, ns=16, sigma=NULL,
                              leaveoneout=FALSE, verbose=TRUE) {
  stopifnot(is.ppp(X))
  nX <- npoints(X)

  ## trap a common error
  if(length(argh <- list(...)) &&
     (is.null(nama <- names(argh)) || !nzchar(nama[[1L]])) &&
     is.numeric(a <- argh[[1L]]) &&
     length(a) == 1L)
    stop("Use argument 'sigma' to specify the maximum bandwidth!",
         call.=FALSE)
  
  ## determine candidate bandwidths
  if(is.numeric(sigma) && length(sigma)) {
    ## sigma is a vector of candidate bandwidths, or a maximum bandwidth
    sigMax <- max(sigma)
    fractions <- if(length(sigma) > 1) sigma/sigMax else 
                 geomseq(from=0.05, to=1, length.out=ns)
  } else if(is.im(sigma)) {
    ## sigma is an image giving the spatially-varying maximum bandwidth
    sigMax <- sigma
    fractions <- seq_len(ns)/ns
  } else if(is.null(sigma)) {
    #' make a sequence of candidate bandwidths
    if(!is.null(srange)) {
      check.range(srange)
    } else {
      nnd <- nndist(X)
      srange <- c(min(nnd[nnd > 0]), diameter(as.owin(X))/2)
    }
    sigMax <- srange[2]
    sigmavalues <- geomseq(from=srange[1L], to=srange[2L], length.out=ns)
    fractions <- sigmavalues/sigMax
  } else stop("Format of sigma is not understood")

  ## set up transition matrix and initial state
  a <- densityHeat.ppp(X, sigMax, ..., internal=list(setuponly=TRUE))
  Y      <- a$Y     # initial state image
  u      <- a$u     # initial state vector (dropping NA)
  Xpos   <- a$Xpos  # location of each data point, index in 'u'
  A      <- a$A     # transition matrix, operates on 'u;
  Nstep  <- a$Nstep # total number of iterations

  ## map desired sigma values to iteration numbers
  nits <- pmax(1L, pmin(Nstep, round(Nstep * fractions^2)))
  nits <- diff(c(0L,nits))

  reciprocalpixelarea <- with(Y, 1/(xstep * ystep))

  ## compute ....
  lambda <- matrix(nrow=ns, ncol=nX)
  if(!leaveoneout) {
    ## usual estimates
    for(k in seq_len(ns)) {
      for(l in seq_len(nits[k])) u <- u %*% A
      lambda[k, ] <- u[Xpos]
    }
  } else {
    ## compute leave-one-out estimates
    if(verbose) {
      cat("Processing", nX, "points ... ")
      pstate <- list()
    }
    for(i in seq_len(nX)) {
      ## initial state = X[-i]
      ui <- u
      Xposi <- Xpos[i]
      ui[Xposi] <- ui[Xposi] - reciprocalpixelarea
      ## run iterations, pausing at each sigma value
      for(k in seq_len(ns)) {
        for(l in seq_len(nits[k])) ui <- ui %*% A
        lambda[k, i] <- ui[Xposi]
      }
      if(verbose) pstate <- progressreport(i, nX, state=pstate)
    }
    if(verbose) cat("Done.\n")
  }

  if(!is.im(sigma)) {
    h <- sigMax * fractions
    hname <- "sigma"
  } else {
    h <- fractions
    hname <- "fract"
  }

  return(list(lambda=lambda, h=h, hname=hname))
}

