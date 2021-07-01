#'
#'   densityHeat.ppp.R
#'
#'   Diffusion estimator of density/intensity
#'

densityHeat <- function(x, sigma, ...) {
  UseMethod("densityHeat")
}

densityHeat.ppp <- function(x, sigma, ..., weights=NULL,
                            connect=8,
                            symmetric=FALSE, sigmaX=NULL, k=1,
                            show=FALSE, se=FALSE,
                            at=c("pixels", "points"),
                            leaveoneout = TRUE,
                            extrapolate = FALSE, coarsen = TRUE, 
                            verbose=TRUE,
                            internal=NULL) {
  stopifnot(is.ppp(x))
  nX <- npoints(x)
  at <- match.arg(at)
  if(length(weights)) check.nvector(weights, nX) else weights <- NULL

  if(extrapolate) {
    ## Richardson extrapolation 
    ## first compute intensity estimate on the desired grid
    cl <- sys.call()
    cl$extrapolate <- FALSE
    L <- eval(cl, sys.parent())
    dimL <- dim(L)
    ## remove all function arguments that control pixel resolution 
    cl$dimyx <- cl$eps <- cl$xy <- NULL
    if(coarsen) {
      ## compute on the desired grid and on a coarser grid
      Lfine <- L
      dimfine <- dimL
      ## compute on coarser grid
      dimcoarse <- round(dimfine/2)
      cl$dimyx <- dimcoarse
      Lcoarse <- eval(cl, sys.parent())
      ## interpolate coarse to fine
      Lcoarse <- as.im(interp.im, W=Window(Lfine), Z=Lcoarse, xy=Lfine,
                       bilinear=TRUE)
    } else {
      ## compute on the desired grid and a finer grid
      Lcoarse <- L
      dimcoarse <- dimL
      ## compute on finer grid
      dimfine <- round(dimcoarse * 2)
      cl$dimyx <- dimfine
      Lfine <- eval(cl, sys.parent())
      ## sample from fine to coarse
      Lfine <- as.im(Lfine, xy=Lcoarse)
    }
    ## Richardson extrapolation, ratio = 2, exponent = 1
    Lextrap <- 2 * Lfine - Lcoarse
    if(se)
      attr(Lextrap, "se") <- attr(L, "se")
    return(Lextrap)
  }

  delayed <- !is.null(sigmaX)
  setuponly <- identical(internal$setuponly, TRUE)
  want.Xpos <- delayed || setuponly

  if(!setuponly && (se || (at == "points" && leaveoneout))) {
    #' NEED INDIVIDUAL HEAT KERNELS FOR EACH DATA POINT
    #' to calculate estimate and standard error,
    #' or leave-one-out estimate
    if(!is.null(sigmaX))
      stop("variance calculation is not implemented for lagged arrivals")
    lambda <- varlam <- switch(at,
                               pixels = as.im(0, W=Window(x), ...),
                               points = numeric(nX))
    if(verbose) {
      pstate <- list()
      cat(paste("Processing", nX, "heat kernels... "))
    }
    if(is.null(weights)) {
      ## unweighted calculation: coded separately for efficiency
      for(i in seq_len(nX)) {
        Heat.i <- densityHeat.ppp(x[i], sigma, ...,
                                  connect=connect, symmetric=symmetric, k=k)
        switch(at,
               pixels = {
                 lambda <- lambda + Heat.i
                 varlam <- varlam + Heat.i^2
               },
               points = {
                 if(leaveoneout) {
                   Heat.ixi <- safelookup(Heat.i,x[-i],warn=FALSE)
                   #'was: Heat.ixi <- Heat.i[ x[-i] ]
                   lambda[-i] <- lambda[-i] + Heat.ixi
                   varlam[-i] <- varlam[-i] + Heat.ixi^2
                 } else {
                   lambda <- lambda + Heat.i[x]
                   varlam <- varlam + Heat.i[x]^2
                 }
               })
        if(verbose) pstate <- progressreport(i, nX, state=pstate)
      }
    } else {
      ## weighted calculation
      for(i in seq_len(nX)) {
        Heat.i <- densityHeat.ppp(x[i], sigma, ...,
                                  connect=connect, symmetric=symmetric, k=k)
        w.i <- weights[i]
        switch(at,
               pixels = {
                 lambda <- lambda + w.i * Heat.i
                 varlam <- varlam + w.i * Heat.i^2
               },
               points = {
                 if(leaveoneout) {
                   Heat.ixi <- Heat.i[ x[-i] ]
                   lambda[-i] <- lambda[-i] + w.i * Heat.ixi
                   varlam[-i] <- varlam[-i] + w.i * Heat.ixi^2
                 } else {
                   lambda <- lambda + w.i * Heat.i[x]
                   varlam <- varlam + w.i * Heat.i[x]^2
                 }
               })
        if(verbose) pstate <- progressreport(i, nX, state=pstate)
      }
    }
    if(verbose) splat("Done.")
    result <- lambda
    attr(result, "se") <- sqrt(varlam)
    return(result)
  }
  check.1.integer(k)
  stopifnot(k >= 1)
  if(!(connect %in% c(4,8)))
    stop("connectivity must be 4 or 8")

  ## initial state for diffusion
  if(delayed) {
    #' smoothing bandwidths attributed to each data point
    check.nvector(sigmaX, nX)
    stopifnot(all(is.finite(sigmaX)))
    stopifnot(all(sigmaX >= 0))
    if(missing(sigma)) sigma <- max(sigmaX) else check.1.real(sigma)
    #' sort in decreasing order of bandwidth
    osx <- order(sigmaX, decreasing=TRUE)
    sigmaX <- sigmaX[osx]
    x <- x[osx]
    #' discretise window
    W <- do.call.matched(as.mask,
                         resolve.defaults(list(...),
                                          list(w=Window(x))))
    #' initial state is zero
    Y <- as.im(W, value=0)
    #' discretised coordinates
    Xpos <- nearest.valid.pixel(x$x, x$y, Y)
  } else {
    #' pixellate pattern
    Y <- pixellate(x, ..., weights=weights, preserve=TRUE, savemap=want.Xpos)
    Xpos <- attr(Y, "map")
  } 

  #' validate sigma
  if(is.im(sigma)) {
    # ensure Y and sigma are on the same grid
    A <- harmonise(Y=Y, sigma=sigma)
    Y <- A$Y
    sigma <- A$sigma
  } else if(is.function(sigma)) {
    sigma <- as.im(sigma, as.owin(Y))
  } else {
    sigma <- as.numeric(sigma)
    check.1.real(sigma)
  }

  #' normalise as density 
  pixelarea <- with(Y, xstep * ystep)
  Y <- Y / pixelarea
  v <- as.matrix(Y)
  #' initial state
  u <- as.vector(v)

  if(want.Xpos) {
    #' map (row, col) to serial number 
    serial <- matrix(seq_len(length(v)), nrow(v), ncol(v))
    Xpos <- serial[as.matrix(as.data.frame(Xpos))]
  }
  
  #' symmetric random walk?
  if(symmetric) {
    asprat <- with(Y, ystep/xstep)
    if(abs(asprat-1) > 0.01)
      warning(paste("Symmetric random walk on a non-square grid",
                    paren(paste("aspect ratio", asprat))),
              call.=FALSE)
  }
  #' determine appropriate jump probabilities & time step
  pmax <- 1/(connect+1) # maximum permitted jump probability
  xstep <- Y$xstep
  ystep <- Y$ystep
  minstep <- min(xstep, ystep)
  if(symmetric) {
    #' all permissible transitions have the same probability 'pjump'.
    #' Determine Nstep, and dt=sigma^2/Nstep, such that
    #' Nstep >= 16 and M * pjump * minstep^2 = dt
    M <- if(connect == 4) 2 else 6
    Nstep <- max(16, ceiling(max(sigma)^2/(M * pmax * minstep^2)))    
    dt <- sn <- (sigma^2)/Nstep
    px <- py <- pxy <- sn/(M * minstep^2)
  } else {
    #' px is the probability of jumping 1 step to the right
    #' py is the probability of jumping 1 step up
    #' if connect=4, horizontal and vertical jumps are exclusive.
    #' if connect=8, horizontal and vertical increments are independent
    #' Determine Nstep, and dt = sigma^2/Nstep, such that
    #' Nstep >= 16 and 2 * pmax * minstep^2 = dt
    Nstep <- max(16, ceiling(max(sigma)^2/(2 * pmax * minstep^2)))
    dt <- sn <- (sigma^2)/Nstep
    px <- sn/(2 * xstep^2)
    py <- sn/(2 * ystep^2)
    if(max(px) > pmax) stop("Internal error: px exceeds pmax")
    if(max(py) > pmax) stop("Internal error: py exceeds pmax")
    if(connect == 8) pxy <- px * py
  }
  #' arrival times
  if(!is.null(sigmaX)) 
    iarrive <- pmax(1, pmin(Nstep, Nstep - round((sigmaX^2)/sn)))
  #' construct adjacency matrices
  dimv <- dim(v)
  my <- gridadjacencymatrix(dimv, across=FALSE, down=TRUE, diagonal=FALSE)
  mx <- gridadjacencymatrix(dimv, across=TRUE,  down=FALSE, diagonal=FALSE)
  if(connect == 8)
    mxy <- gridadjacencymatrix(dimv, across=FALSE,  down=FALSE, diagonal=TRUE)
  #' restrict to window
  if(anyNA(u)) {
    ok <- !is.na(u)
    u <- u[ok]
    if(want.Xpos) {
      #' adjust serial numbers
      Xpos <- cumsum(ok)[Xpos]
    }
    mx <- mx[ok,ok,drop=FALSE]
    my <- my[ok,ok,drop=FALSE]
    if(connect == 8) 
      mxy <- mxy[ok,ok,drop=FALSE]
    if(is.im(sigma)) {
      px <- px[ok]
      py <- py[ok]
      if(connect == 8) 
        pxy <- pxy[ok]
    }
  } else {
    ok <- TRUE
    if(is.im(sigma)) {
      px <- px[]
      py <- py[]
      if(connect == 8) pxy <- pxy[]
    }
  }
  #' construct iteration matrix
  if(connect == 4) {
    A <- px * mx + py * my
  } else {
    A <- px * (1 - 2 * py) * mx + py * (1 - 2 * px) * my + pxy * mxy
  }
  #' debug
  stopifnot(min(rowSums(A)) >= 0)
  stopifnot(max(rowSums(A)) <= 1)
  #' 
  diag(A) <- 1 - rowSums(A)
  #' k-step transition probabilities
  if(k > 1) {
    Ak <- A
    for(j in 2:k) Ak <- Ak %*% A
  } else Ak <- A
  k <- as.integer(k)
  Nstep <- as.integer(Nstep)
  Nblock <- Nstep/k
  Nrump  <- Nstep - Nblock * k
  #' secret exit - return setup data only
  if(setuponly)
    return(list(Y=Y, u=u, Xpos=Xpos, sigma=sigma, A=A, Ak=Ak, k=k,
                Nstep=Nstep, Nblock=Nblock, Nrump=Nrump,
                dx=xstep, dy=ystep, dt=dt))
  #' run
  U <- u
  Z <- Y
  if(!delayed) {
    if(!show) {
      for(iblock in 1:Nblock) U <- U %*% Ak
    } else {
      opa <- par(ask=FALSE)
      each <- max(1, round(Nblock/60))
      for(iblock in 1:Nblock) {
        U <- U %*% Ak
        if(iblock %% each == 0) {
          Z[] <- as.vector(U)
          f <- sqrt((iblock * k)/Nstep)
          main <- if(is.im(sigma)) paste(signif(f, 3), "* sigma") else
                  paste("sigma =", signif(f * sigma, 3))
          plot(Z, main=main)
          Sys.sleep(0.4)
        }
      }
      par(opa)
    }
    if(Nrump > 0) for(istep in 1:Nrump) U <- U %*% A
  } else {
    #' lagged arrivals
    used <- rep(FALSE, nX)
    contrib <- (weights %orifnull% rep(1,nX))/pixelarea
    if(!show) {
      for(iblock in 1:Nblock) {
        U <- U %*% Ak
        if(any(ready <- (!used & (iarrive <= iblock * k)))) {
          #' add points
          for(i in which(ready)) {
            j <- Xpos[i]
            U[j] <- U[j] + contrib[i]
            used[i] <- TRUE
          }
        }
      }
    } else {
      opa <- par(ask=FALSE)
      each <- max(1, round(Nblock/60))
      for(iblock in 1:Nblock) {
        U <- U %*% Ak
        if(any(ready <- (!used & (iarrive <= iblock * k)))) {
          #' add points
          for(i in which(ready)) {
            j <- Xpos[i]
            U[j] <- U[j] + contrib[i]
            used[i] <- TRUE
          }
        }
        if(iblock %% each == 0) {
          Z[] <- as.vector(U)
          f <- sqrt((iblock * k)/Nstep)
          main <- if(is.im(sigma)) paste(signif(f, 3), "* sigma") else
                  paste("sigma =", signif(f * sigma, 3))
          plot(Z, main=main)
          Sys.sleep(0.4)
        }
      }
      par(opa)
    }
    if(Nrump > 0) for(istep in 1:Nrump) U <- U %*% A
  }
  #' pack up
  Z[] <- as.vector(U)
  if(at == "points") Z <- Z[x]
  return(Z)
}
