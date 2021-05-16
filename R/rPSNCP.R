#'   simulation of product shot-noise Cox process
#'   Original: (c) Abdollah Jalilian 2021
#'   Adapted to spatstat by Adrian Baddeley
#'   $Revision: 1.5 $ $Date: 2021/05/16 02:24:17 $

rPSNCP <- local({
  
  ## ===================================================================
  ## kernel functions
  ## ===================================================================
  
  bkernels <- list(
    ## Gaussian kernel with bandwidth omega
    Thomas = function(r, omega, ...){ 
      exp(- r^2/(2 * omega^2)) / (2 * pi * omega^2)
    },
    ## Variance-Gamma (Bessel) kernel
    ## with bandwidth omega and shape parameter nu.ker
    VarGamma = function(r, omega, nu.ker){
      stopifnot(nu.ker > -1/2)
      sigma2 <- 1 / (4 * pi * nu.ker * omega^2)
      u <- r/omega
      u <- ifelse(u > 0,
      (u^nu.ker) * besselK(u, nu.ker) / (2^(nu.ker - 1) * gamma(nu.ker)),
      1)
      return(abs(sigma2 * u))
    },
    ## Cauchy kernel with bandwith omega
    Cauchy = function(r, omega, ...){
      ((1 + (r / omega)^2)^(-1.5)) / (2 * pi * omega^2)
    }
    ## end of 'bkernels' list
    )
  
  ## ===================================================================
  ## simulating from the product shot-noise Cox processes
  ## ===================================================================

  ## simulation from the null model of independent shot-noise components
  rPSNCP0 <- function(lambda, kappa, omega, kernels=NULL, nu.ker=NULL, 
                      win=owin(), nsim=1, drop=TRUE,
                      ...,
                      cnames=NULL, 
                      epsth=0.001
                      # , mc.cores=1L
                      ) {
    m <- length(lambda)
    if ((length(kappa) != m) || length(omega) != m ) 
      stop("arguments kappa and omega must have the same length as lambda")
    if (is.null(kernels))
      kernels <- rep("Thomas", m)
    else if(length(kernels) != m)
      stop("length of argument 'kernels' must equal the number of components")
    if(is.null(nu.ker))
      nu.ker <- rep(-1/4, m)
    lambda <- as.list(lambda)
    if (is.null(cnames))
      cnames <- 1:m
  ## simulation from the null model of independent shot-noise components
    corefun0 <- function(dumm) {
      xp <- yp <- numeric(0)
      mp <- integer(0)
      for (i in 1:m) {
        mui <- lambda[[i]]/kappa[i]
        Xi <- switch(kernels[i], 
                     Thomas   = rThomas(kappa[i], scale=omega[i],
                                        mu=mui, win=win,
                                        ...),
                     Cauchy   =  rCauchy(kappa[i], scale=omega[i],
                                         mu=mui, win=win,
                                         thresh=epsth,
                                         ...),
                     VarGamma = rVarGamma(kappa[i], scale=omega[i],
                                          mu=mui, win=win,
                                          nu.ker=nu.ker[i], nu.pcf=NULL,
                                          thresh=epsth,
                                          ...))
        xp <- c(xp, Xi$x)
        yp <- c(yp, Xi$y)
        mp <- c(mp, rep.int(i, Xi$n))
      }
      mp <- factor(mp, labels=cnames)
      out <- ppp(xp, yp, window=win, marks=mp, check=FALSE)
      return(out)
    }
    ## outlist <- if (mc.cores == 1) lapply(1:nsim, corefun0) 
    ## else parallel::mclapply(1:nsim, corefun0, mc.cores=mc.cores)
    outlist <- lapply(1:nsim, corefun0)
    outlist <- simulationresult(outlist, nsim, drop)
    return(outlist)
  }
  
  # ===================================================================
  # simulation from the model
  rPSNCP <- function(lambda=rep(100, 4), kappa=rep(25, 4), omega=rep(0.03, 4), 
                     alpha=matrix(runif(16, -1, 3), nrow=4, ncol=4), 
                     kernels=NULL, nu.ker=NULL, win=owin(),
                     nsim=1, drop=TRUE, 
                     ...,
                     cnames=NULL, epsth=0.001
                                        # , mc.cores=1L
                     ) {
    m <- length(lambda)
    if ((length(kappa) != m) || length(omega) != m ) 
      stop("Arguments kappa and omega must have the same length as lambda")
    if (!all(dim(alpha) == c(m, m)))
      stop("Dimensions of matrix alpha are not correct")
    if (is.null(kernels))
      kernels <- rep("Thomas", m)
    else if(length(kernels) != m)
      stop("Length of argument kernels must equal the number of components")
    if (is.null(nu.ker))
      nu.ker <- rep(-1/4, m)
    diag(alpha) <- 0
    if(all(alpha == 0))
      return(rPSNCP0(lambda=lambda, kappa=kappa, omega=omega, kernels=kernels, 
                     nu.ker=nu.ker, win=win, nsim=nsim, cnames=cnames, 
                     ..., 
                     epsth=epsth
                                        # , mc.cores=mc.cores
                     ))
    
    lambda <- as.list(lambda)
    frame <- boundingbox(win)
    dframe <- diameter(frame)
    W <- as.mask(win, ...)
    Wdim <- dim(W)
    wx <- as.vector(raster.x(W))
    wy <- as.vector(raster.y(W))
    
    sigma <- rmax <- numeric(m)
    for (i in 1:m) {
      if(is.im(lambda[[i]])) 
        lambda[[i]] <- as.im(lambda[[i]], dimyx=Wdim, W=W)
      keri <- function(r){ bkernels[[kernels[i]]](r, omega[i], nu.ker[i]) }
      keri0 <- keri(0)
      sigma[i] <- kappa[i] / keri0
      kerithresh <- function(r){ keri(r) / keri0 - epsth}
      rmax[i] <- uniroot(kerithresh, lower = omega[i] / 2, 
                         upper = 5 * dframe)$root # 4 * omega[i] #
    }
    dilated <- grow.rectangle(frame, max(rmax))
    
    corefun <- function(idumm)
    {
      Phi <- lapply(kappa, rpoispp, win=dilated)
      fr <- vector("list", length=m)
      for (i in 1:m) {
        keri <- function(r){ bkernels[[kernels[i]]](r, omega[i], nu.ker[i]) }
        keri0 <- keri(0)
        Phii <- Phi[[i]]
        fr[[i]] <- keri(crossdist.default(wx, wy, Phii$x, Phii$y)) / keri0
      }
      
      if (is.null(cnames))
        cnames <- 1:m
      xp <- yp <- numeric(0)
      mp <- integer(0)
      for (i in 1:m) {
        Si <- rowSums(fr[[i]])  / sigma[i]
        E <- matrix(1, nrow=length(wx), ncol=m)
        for (j in (1:m)[-i]) {
          E[, j] <- apply(1 + alpha[j, i] * fr[[j]], 1, prod) * exp(-alpha[j, i] * sigma[j])
        }
        values <-  Si * apply(E, 1, prod)
        Lam <- im(values, xcol=W$xcol, yrow=W$yrow, unitname = unitname(W))
        rhoi <- lambda[[i]]
        Xi <- rpoispp(rhoi * Lam)
        xp <- c(xp, Xi$x)
        yp <- c(yp, Xi$y)
        mp <- c(mp, rep.int(i, Xi$n))
      }
      mp <- factor(mp, labels=cnames)
      simout <- ppp(xp, yp, window=win, marks=mp, check=FALSE)
      # attr(simout, "parents") <- Phi
      return(simout)
    }
    ## outlist <- if (mc.cores == 1) lapply(1:nsim, corefun) 
    ## else parallel::mclapply(1:nsim, corefun, mc.cores=mc.cores)
    outlist <- lapply(1:nsim, corefun)
    outlist <- simulationresult(outlist, nsim, drop)
    return(outlist)
  }
  
  rPSNCP
})

