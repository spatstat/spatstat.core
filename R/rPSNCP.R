#'   simulation of product shot-noise Cox process 

rPSNCP <- local({
  
  # ===================================================================
  # kernels functions
  # ===================================================================
  
  bkernels <- list()
  # Gaussian kernel with bandwidth omega
  bkernels$Thomas <- function(r, omega, ...){ 
    exp(- r^2/(2 * omega^2)) / (2 * pi * omega^2)
  }
  # Variance-Gamma (Bessel) kernel with bandwidth omega and shape parameter nu.ker
  bkernels$VarGamma <- function(r, omega, nu.ker){
    stopifnot(nu.ker > -1/2)
    sigma2 <- 1 / (4 * pi * nu.ker * omega^2)
    u <- r/omega
    u <- ifelse(u > 0, (u^nu.ker) * besselK(u, nu.ker) / (2^(nu.ker - 1) * gamma(nu.ker)), 1)
    return(abs(sigma2 * u))
  }
  # Cauchy kernel with bandwith omega
  bkernels$Cauchy <- function(r, omega, ...){
    ((1 + (r / omega)^2)^(-1.5)) / (2 * pi * omega^2)
  }
  
  # ===================================================================
  # simulating from the product shot-noise Cox processes
  # ===================================================================
  
  # simulation from the null model of independent shot-noise components
  rPSNCP0 <- function(lambda, kappa, omega, kernels=NULL, nu.ker=NULL, 
                      win=owin(), nsim=1, names=NULL, 
                      eps=NULL, dimyx=NULL, xy=NULL, epsth=0.001, mc.cores=1L)
  {
    m <- length(lambda)
    if ((length(kappa) != m) || length(omega) != m ) 
      stop("kappa and omega paramters must be of the same size.")
    if (is.null(kernels))
      kernels <- rep("Thomas", m)
    else if(length(kernels) != m)
      stop("kernels must be a vector of the size of the number of components.")
    if (is.null(nu.ker))
      nu.ker <- rep(-1/4, m)
    lambda <- as.list(lambda)
    if (is.null(names))
      names <- 1:m
    corefun0 <- function(dumm)
    {
      xp <- yp <- mp <- NULL
      for (i in 1:m)
      {
        rhoi <- lambda[[i]]
        if (is.numeric(rhoi))
          mui <- rhoi / kappa[i]
        else if (is.im(rhoi)) 
          mui <- eval.im(rhoi / kappa[i])
        Xi <- switch(kernels[i], 
                     Thomas = rThomas(kappa[i], omega[i], mui, win=win),
                     Cauchy =  rCauchy(kappa[i], omega[i], mui, win=win, eps=epsth),
                     VarGamma = rVarGamma(kappa[i], nu.ker=nu.ker[i], omega[i], mui, win=win, eps=epsth, nu.pcf=NULL))
        xp <- c(xp, Xi$x)
        yp <- c(yp, Xi$y)
        mp <- c(mp, rep(names[i], Xi$n))
      }
      
      out <- ppp(xp, yp, window=win, marks=as.factor(mp))
      #  attr(out, "parents") <- Phi
    }
    outlist <- if (mc.cores == 1) lapply(1:nsim, corefun0) 
    else parallel::mclapply(1:nsim, corefun0, mc.cores=mc.cores)
    if (nsim == 1) 
      return(outlist[[1]])
    names(outlist) <- paste("Simulation", 1:nsim)
    return(outlist)
  }
  
  # ===================================================================
  # simulation from the model
  rPSNCP <- function(lambda=rep(100, 4), kappa=rep(25, 4), omega=rep(0.03, 4), 
                     alpha=matrix(runif(16, -1, 3), nrow=4, ncol=4), 
                     kernels=NULL, nu.ker=NULL, win=owin(), nsim=1, names=NULL, 
                     eps=NULL, dimyx=NULL, xy=NULL, epsth=0.001, mc.cores=1L)
  {
    m <- length(lambda)
    if ((length(kappa) != m) || length(omega) != m ) 
      stop("kappa and omega paramters must be of the same size.")
    if (!all(dim(alpha) == c(m, m)))
      stop("alpha paramter is not a matrix of correct dimensions.")
    if (is.null(kernels))
      kernels <- rep("Thomas", m)
    else if(length(kernels) != m)
      stop("kernels must be a vector of the size of the number of components.")
    if (is.null(nu.ker))
      nu.ker <- rep(-1/4, m)
    diag(alpha) <- 0
    if (all(alpha == 0))
      return(rPSNCP0(lambda=lambda, kappa=kappa, omega=omega, kernels=kernels, 
                     nu.ker=nu.ker, win=win, nsim=nsim, names=names, 
                     eps=eps, dimyx=dimyx, xy=xy, epsth=epsth, mc.cores=mc.cores))
    
    lambda <- as.list(lambda)
    frame <- boundingbox(win)
    dframe <- diameter(frame)
    W <- as.mask(win, eps = eps, dimyx = dimyx, xy = xy)
    wx <- as.vector(raster.x(W))
    wy <- as.vector(raster.y(W))
    
    sigma <- rmax <- numeric(m)
    for (i in 1:m)
    {
      if(is.im(lambda[[i]])) 
        lambda[[i]] <- as.im(lambda[[i]], dimyx=W$dim, W=W)
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
      for (i in 1:m)
      {
        keri <- function(r){ bkernels[[kernels[i]]](r, omega[i], nu.ker[i]) }
        keri0 <- keri(0)
        fr[[i]] <- keri(crossdist.default(wx, wy, Phi[[i]]$x, Phi[[i]]$y)) / keri0
      }
      
      if (is.null(names))
        names <- 1:m
      xp <- yp <- mp <- NULL
      for (i in 1:m)
      {
        Si <- rowSums(fr[[i]])  / sigma[i]
        E <- matrix(1, nrow=length(wx), ncol=m)
        for (j in (1:m)[-i])
        {
          E[, j] <- apply(1 + alpha[j, i] * fr[[j]], 1, prod) * exp(-alpha[j, i] * sigma[j])
        }
        values <-  Si * apply(E, 1, prod)
        Lam <- im(values, xcol=W$xcol, yrow=W$yrow, unitname = unitname(W))
        rhoi <- lambda[[i]]
        Xi <- rpoispp(eval.im(rhoi * Lam))
        xp <- c(xp, Xi$x)
        yp <- c(yp, Xi$y)
        mp <- c(mp, rep(names[i], Xi$n))
      }
      
      simout <- ppp(xp, yp, window=win, marks=as.factor(mp))
      # attr(simout, "parents") <- Phi
    }
    outlist <- if (mc.cores == 1) lapply(1:nsim, corefun) 
    else parallel::mclapply(1:nsim, corefun, mc.cores=mc.cores)
    if (nsim == 1) 
      return(outlist[[1]])
    names(outlist) <- paste("Simulation", 1:nsim)
    return(outlist)
  }
  
  rPSNCP
})
