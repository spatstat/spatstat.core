#'   simulation of product shot-noise Cox process 
# as described in the following reference:
# Jalilian, A., Guan, Y., Mateu, J., & Waagepetersen, R. (2015). Multivariate product‐shot‐noise Cox point process models. Biometrics, 71(4), 1022-1033.

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
rPSNCP0 <- function(rho, kappa, omega, win=owin(), clusters=NULL, nsim=1,
                    nu.ker=NULL, ij=NULL, eps = NULL, dimyx = NULL, xy = NULL, epsth=0.001, mc.cores=1L)
{
  m <- length(rho)
  if ((length(kappa) != m) || length(omega) != m ) 
    stop("kappa and omega paramters must be of the same size.")
  if (is.null(clusters))
    clusters <- rep("Thomas", m)
  else if(length(clusters) != m)
    stop("clusters must be a vector of the size of the number of components.")
  if (is.null(nu.ker))
    nu.ker <- rep(-1/4, m)
  rho <- as.list(rho)
  if (is.null(ij))
    ij <- 1:m
  corefun0 <- function(dumm)
  {
    xp <- yp <- mp <- NULL
    for (i in 1:m)
    {
      rhoi <- rho[[i]]
      if (is.numeric(rhoi))
        mui <- rhoi / kappa[i]
      else if (is.im(rhoi)) 
        mui <- eval.im(rhoi / kappa[i])
      Xi <- switch(clusters[i], 
                   Thomas = rThomas(kappa[i], omega[i], mui, win=win),
                   Cauchy =  rCauchy(kappa[i], omega[i], mui, win=win, eps=epsth),
                   VarGamma = rVarGamma(kappa[i], nu.ker=nu.ker[i], omega[i], mui, win=win, eps=epsth, nu.pcf=NULL))
      xp <- c(xp, Xi$x)
      yp <- c(yp, Xi$y)
      mp <- c(mp, rep(ij[i], Xi$n))
    }
    
    out <- ppp(xp, yp, window=win, marks=as.factor(mp))
    #  attr(out, "parents") <- Phi
  }
  outlist <- if (mc.cores == 1) lapply(1:nsim, corefun0) 
  else mclapply(1:nsim, corefun0, mc.cores=mc.cores)
  if (nsim == 1) 
    return(outlist[[1]])
  names(outlist) <- paste("Simulation", 1:nsim)
  return(outlist)
}

# ===================================================================
# simulation from the model
rPSNCP <- function(rho, kappa, omega, alpha, win=owin(), clusters=NULL, nsim=1,
                   nu.ker=NULL, ij=NULL, eps = NULL, dimyx = NULL, xy = NULL, epsth=0.001, mc.cores=1L)
{
  m <- length(rho)
  if ((length(kappa) != m) || length(omega) != m ) 
    stop("kappa and omega paramters must be of the same size.")
  if (!all(dim(alpha) == c(m, m)))
    stop("alpha paramter is not a matrix of correct dimensions.")
  if (is.null(clusters))
    clusters <- rep("Thomas", m)
  else if(length(clusters) != m)
    stop("clusters must be a vector of the size of the number of components.")
  if (is.null(nu.ker))
    nu.ker <- rep(-1/4, m)
  diag(alpha) <- 0
  if (all(alpha == 0))
    return(rPSNCP0(rho=rho, kappa=kappa, omega=omega, win=win, clusters=clusters, nsim=nsim,
                   nu.ker=nu.ker, ij=ij, eps = eps, dimyx = dimyx, xy = xy, epsth=epsth, mc.cores=mc.cores))
  
  rho <- as.list(rho)
  frame <- boundingbox(win)
  dframe <- diameter(frame)
  W <- as.mask(win, eps = eps, dimyx = dimyx, xy = xy)
  wx <- as.vector(raster.x(W))
  wy <- as.vector(raster.y(W))
  
  sigma <- rmax <- numeric(m)
  for (i in 1:m)
  {
    if(is.im(rho[[i]])) 
      rho[[i]] <- as.im(rho[[i]], dimyx=W$dim, W=W)
    keri <- function(r){ bkernels[[clusters[i]]](r, omega[i], nu.ker[i]) }
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
      keri <- function(r){ bkernels[[clusters[i]]](r, omega[i], nu.ker[i]) }
      keri0 <- keri(0)
      fr[[i]] <- keri(crossdist.default(wx, wy, Phi[[i]]$x, Phi[[i]]$y)) / keri0
    }
    
    if (is.null(ij))
      ij <- 1:m
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
      rhoi <- rho[[i]]
      Xi <- rpoispp(eval.im(rhoi * Lam))
      xp <- c(xp, Xi$x)
      yp <- c(yp, Xi$y)
      mp <- c(mp, rep(ij[i], Xi$n))
    }
    
    simout <- ppp(xp, yp, window=win, marks=as.factor(mp))
    # attr(simout, "parents") <- Phi
  }
  outlist <- if (mc.cores == 1) lapply(1:nsim, corefun) 
  else mclapply(1:nsim, corefun, mc.cores=mc.cores)
  if (nsim == 1) 
    return(outlist[[1]])
  names(outlist) <- paste("Simulation", 1:nsim)
  return(outlist)
}

# ===================================================================
# Example 1: homogeneous
rho <- c(250, 300, 180, 400)
kappa <- c(30, 25, 20, 25)
omega <- c(0.02, 0.025, 0.03, 0.02)
alpha <- matrix(runif(16, -1, 1), nrow=4, ncol=4)
X <- rPSNCP(rho, kappa, omega, alpha)
plot(X)
plot(split.ppp(X))

#Example 2: inhomogeneous
rho <- list()
z1 <- scaletointerval.im(bei.extra$elev, from=0, to=1)
z2 <- scaletointerval.im(bei.extra$grad, from=0, to=1)
rho[[1]] <- exp(-8 + 1.5 * z1 + 0.5 * z2)
rho[[2]] <- exp(-7.25 + 1 * z1  - 1.5 * z2)
rho[[3]] <- exp(-6 - 1.5 * z1 + 0.5 * z2)
rho[[4]] <- exp(-7.5 + 2 * z1 - 3 * z2)
lapply(rho, function(o){ integral.im(o) })
kappa <- c(35, 30, 20, 25) / (1000 * 500)
omega <- c(15, 35, 40, 25)
alpha <- matrix(runif(16, -1, 1), nrow=4, ncol=4)
X <- rPSNCP(rho, kappa, omega, alpha, win = bei$window, dimyx=c(101, 201))
plot(X)
plot(split.ppp(X), cex=0.5)
