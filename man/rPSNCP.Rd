\name{rPSNCP}
\alias{rPSNCP}
\title{Simulate Product Shot-noise Cox Process}
\description{
  Generate a random multitype point pattern, a realisation of the
  product shot-noise Cox process.
}
\usage{
 rPSNCP(lambda=rep(100, 4), kappa=rep(25, 4), omega=rep(0.03, 4), 
        alpha=matrix(runif(16, -1, 3), nrow=4, ncol=4), 
        kernels=NULL, nu.ker=NULL, win=owin(), nsim=1, drop=TRUE,
        \dots,
        cnames=NULL, epsth=0.001)
%        , mc.cores=1L
}
\arguments{
  \item{lambda}{
    List of intensities of component processes. Either a
    numberic vector determining the constant (homogeneous) intensities 
    or a list of pixel images (objects of class \code{"im"}) determining 
    the (inhomogeneous) intensity functions of  component processes. 
    The length of \code{lambda} determines the number of component processes.
  }
  \item{kappa}{
    Numeric vector of intensities of the Poisson process of cluster centres
    for component processes. Must have the same size as \code{lambda}.
  }
  \item{omega}{
    Numeric vector of bandwidths of cluster dispersal kernels
    for component processes. Must have the same size as \code{lambda} 
    and \code{kappa}.
  }
  \item{alpha}{
    Matrix of interaction parameters. Square numeric matrix with the same 
    number of rows and columns as the length of \code{lambda},
    \code{kappa} and \code{omega}.
    All entries of \code{alpha} must be greater than -1.
  }
  \item{kernels}{
    Vector of character string determining the cluster dispersal kernels
    of component processes. Impleneted kernels are Gaussian 
    kernel (\code{"Thomas"}) with bandwidth \code{omega}, 
    Variance-Gamma (Bessel) kernel (\code{"VarGamma"}) with 
    bandwidth \code{omega} and shape parameter \code{nu.ker} 
    and Cauchy kernel (\code{"Cauchy"}) with bandwidth \code{omega}.
    Must have the same length as \code{lambda}, \code{kappa} and \code{omega}.
  }
  \item{nu.ker}{
    Numeric vector of bandwidths of shape parameters for Varaince-Gamma
    kernels. 
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{cnames}{
    Optional vector of character strings giving the names of
    the component processes.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[spatstat.geom]{as.mask}}
    to determine the pixel array geometry.
    See \code{\link[spatstat.geom]{as.mask}}.
  }
  \item{epsth}{
    Numerical threshold to determine the maximum interaction range for 
    cluster kernels.
    % See Details.
    % NO DETAILS ARE GIVEN!
  }
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
%  \item{mc.cores}{
%    Integer value indicating the number of cores for parallel computing using
%    \code{"mclapply"} function in the \pkg{parallel} package.
%  }
}
\value{
  A point pattern (an object of class \code{"ppp"}) if \code{nsim=1}, or a
     list of point patterns if \code{nsim > 1}.  Each point pattern is
     multitype (it carries a vector of marks which is a factor).
}
\details{
  This function generates a realisation of a product shot-noise Cox
  process (PSNCP). This is a multitype (multivariate) Cox point process 
  in which each element of the multivariate random intensity
  \eqn{\Lambda(u)} 
  of the process is obtained by 
  \deqn{
    \Lambda_i(u) = \lambda_i(u) S_i(u) \prod_{j \neq i} E_{ji}(u)
  }{
    Lambda[i](u) = lambda[i](u) S[i](u) prod[j != i] E[ji](u)
  }
  where \eqn{\lambda_i(u)}{\lambda[i](u)} is the intensity
  of component \eqn{i} of the process,
  \deqn{
    S_i(u) = \frac{1}{\kappa_{i}} \sum_{v \in \Phi_i} k_{i}(u - v)
  }{
    S[i](u) = 1 / (kappa[i]) sum[v in Phi[i]] k[i](u - v)
  }
  is the shot-noise random feild for component \eqn{i} and 
  \deqn{
      E_{ji}(u) = \exp(-\kappa_{j} \alpha_{ji} / k_{j}(0)) \prod_{v \in \Phi_{j}} {1 + \alpha_{ji} \frac{k_j(u-v)}{k_j(0)}}
    }{
      E[ji](u) = exp(-\kappa[j] \alpha[ji] / k[j](0)) prod[v in Phi[j]] (1 + alpha[ji] k[j](u-v) / k[j](0))
    }
  is a product field controlling impulses from the parent Poisson process 
  \eqn{\Phi_j}{\Phi[j]} with constant intensity \eqn{\kappa_j}{\kappa[j]} of 
  component process \eqn{j} on \eqn{\Lambda_i(u)}{\Lambda[i](u)}.
  Here \eqn{k_i(u)}{k[i](u)} is an isotropic kernel (probability
  density) function on \eqn{R^2} with bandwidth \eqn{\omega_i}{\omega[i]} and 
  shape parameter \eqn{\nu_i}{\nu[i]},
  and \eqn{\alpha_{ji}>-1}{\alpha[j,i] > -1} is the interaction parameter.
}
\seealso{
\code{\link{rmpoispp}},
\code{\link{rThomas}},
\code{\link{rVarGamma}},
\code{\link{rCauchy}},
\code{\link{rNeymanScott}}
}
\references{
  Jalilian, A., Guan, Y., Mateu, J. and Waagepetersen, R. (2015)
  Multivariate product-shot-noise Cox point process models. 
  \emph{Biometrics}  \bold{71}(4), 1022--1033.
}
\examples{
  online <- interactive()
  # Example 1: homogeneous components
  lambda <- c(250, 300, 180, 400)
  kappa <- c(30, 25, 20, 25)
  omega <- c(0.02, 0.025, 0.03, 0.02)
  alpha <- matrix(runif(16, -1, 1), nrow=4, ncol=4)
  if(!online) lambda <- lambda/10
  X <- rPSNCP(lambda, kappa, omega, alpha)
  if(online) {
    plot(X)
    plot(split(X))
  }

  #Example 2: inhomogeneous components
  z1 <- scaletointerval.im(bei.extra$elev, from=0, to=1)
  z2 <- scaletointerval.im(bei.extra$grad, from=0, to=1)
  if(!online) {
    ## reduce resolution to reduce check time
    z1 <- as.im(z1, dimyx=c(40,80))
    z2 <- as.im(z2, dimyx=c(40,80))
  } 
  lambda <- list(
         exp(-8 + 1.5 * z1 + 0.5 * z2),
         exp(-7.25 + 1 * z1  - 1.5 * z2),
         exp(-6 - 1.5 * z1 + 0.5 * z2),
         exp(-7.5 + 2 * z1 - 3 * z2))
  kappa <- c(35, 30, 20, 25) / (1000 * 500)
  omega <- c(15, 35, 40, 25)
  alpha <- matrix(runif(16, -1, 1), nrow=4, ncol=4)
  if(!online) lambda <- lapply(lambda, "/", e2=10)
  sapply(lambda, integral)
  X <- rPSNCP(lambda, kappa, omega, alpha, win = bei$window, dimyx=dim(z1))
  if(online) {
    plot(X)
    plot(split(X), cex=0.5)
  }
}
\author{Abdollah Jalilian.
  Modified by \spatstatAuthors.
}
\keyword{spatial}
\keyword{datagen}
