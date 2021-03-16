\name{densityHeat.ppp}
\alias{densityHeat.ppp}
\title{
  Diffusion Estimate of Point Pattern Intensity
}
\description{
  Computes the diffusion estimate of the intensity of a point pattern.
}
\usage{
\method{densityHeat}{ppp}(x, sigma, \dots, weights=NULL,
          connect=8, symmetric=FALSE,
          sigmaX=NULL, k=1, show=FALSE, se=FALSE,
          at=c("pixels", "points"),
          leaveoneout = TRUE,
          extrapolate = FALSE, coarsen = TRUE,
          verbose=TRUE, internal=NULL)
}
\arguments{
  \item{x}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{sigma}{
    Smoothing bandwidth. A single number giving the equivalent
    standard deviation of the smoother.
    Alternatively, a pixel image (class \code{"im"}) or a
    \code{function(x,y)} giving the spatially-varying bandwidth.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{pixellate.ppp}}
    controlling the pixel resolution.
  }
  \item{weights}{
    Optional numeric vector of weights associated with each point of
    \code{x}.
  }
  \item{connect}{
    Grid connectivity: either 4 or 8.
  }
  \item{symmetric}{
    Logical value indicating whether to \emph{force} the algorithm
    to use a symmetric random walk.
  }
  \item{sigmaX}{
    Numeric vector of bandwidths, one associated with each data point in
    \code{x}. See Details.
  }
  \item{k}{
    Integer. Calculations will be performed by repeatedly multiplying
    the current state by the \code{k}-step transition matrix.
  }
  \item{show}{
    Logical value indicating whether to plot successive iterations.
  }
  \item{se}{
    Logical value indicating whether to compute standard errors.
  }
  \item{at}{
    Character string specifying whether to compute values
    at a grid of pixels (\code{at="pixels"}, the default)
    or at the data points of \code{x} (\code{at="points"}).
  }
  \item{leaveoneout}{
    Logical value specifying whether to compute a leave-one-out
    estimate at each data point, when \code{at="points"}.
  }
  \item{extrapolate}{
    Logical value specifying whether to use Richardson extrapolation
    to improve the accuracy of the computation.
  }
  \item{coarsen}{
    Logical value, controlling the calculation performed when
    \code{extrapolate=TRUE}. See Details.
  }
  \item{verbose}{
    Logical value specifying whether to print progress reports.
  }
  \item{internal}{
    Developer use only.
  }
}
\details{
  This command computes a diffusion kernel estimate 
  of point process intensity from the observed point pattern \code{x}.
  
  The function \code{\link{densityHeat}} is generic,
  with methods for point patterns in two dimensions
  (class \code{"ppp"}) and point patterns on a linear network
  (class \code{"lpp"}). The function \code{densityHeat.ppp} described
  here is the method for class \code{"ppp"}. Given a two-dimensional
  point pattern \code{x}, it computes a diffusion kernel estimate
  of the intensity of the point process which generated \code{x}.

  Diffusion kernel estimates were developed
  by Botev et al (2010), Barry and McIntyre (2011) and
  Baddeley et al (2021). 
  
  Barry and McIntyre (2011) proposed an estimator for point process
  intensity based on a random walk on the pixel grid inside the
  observation window. Baddeley et al (2021) showed that the
  Barry-McIntyre method is a special case of the \emph{diffusion estimator} 
  proposed by Botev et al (2010).

  The original Barry-McIntyre algorithm assumes a symmetric random walk
  (i.e. each possible transition has the same probability \eqn{p})
  and requires a square pixel grid (i.e. equal
  spacing in the \eqn{x} and \eqn{y} directions). Their original
  algorithm is used if \code{symmetric=TRUE}. Use the \code{\dots}
  arguments to ensure a square grid: for example, the argument
  \code{eps} specifies a square grid with spacing \code{eps} units.

  The more general algorithm used here (Baddeley et al, 2021)
  does not require a square grid of pixels.
  If the pixel grid is not square, and if \code{symmetric=FALSE}
  (the default), then the random walk is not symmetric,
  in the sense that the probabilities of different jumps will be
  different, in order to ensure that the smoothing is isotropic.

  This implementation also includes two generalizations to
  the case of adaptive smoothing (Baddeley et al, 2021).

  In the first version of adaptive smoothing, the bandwidth is
  spatially-varying.
  The argument \code{sigma} should be a pixel image (class \code{"im"})
  or a \code{function(x,y)} specifying the bandwidth at each spatial
  location. The smoothing is performed by solving the 
  heat equation with spatially-varying parameters.

  In the second version of adaptive smoothing, each data point in
  \code{x} is smoothed using a separate bandwidth.
  The argument \code{sigmaX} should be a numeric vector
  specifying the bandwidth for each point of \code{x}.
  The smoothing is performed using the lagged arrival algorithm.
  The argument \code{sigma} can be omitted.

  If \code{extrapolate=FALSE} (the default), calculations are performed
  using the Euler scheme for the heat equation. 
  If \code{extrapolate=TRUE}, the accuracy of the result will be
  improved by applying Richardson extrapolation (Baddeley et al, 2021, Section
  4). After computing the intensity estimate using the Euler scheme
  on the desired pixel grid, another estimate is computed using the same
  method on another pixel grid, and the two estimates are combined by
  Richardson extrapolation to obtain a more accurate result.
  The second grid is coarser than the original grid if
  \code{coarsen=TRUE} (the default), and finer than the original grid
  if \code{coarsen=FALSE}. Setting \code{extrapolate=TRUE} increases
  computation time by 35\% if \code{coarsen=TRUE} and by 400\% if
  \code{coarsen=FALSE}.
}
\value{
  Pixel image (object of class \code{"im"}) giving the estimated
  intensity of the point process.

  If \code{se=TRUE}, the result has an attribute \code{"se"}
  which is another pixel image giving the estimated standard error.

  If \code{at="points"} then the result is a numeric vector
  with one entry for each point of \code{x}.
}
\seealso{
  \code{\link[spatstat.core]{density.ppp}} for the usual kernel estimator,
  and \code{\link[spatstat.core]{adaptive.density}} for the
  tessellation-based estimator.
}
\references{
  Baddeley, A., Davies, T., Rakshit, S., Nair, G. and McSwiggan, G. (2021)
  Diffusion smoothing for spatial point patterns.
  \emph{Statistical Science}, in press.
  
  Barry, R.P. and McIntyre, J. (2011)
  Estimating animal densities and home range in regions with irregular
  boundaries and holes: a lattice-based alternative to the kernel
  density estimator. \emph{Ecological Modelling} \bold{222}, 1666--1672.

  Botev, Z.I., Grotowski, J.F. and Kroese, D.P. (2010)
  Kernel density estimation via diffusion.
  \emph{Annals of Statistics} \bold{38}, 2916--2957.
}
\author{
  Adrian Baddeley and Tilman Davies.
}
\examples{
   online <- interactive()
   if(!online) op <- spatstat.options(npixel=32)

   X <- runifpoint(25, letterR)
   Z <- densityHeat(X, 0.2)
   if(online) {
     plot(Z, main="Diffusion estimator")
     plot(X, add=TRUE, pch=16)
     integral(Z) # should equal 25
   }

   Z <- densityHeat(X, 0.2, se=TRUE)
   Zse <- attr(Z, "se")
   if(online) plot(solist(estimate=Z, SE=Zse), main="")

   Zex <- densityHeat(X, 0.2, extrapolate=TRUE)

   ZS <- densityHeat(X, 0.2, symmetric=TRUE, eps=0.125)
   if(online) {
     plot(ZS, main="fixed bandwidth")
     plot(X, add=TRUE, pch=16)
   }

   sig <- function(x,y) { (x-1.5)/10 }
   ZZ <- densityHeat(X, sig)
   if(online) {
     plot(ZZ, main="adaptive (I)")
     plot(X, add=TRUE, pch=16)
   }

   sigX <- sig(X$x, X$y)
   AA <- densityHeat(X, sigmaX=sigX)
   if(online) {
     plot(AA, main="adaptive (II)")
     plot(X, add=TRUE, pch=16)
   }
   if(!online) spatstat.options(op)
}
\keyword{spatial}
\keyword{smooth}
