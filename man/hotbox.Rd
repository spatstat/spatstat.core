\name{hotbox}
\alias{hotbox}
\title{
  Heat Kernel for a Two-Dimensional Rectangle
}
\description{
  Calculate values of the heat kernel in a rectangle
  with insulated edges.
}
\usage{
hotbox(Xsource, Xquery, sigma,
       \dots, W=NULL, squared=FALSE, nmax=20)
}
\arguments{
  \item{Xsource}{
    Point pattern of sources of heat.
    Object of class \code{"ppp"}
    or convertible to a point pattern using \code{as.ppp(Xsource, W)}.
  }
  \item{Xquery}{
    Locations where the heat kernel value is required.
    An object of class \code{"ppp"} specifying query location points,
    or an object of class \code{"im"} or \code{"owin"} specifying a
    grid of query points.
  }
  \item{sigma}{
    Bandwidth for kernel. A single number.
  }
  \item{\dots}{
    Extra arguments (passed to \code{\link[spatstat.geom]{as.mask}})
    controlling the
    pixel resolution of the result, when \code{Xquery} is a window
    or an image.
  }
  \item{W}{
    Window (object of class \code{"owin"}) used to define the spatial domain
    when \code{Xsource} is not of class \code{"ppp"}.
  }
  \item{squared}{
    Logical value indicating whether to take the square
    of each heat kernel value, before summing over the source points.
  }
  \item{nmax}{
    Number of terms to be used from the infinite-sum expression
    for the heat kernel. A single integer.
  }
}
\details{
  This function computes the sum of heat kernels
  associated with each of the source points, evaluating them
  at each query location.

  The window for evaluation of the heat kernel must be a rectangle.
  
  The heat kernel in any region can be expressed as an infinite sum of terms
  associated with the eigenfunctions of the Laplacian. The heat kernel
  in a rectangle is the product of heat kernels for
  one-dimensional intervals on the horizontal and vertical axes. This
  function uses \code{\link[spatstat.geom]{hotrod}} to compute the
  one-dimensional heat kernels, truncating the infinite sum to the
  first \code{nmax} terms, and then calculates the two-dimensional heat
  kernel from each source point to each query location. If
  \code{squared=TRUE} these values are squared. Finally the values are
  summed over all source points to obtain a single value for each
  query location.
}
\value{
  If \code{Xquery} is a point pattern,
  the result is a numeric vector with one entry for each query point.

  If \code{Xquery} is an image or window, the result is
  a pixel image.
}
\seealso{
  \code{\link{densityHeat.ppp}}
}
\references{
  Baddeley, A., Davies, T., Rakshit, S., Nair, G. and McSwiggan, G. (2021)
  Diffusion smoothing for spatial point patterns.
  \emph{Statistical Science}, in press.
}
\author{
  Adrian Baddeley and Greg McSwiggan.
}
\examples{
  X <- runifpoint(10)

  Y <- runifpoint(5)
  hotbox(X, Y, 0.1)

  plot(hotbox(X, Window(X), 0.1))
  points(X, pch=16)
}
\keyword{math}
