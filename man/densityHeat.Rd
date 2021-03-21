\name{densityHeat}
\alias{densityHeat}
\title{
  Diffusion Estimate of Point Pattern Intensity
}
\description{
  Computes a diffusion estimate of intensity for a point pattern.
}
\usage{
densityHeat(x, sigma, \dots)
}
\arguments{
  \item{x}{
    Point pattern (object of class \code{"ppp"} or another class).
  }
  \item{sigma}{
    Smoothing bandwidth. Usually a single number giving the equivalent
    standard deviation of the smoother. 
  }
  \item{\dots}{
    Additional arguments depending on the method.
  }
}
\details{
  The generic function \code{densityHeat} computes an
  estimate of point process intensity using a diffusion kernel method.

  Further details depend on the class of point pattern \code{x}.
  See the help file for the appropriate method.
}
\value{
  Depends on the class of \code{x}.
}
\seealso{
  For two-dimensional point patterns (objects of class \code{"ppp"}),
  the diffusion kernel estimator is \code{\link{densityHeat.ppp}}.
  The usual kernel estimator is \code{\link{density.ppp}},
  and the tessellation-based estimator is
  \code{\link{adaptive.density}}.

%%  For point patterns on a linear network (objects of class
%%  \code{"lpp"}), see \code{\link[spatstat.linnet]{densityHeat.lpp}}.
}
\author{
  Adrian Baddeley and Tilman Davies.
}
\keyword{spatial}
\keyword{smooth}
