\name{rjitter.psp}
\alias{rjitter.psp}
\title{Random Perturbation of Line Segment Pattern}
\description{
  Randomly pertubs a spatial pattern of line segments by applying
  independent random displacements to the segment endpoints.
}
\usage{
\method{rjitter}{psp}(X, radius, \dots, clip=TRUE, nsim=1, drop=TRUE)
}
\arguments{
  \item{X}{
    A point pattern on a linear network
    (object of class \code{"psp"}).
  }
  \item{radius}{
    Scale of perturbations. A positive numerical value.
    Each point will be displaced by a random distance,
    with maximum displacement equal to this value.
  }
  \item{\dots}{
    Ignored.
  }
  \item{clip}{
    Logical value specifying what to do if segments cross the boundary
    of the window. See Details.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a spatial pattern of line segments (class
    \code{"psp"}) rather than a list of length 1 containing this pattern.
  }
}
\details{
  The function \code{\link[spatstat.geom]{rjitter}} is generic. This function is
  the method for the class \code{"psp"} of line segment patterns.

  Each of the endpoints of each segment in \code{X} will be
  subjected to an independent random displacement. The displacement vectors are
  uniformly distributed in a circle of radius \code{radius}.

  If \code{clip=TRUE} (the default), segment endpoints
  are permitted to move to locations slightly outside the window of \code{X},
  and the resulting segments will be clipped to the window.
  If \code{clip=FALSE}, segment endpoints are conditioned to fall
  inside the window.

  If \code{nsim=1} and \code{drop=TRUE}, the result is another
  spatial pattern of line segments (object of class \code{"psp"}).
  Otherwise, the result is a list of \code{nsim} line segment patterns.
}
\value{
  A spatial pattern of line segments
  (object of class \code{"psp"})
  or a list of such patterns.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.geom]{rjitter}} for point patterns in two dimensions.
}
\examples{
   E <- edges(letterR)
   Window(E) <- owin(c(1.9, 4.1), c(0.5, 3.5))
   plot(rjitter(E, 0.1))
}
\keyword{spatial}
\keyword{datagen}
