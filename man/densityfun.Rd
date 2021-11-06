\name{densityfun.ppp}
\alias{densityfun}
\alias{densityfun.ppp}
\title{
  Kernel Estimate of Intensity as a Spatial Function
}
\description{
  Compute a kernel estimate of intensity for a point pattern,
  and return the result as a function of spatial location.  
}
\usage{
densityfun(X, \dots)

\method{densityfun}{ppp}(X, sigma = NULL, \dots,
                        weights = NULL, edge = TRUE, diggle = FALSE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{sigma}{
    Smoothing bandwidth, or bandwidth selection function,
    passed to \code{\link{density.ppp}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.ppp}}.
  }
  \item{weights}{
    Optional vector of weights associated with the points of \code{X}.
  }
  \item{edge,diggle}{
    Logical arguments controlling the edge correction.
    Arguments passed to \code{\link{density.ppp}}.
  }
}
\details{
  The commands \code{densityfun} and \code{\link{density}}
  both perform kernel estimation of the intensity of a point pattern.
  The difference is that \code{\link{density}} returns a pixel image,
  containing the estimated intensity values at a grid of locations, while
  \code{densityfun} returns a \code{function(x,y)} which can be used
  to compute the intensity estimate at \emph{any} spatial locations
  with coordinates \code{x,y}.
  For purposes such as model-fitting it is more accurate to
  use \code{densityfun}.
}
\section{Using the result of \code{densityfun}}{
  If \code{f <- densityfun(X)}, where \code{X} is a two-dimensional
  point pattern, the resulting object \code{f}
  is a \code{function} in the \R language.

  By calling this function,
  the user can evaluate the estimated intensity at any desired spatial
  locations. 

  Additionally \code{f} belongs to other
  classes which allow it to be printed and plotted easily.

  The function \code{f} has arguments \code{x,y,drop}.
  \itemize{
    \item The arguments \code{x,y} of \code{f}
    specify the query locations. They can be numeric vectors
    of coordinates. Alternatively \code{x} can be a point pattern
    (or data acceptable to \code{\link{as.ppp}}) and \code{y} is omitted.
    The result of \code{f(x,y)} is a numeric vector giving the values of
    the intensity.
    \item 
    The argument \code{drop} of \code{f} specifies how to handle query locations
    which are outside the window of the original data.
    If \code{drop=TRUE} (the default), such locations are ignored.
    If \code{drop=FALSE}, a value of \code{NA} is returned
    for each such location.
  }
  Note that the smoothing parameters, such as the bandwidth
  \code{sigma}, are assigned when \code{densityfun} is executed.
  Smoothing parameters are fixed inside the function \code{f}
  and cannot be changed by arguments of \code{f}.
}
\value{
  A \code{function} with arguments \code{x,y,drop}.
  The function also belongs to the class \code{"densityfun"} which has
  methods for \code{print} and \code{\link{as.im}}.
  It also belongs to the class \code{"funxy"} which has methods
  for \code{plot}, \code{contour} and \code{persp}.
}
\seealso{
  \code{\link{density}}.

  To interpolate values observed at the points, use \code{\link{Smoothfun}}.
}
\examples{
  f <- densityfun(swedishpines)
  f
  f(42, 60)
  X <- runifpoint(2, Window(swedishpines))
  f(X)
  plot(f)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
