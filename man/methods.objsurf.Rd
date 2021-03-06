\name{methods.objsurf}
\Rdversion{1.1}
\alias{methods.objsurf} %DoNotExport
\alias{print.objsurf}
\alias{plot.objsurf}
\alias{persp.objsurf}
\alias{image.objsurf}
\alias{contour.objsurf}
\alias{summary.objsurf}
\alias{print.summary.objsurf}
\title{
  Methods for Objective Function Surfaces
}
\description{
  Methods for printing and plotting an objective function surface.
}
\usage{
\method{print}{objsurf}(x, ...)
\method{plot}{objsurf}(x, ...)
\method{image}{objsurf}(x, ...)
\method{contour}{objsurf}(x, ...)
\method{persp}{objsurf}(x, ...)
\method{summary}{objsurf}(object, ...)
\method{print}{summary.objsurf}(x, ...)
}
\arguments{
  \item{x,object}{
    Object of class \code{"objsurf"} representing an objective function surface.
  }
  \item{\dots}{
    Additional arguments passed to plot methods.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{plot}},
  \code{\link{image}}, \code{\link{contour}}, \code{\link{persp}}
  and \code{\link{summary}} for the class \code{"objsurf"}.
}
\value{
  For \code{print.objsurf}, \code{print.summary.objsurf},
  \code{plot.objsurf} and \code{image.objsurf} the value is \code{NULL}.

  For \code{contour.objsurf} and \code{persp.objsurf}
  the value is described in the help for
  \code{\link{contour.default}} and \code{\link{persp.default}}
  respectively.

  For \code{summary.objsurf} the result is a list, of class
  \code{summary.objsurf}, containing summary information. This list is
  printed in sensible format by \code{print.summary.objsurf}.
}
\author{
  \adrian and \ege.
}
\seealso{
  \code{\link{objsurf}}
}
\examples{
 fit <- kppm(redwood ~ 1, "Thomas")
 os <- objsurf(fit)
  os
  summary(os)
  plot(os)
  contour(os, add=TRUE)
  persp(os)
}
\keyword{spatial}
\keyword{hplot}
