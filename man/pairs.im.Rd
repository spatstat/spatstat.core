\name{pairs.im}
\alias{pairs.im}
\title{
  Scatterplot Matrix for Pixel Images
}
\description{
  Produces a scatterplot matrix of the pixel values
  in two or more pixel images.
}
\usage{
\method{pairs}{im}(..., plot=TRUE, drop=TRUE)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is either
    a pixel image (object of class \code{"im"})
    or a named argument to be passed to \code{\link{pairs.default}}.
    Alternatively, a single argument which is a list of pixel images.
  }
  \item{plot}{
    Logical. If \code{TRUE}, the scatterplot matrix is plotted.
  }
  \item{drop}{
    Logical value specifying whether pixel values that are \code{NA} 
    should be removed from the data frame that is returned by the function.
    This does not affect the plot.
  }
}
\details{
  This is a method for the generic function \code{\link{pairs}}
  for the class of pixel images.
  
  It produces a square array of plot panels, in which each panel shows
  a scatterplot of the pixel values of one image against the
  corresponding pixel values of another image.

  At least two of the arguments \code{\dots}
  should be pixel images (objects of class \code{"im"}).
  Their spatial domains must overlap, but need not have the same
  pixel dimensions.

  First the pixel image domains are intersected, and converted to a
  common pixel resolution. Then the corresponding pixel values of each image are
  extracted. Then \code{\link{pairs.default}} is called to 
  plot the scatterplot matrix.

  Any arguments in \code{\dots} which are not pixel images will be
  passed to \code{\link{pairs.default}} to control the plot.

  The return value of \code{pairs.im} is a data frame, returned
  invisibly. The data frame has one column
  for each image. Each row contains the pixel values of
  the different images for one pixel in the raster.
  If \code{drop=TRUE} (the default), any row which contains \code{NA}
  is deleted. The plot is not affected by the value of \code{drop}.
}
\section{Image or Contour Plots}{
  Since the scatterplots may show very dense concentrations of points,
  it may be useful to set \code{panel=panel.image}
  or \code{panel=panel.contour} to draw a colour image or contour plot
  of the kernel-smoothed density of the scatterplot in each panel.
  The argument \code{panel} is passed
  to \code{\link{pairs.default}}. See the help for 
  \code{\link{panel.image}} and \code{\link{panel.contour}}.
}
\section{Low Level Control of Graphics}{
  To control the appearance of the individual scatterplot panels,
  see \code{\link{pairs.default}}, \code{\link{points}}
  or \code{\link{par}}.
  To control the plotting symbol for the points in the scatterplot,
  use the arguments \code{pch}, \code{col}, \code{bg} as described
  under \code{\link{points}} (because 
  the default panel plotter is the function \code{\link{points}}).
  To suppress the tick marks on the plot axes,
  type \code{par(xaxt="n", yaxt="n")} before calling \code{pairs}.
}
\value{
  Invisible. A \code{data.frame} containing the
  corresponding pixel values for each image.
  The return value also belongs to the class \code{plotpairsim} which has
  a plot method, so that it can be re-plotted.
}
\seealso{
  \code{\link{pairs}},
  \code{\link{pairs.default}},
  \code{\link{panel.contour}},
  \code{\link{panel.image}},
  \code{\link{plot.im}},
  \code{\link{cov.im}},
  \code{\link{im}},
  \code{\link{par}}
}
\examples{
  X <- density(rpoispp(30))
  Y <- density(rpoispp(40))
  Z <- density(rpoispp(30))
  p <- pairs(X,Y,Z)
  p
  plot(p)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}

