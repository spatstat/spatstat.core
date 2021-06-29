\name{cov.im}
\alias{cov.im}
\alias{cor.im}
\title{
  Covariance and Correlation between Images
}
\description{
  Compute the covariance or correlation between (the corresponding pixel
  values in) several images.
}
\usage{
cov.im(\dots, use = "everything", method = c("pearson", "kendall", "spearman"))
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is 
    a pixel image (object of class \code{"im"}).
    Alternatively, a single argument which is a list of pixel images.
  }
  \item{use}{
    Argument passed to \code{\link[stats]{cov}}
    or \code{\link[stats]{cor}} determining how to handle
    \code{NA} values in the data.
  }
  \item{method}{
    Argument passed to \code{\link[stats]{cov}}
    or \code{\link[stats]{cor}} determining the type of correlation
    that will be computed.
  }
}
\details{
  The arguments \code{\dots}
  should be pixel images (objects of class \code{"im"}).
  Their spatial domains must overlap, but need not have the same
  pixel dimensions.

  These functions compute the covariance or correlation between
  the corresponding pixel values in the images given.

  The pixel image domains are intersected, and converted to a
  common pixel resolution. Then the corresponding pixel values of each image are
  extracted. Finally the correlation or covariance between the
  pixel values of each pair of images, at corresponding pixels, is
  computed.

  The result is a symmetric matrix with one row and column for each image.
  The \code{[i,j]} entry is the correlation or covariance between
  the \code{i}th and \code{j}th images in the argument list.
  The row names and column names of the matrix are copied from the
  argument names if they were given (i.e. if the arguments were given
  as \code{name=value}).

  Note that \code{\link[stats]{cor}} and \code{\link[stats]{cov}}
  are not generic, so you have to type \code{cor.im}, \code{cov.im}.
}
\value{
  A symmetric matrix.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[stats]{cor}},
  \code{\link[stats]{cov}}

  \code{\link{pairs.im}}
}
\examples{
  cor.im(bei.extra)
}
\keyword{spatial}
\keyword{univar}
\keyword{nonparametric}

