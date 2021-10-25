\name{pairMean}
\alias{pairMean}
\title{
  Mean of a Function of Interpoint Distance
}
\description{
  Computes the mean value, or the double integral, of a specified
  function of the distance between two independent random points
  in a given window or windows.
}
\usage{
pairMean(fun, W, V = NULL, ..., normalise = TRUE)
}
\arguments{
  \item{fun}{
    A function in the \R language which takes one argument.
  }
  \item{W}{
    A window (object of class \code{"owin"}) containing the
    first random point.
  }
  \item{V}{
    Optional. Another window containing the second random point.
    Defaults to \code{W}.
  }
  \item{\dots}{
    Further optional arguments passed to \code{\link{distcdf}} to determine the
    pixel resolution for the calculation and the
    probability distributions of the random points.
  }
  \item{normalise}{
    Logical value specifying whether to calculate the mean
    value (\code{normalise=TRUE}, the default) or the
    double integral (\code{normalise=FALSE}).
  }
}
\details{
  This command computes the mean value of
  \code{fun(T)} where \code{T} is 
  the Euclidean distance \eqn{T = \|X_1 - X_2\|}{T = |X1-X2|}
  between two independent random points \eqn{X_1}{X1} and \eqn{X_2}{X2}.

  In the simplest case, the command \code{pairMean(fun, W)},
  the random points are assumed to be uniformly distributed in the same
  window \code{W}.  Alternatively the two random points may be 
  uniformly distributed in two different windows \code{W} and \code{V}.
  Other options are described in \code{\link{distcdf}}.

  The algorithm uses \code{\link{distcdf}} to compute the cumulative
  distribution function of \code{T}, and \code{\link{stieltjes}} to
  compute the mean value of \code{fun(T)}.

  If \code{normalise=TRUE} (the default) the result is the mean value of
  \code{fun(T)}. If \code{normalise=FALSE} the result is the double
  integral.
}
\value{
  A single numeric value.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{distcdf}}
}
\examples{
   pairMean(function(d) { d^2 }, disc())
}
\keyword{spatial}
\keyword{math}
