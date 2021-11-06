\name{spatcov}
\alias{spatcov}
\title{
  Estimate the Spatial Covariance Function of a Random Field
}
\description{
  Given a pixel image, calculate an estimate of the spatial
  covariance function. Given two pixel images, calculate an estimate of their
  spatial cross-covariance function.
}
\usage{
spatcov(X, Y=X, \dots, correlation=FALSE, isotropic = TRUE,
        clip = TRUE, pooling=TRUE)
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"}).
  }
  \item{Y}{
    Optional. Another pixel image.
  }
  \item{correlation}{
    Logical value specifying whether to standardise
    so that the spatial correlation function is returned.
  }
  \item{isotropic}{
    Logical value specifying whether to assume the covariance
    is isotropic, so that the result is a function of the lag distance.
  }
  \item{clip}{
    Logical value specifying whether to restrict the results to
    the range of spatial lags where the estimate is reliable.
  }
  \item{pooling}{
    Logical value specifying the estimation method when \code{isotropic=TRUE}.
  }
  \item{\dots}{Ignored.}
}
\details{
  In normal usage, only the first argument \code{X} is given.
  Then the pixel image \code{X} is treated as a realisation of a stationary
  random field, and its spatial covariance function is estimated.

  Alternatively if \code{Y} is given,
  then \code{X} and \code{Y} are assumed to be
  jointly stationary random fields, and their spatial cross-covariance
  function is estimated.

  For any random field \code{X}, the spatial covariance
  is defined for any two spatial locations \eqn{u} and \eqn{v} by
  \deqn{
    C(u,v) = \mbox{cov}(X(u), X(v))
  }{
    C(u,v) = cov(X(u), X(v))
  }
  where \eqn{X(u)} and \eqn{X(v)} are the values of the random field
  at those locations. Here\eqn{\mbox{cov}}{cov} denotes the
  statistical covariance, defined for any random variables
  \eqn{A} and \eqn{B} by
  \eqn{\mbox{cov}(A,B) = E(AB) - E(A) E(B)}{cov(A,B) = E(AB) - E(A) E(B)}
  where \eqn{E(A)} denotes the expected value of \eqn{A}.

  If the random field is assumed to be stationary (at least second-order
  stationary) then the spatial covariance \eqn{C(u,v)}
  depends only on the lag vector \eqn{v-u}:
  \deqn{
    C(u,v) = C_2(v-u)
  }
  \deqn{
    C(u,v) = C2(v-u)
  }
  where \eqn{C_2}{C2} is a function of a single vector argument.
  
  If the random field is stationary and isotropic, then the spatial
  covariance depends only on the lag distance
  \eqn{\| v - u \|}{||v-u||}:
  \deqn{
    C_2(v-u) = C_1(\|v-u\|)
  }{
    C2(v-u) = C1(||v-u||)
  }
  where \eqn{C_1}{C1} is a function of distance.
  
  The function \code{spatcov} computes estimates of the
  covariance function \eqn{C_1}{C1} or \eqn{C_2}{C2} as follows:
  \itemize{
    \item 
    If \code{isotropic=FALSE}, an estimate of the 
    covariance function \eqn{C_2}{C2} is computed,
    assuming the random field is stationary, using the naive
    moment estimator,
    \code{C2 = imcov(X-mean(X))/setcov(Window(X))}.
    The result is a pixel image.
    \item 
    If \code{isotropic=TRUE} (the default)
    an estimate of the covariance function \eqn{C_1}{C1}
    is computed, assuming the random field is stationary and isotropic.
    \itemize{
      \item 
      When \code{pooling=FALSE}, the estimate of \eqn{C_1}{C1}
      is the rotational average of the naive estimate of \eqn{C_2}{C2}.
      \item 
      When \code{pooling=TRUE} (the default), the  estimate of \eqn{C_1}{C1}
      is the ratio of the rotational averages of the numerator and
      denominator which form the naive estimate of  \eqn{C_2}{C2}.
    }
    The result is a function object (class \code{"fv"}).
  }

  If the argument \code{Y} is given, it should be a pixel image
  compatible with \code{X}. An estimate of the spatial cross-covariance function
  between \code{X} and \code{Y} will be computed.
}
\value{
  If \code{isotropic=TRUE} (the default), the result is a function value
  table (object of class \code{"fv"}) giving the estimated values of the
  covariance function or spatial correlation function
  for a sequence of values of the spatial lag
  distance \code{r}.

  If \code{isotropic=FALSE}, the result is a pixel image
  (object of class \code{"im"}) giving the estimated values of the
  spatial covariance function or spatial correlation function
  for a grid of values of the spatial lag vector.
}
\author{
  \adrian
}
\seealso{
  \code{\link{imcov}}, \code{\link{setcov}}
}
\examples{
if(offline <- !interactive()) op <- spatstat.options(npixel=32)

  D <- density(cells)
  plot(spatcov(D))

if(offline) spatstat.options(op)
}
\keyword{spatial}
\keyword{nonparametric}
