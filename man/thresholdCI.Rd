\name{thresholdCI}
\alias{thresholdCI}
\title{
  Confidence Interval for Threshold of Numerical Predictor
}
\description{
  Given a point pattern and a spatial covariate that has some predictive
  value for the point pattern,
  compute a confidence interval for the optimal value of the threshold
  that should be used to convert the covariate to a binary predictor.
}
\usage{
thresholdCI(X, Z, confidence = 0.95, nsim = 1000, parametric = FALSE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{Z}{
    Spatial covariate with numerical values.
    Either a pixel image (object of class \code{"im"}),
    a distance function (object of class \code{"distfun"})
    or a \code{function(x,y)} in the \R language.
  }
  \item{confidence}{
    Confidence level. A number between 0 and 1.
  }
  \item{nsim}{
    Number of bootstrap simulations to perform.
  }
  \item{parametric}{
    Logical value specifying whether to use the parametric bootstrap.
  }
}
\details{
  The spatial covariate \code{Z} is assumed to have some utility as a
  predictor of the point pattern \code{X}.
  
  This code computes a bootstrap confidence interval
  for the best threshold value \eqn{z} for converting the
  numerical predictor to a binary predictor, for use in
  techniques such as Weights of Evidence.
}
\value{
  A matrix containing upper and lower limits for the
  threshold \code{z} and the corresponding upper and lower limits for
  the fraction of area of the study region.
}
\references{
Baddeley, A., Brown, W., Milne, R.K., Nair, G.,
Rakshit, S., Lawrence, T., Phatak, A. and Fu, S.C. (2021)
Optimal thresholding of predictors in mineral prospectivity analysis.
\emph{Natural Resources Research} \bold{30} 923--969.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{thresholdSelect}}
}
\examples{
  gold <- rescale(murchison$gold, 1000, "km")
  faults <- rescale(murchison$faults, 1000, "km")
  distfault <- distfun(faults)
  thresholdCI(gold, distfault, nsim=100)
}
\keyword{spatial}
\keyword{models}
\keyword{nonparametric}
