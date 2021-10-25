\name{thresholdSelect}
\alias{thresholdSelect}
\title{
  Select Threshold to Convert Numerical Predictor to Binary Predictor
}
\description{
  Given a point pattern and a spatial covariate that has some predictive
  value for the point pattern,
  determine the optimal value of the threshold for converting
  the covariate to a binary predictor.
}
\usage{
thresholdSelect(X, Z, method = c("Y", "LL", "AR", "t", "C"), Zname)
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
  \item{method}{
    Character string (partially matched)
    specifying the method to be used to select the
    optimal threshold value. See Details.
  }
  \item{Zname}{
    Optional character string giving a short name for the covariate.
  }
}
\details{
  The spatial covariate \code{Z} is assumed to have some utility as a
  predictor of the point pattern \code{X}.
  
  This code chooses the best threshold value \eqn{v} for converting the
  numerical predictor \code{Z} to a binary predictor, for use in
  techniques such as Weights of Evidence.

  The best threshold is selected by maximising the criterion
  specified by the argument \code{method}. Options are:
  \itemize{
    \item \code{method="Y"} (the default): the Youden criterion 
    \item \code{method="LL"}: log-likelihood 
    \item \code{method="AR"}: the Akman-Raftery criterion
    \item \code{method="t"}: the Studentised Weights-of-Evidence contrast
    \item \code{method="C"}: the Weights-of-Evidence contrast
  }
  These criteria are explained in Baddeley et al (2021).
}
\value{
  A numerical value giving the selected threshold.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted (the plot shows the criterion used to select
  the threshold).
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
  \code{\link{thresholdCI}}
}
\examples{
  gold <- rescale(murchison$gold, 1000, "km")
  faults <- rescale(murchison$faults, 1000, "km")
  distfault <- distfun(faults)
  z <- thresholdSelect(gold, distfault)
  z
  plot(z, xlim=c(0, 20))
}
\keyword{spatial}
\keyword{models}
\keyword{nonparametric}
