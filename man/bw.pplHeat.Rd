\name{bw.pplHeat}
\alias{bw.pplHeat}
\title{
  Bandwidth Selection for Diffusion Smoother by Likelihood Cross-Validation
}
\description{
  Selects an optimal bandwidth for diffusion smoothing
  by point process likelihood cross-validation.
}
\usage{
bw.pplHeat(X, \dots, srange=NULL, ns=16, sigma=NULL,
         leaveoneout=TRUE, verbose = TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{densityHeat.ppp}}.
  }
  \item{srange}{
    Numeric vector of length 2 specifying a range of bandwidths to be
    considered.
  }
  \item{ns}{
    Integer. Number of candidate bandwidths to be considered.
  }
  \item{sigma}{
    Maximum smoothing bandwidth.
    A numeric value, or a pixel image, or a \code{function(x,y)}.
    Alternatively a numeric vector containing a sequence of
    candidate bandwidths.
  }
  \item{leaveoneout}{
    Logical value specifying whether intensity values at data points
    should be estimated using the leave-one-out rule.
  }
  \item{verbose}{
    Logical value specifying whether to print progress reports.
  }
}
\details{
  This algorithm selects the optimal global bandwidth for
  kernel estimation of intensity for the dataset \code{X}
  using diffusion smoothing \code{\link{densityHeat.ppp}}.
  
  If \code{sigma} is a numeric value, the algorithm finds the
  optimal bandwidth \code{tau <= sigma}.

  If \code{sigma} is a pixel image or function, the algorithm
  finds the optimal fraction \code{0 < f <= 1} such that
  smoothing with \code{f * sigma} would be optimal.
}
\value{
  A numerical value giving the selected bandwidth
  (if \code{sigma} was a numeric value)
  or the selected fraction of the maximum bandwidth
  (if \code{sigma} was a pixel image or function).
  The result also belongs to the class \code{"bw.optim"} which can be
  plotted. 
}
\author{
  Adrian Baddeley and Tilman Davies.
}
\seealso{
  \code{\link{bw.CvLHeat}} for an alternative method.
  
  \code{\link{densityHeat.ppp}}
}
\examples{
  online <- interactive()
  if(!online) op <- spatstat.options(npixel=32)
  f <- function(x,y) { dnorm(x, 2.3, 0.1) * dnorm(y, 2.0, 0.2) }
  X <- rpoint(15, f, win=letterR)
  plot(X)
  b <- bw.pplHeat(X, sigma=0.25)
  b
  plot(b)
  if(!online) spatstat.options(op)
}
\keyword{spatial}
\keyword{smooth}
