\name{intensity.interact}  
\alias{intensity.interact}
\title{
  Intensity of pairwise-interaction models
}
\description{
  Approximates the intensity of a Gibbs pairwise-interaction model.
}
\usage{
  \method{intensity}{interact}(X, beta=NULL,interaction.coef=NULL,method=c("Poisson","DPP"), \dots)
}
\arguments{
  \item{X}{
    A pairwise-interaction function (object of class \code{"interact"}).
  }
  \item{beta}{
    The spatial trend. A numerical value or a pixel image. 
  }
  \item{interaction.coef}{
    Vector of interaction coefficients. 
  }
  \item{method}{
    The approximation method. By defaut the Poisson-saddlepoint approximation is used.
  }
  \item{\dots}{
    Arguments passed to \code{\link{predict.ppm}} in some cases.
    See Details.
  }
}


\details{
  This is a method for the generic function \code{\link{intensity}}
  for pairwise-interaction functions (class \code{"interact"}).
  
  For a pairwise-interaction function, each pair of points in the point pattern contributes a factor
  \eqn{h(d)} to the probability density, where \eqn{d} is the distance
  between the two points. The factor term \eqn{h(d)} is
  \deqn{h(d) = \exp(-\theta \mbox{pot}(d))}{h(d) = exp(-theta * pot(d))}
  provided \eqn{\mbox{pot}(d)}{pot(d)} is finite, 
  where \eqn{\theta}{theta} is the coefficient vector in the model. 
  The probability density of a point pattern is
  \deqn{
    f(x_1,\ldots,x_n) =
      \alpha \beta^{n(x)} \prod_{i} h(d_i)
  }{
    f(x_1,\ldots,x_n) =
      alpha . beta^n(x) . product{ h(d[i])}
  }
  where \eqn{x_1,\ldots,x_n}{x[1],\ldots,x[n]} represent the 
  points of the pattern, \eqn{n(x)} is the number of points in the
  pattern, the product is over all pairwise distances \eqn{d_i}{d[i]},
  and \eqn{\alpha}{alpha} is the normalising constant.
  
  The function \eqn{h(d)} is given by the interaction object \code{X}, \eqn{\theta}{theta} corresponds to \code{interaction.coef} and \eqn{\beta}{beta} is the spatial trend. 
  
  
  The intensity of a point process model is the expected
  number of random points per unit area. 
  
  If \code{X} is a Poisson point process model, \code{beta} is returned. 
  
  If \code{X} is another pairwise-interaction Gibbs point process model and \code{method="Poisson"}, the intensity is
  computed approximately using the Poisson-saddlepoint approximation
  (Baddeley and Nair, 2012a, 2012b).
   If \code{method="DPP"}, the approximation is based on (Coeurjolly and Lavancier, 2018), resulting to a better 
  approximation of the intensity for repulsive interactions. This method is only available for stationary models (\code{beta} is a numerical value) and requires both the interaction range and the hardcore parameter. If not supplied, the interaction range is approximated and the hardcore parameter is set to 0. 
 
  For a non-stationary Gibbs model (\code{beta} is a pixel image), the 
  pseudostationary solution (Baddeley and Nair, 2012b;
                             Anderssen et al, 2014) is used. The result is a pixel image,
  whose resolution is controlled by the arguments \code{\dots} which are passed
  to \code{\link{predict.ppm}}.
}
\value{
  A numeric value (if the model is stationary)
  or a pixel image.
}
\references{
  Anderssen, R.S., Baddeley, A., DeHoog, F.R. and Nair, G.M. (2014)
  Solution of an integral equation arising in spatial point process theory.
  \emph{Journal of Integral Equations and Applications} 
  \bold{26} (4) 437--453.
  
  Baddeley, A. and Nair, G. (2012a) 
  Fast approximation of the intensity of Gibbs point processes.
  \emph{Electronic Journal of Statistics} \bold{6} 1155--1169.
  
  Baddeley, A. and Nair, G. (2012b)
  Approximating the moments of a spatial point process.
  \emph{Stat} \bold{1}, 1, 18--30.
  \code{DOI: 10.1002/sta4.5}
  
  
  Coeurjolly, J.-F. and Lavancier, F. (2018)
  Approximation intensity for pairwise interaction Gibbs point processes
  using determinantal point processes.
  \emph{Electronic Journal of Statistics} \bold{12} 3181--3203.
  
}
\seealso{
  \code{\link{intensity}},
  \code{\link{intensity.ppp}}
  \code{\link{intensity.ppm}}
}
\examples{
  inte<-Strauss(r=0.1)
  beta<-1000
  co<- log(0.5)
  intensity(inte,beta,co,method="Poisson")
  intensity(inte,beta,co,method="DPP")
  
  inte<-StraussHard(r=0.1, hc=0.05)
  beta<-1000
  co<-log(0.5)
  intensity(inte,beta,co,method="Poisson")
  intensity(inte,beta,co,method="DPP")
  
}
\author{
  J.-F. Coeurjolly and F. Lavancier. 
}
\keyword{spatial}
\keyword{models}