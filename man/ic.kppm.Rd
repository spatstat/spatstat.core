\name{ic.kppm}
\alias{ic.kppm}

\concept{point process model}
\concept{Cox point process}
\concept{cluster process}
\concept{Neyman-Scott cluster process}
\title{Model selection criteria for the intensity function of a Poisson, cluster or Cox point process}
\description{
  Information criteria for selecting the intensity function model of a Poisson, cluster or Cox point process.
}
\usage{
  kppm(K,X)
}
\arguments{
  \item{K}{
    An object returned by ppm representing a fitted Poisson point process model or an object returned by kppm representing a fitted cluster or Cox point process model.
    \item{X}
    The point pattern data (just used to provide number of points).
}
\details{
  This function returns information criteria for selecting the intensity function model of a Poisson, cluster or Cox point process fitted by first order composite likelihood (i.e. using the Poisson likelihood function). 
  
  Degrees of freedom \eqn{df} for the information criteria are given by the trace of \eqn{S^{-1} \Sigma} where \eqn{S} is the sensitity matrix and \eqn{\Sigma} is the variance matrix for the log composite likelihood score function. In case of a Poisson process, \eqn{df} is the number of parameters in the model for the intensity function.
  
  The composite Bayesian information criterion (cbic) is \eqn{-2l + \log(n) df} where \eqn{l} is the maximal log composite likelihood and \eqn{n} is the observed number of points. It reduces to the BIC criterion in case of a Poisson process.
  
  The composite information criterion (cic) is \eqn{-2l + 2 df} and reduces to the AIC in case of a Poisson process.
  
  NOTE: the information criteria are for selecting the intensity function model (a set of covariates) within a given model class. It can not be used to choose among different types cluster or Cox point process models (e.g. can not choose between Thomas and LGCP).
}
\value{
  A list with entries loglike, cbic, cic and df. Here loglike is the fitted log first order composite likelihood, cbic is composite Bayesian information criterion, cic is is the composite likelihood criterion and df is the adjusted degrees of freedom for the fitted intensity function model.
 } 
\seealso{
  \code{\link{kppm}},
  \code{\link{kppm.object}},
 
}
\references{
  Choiruddin, A., Coeurjolly, J.F. and Waagepetersen, R. (2020) 
  Information criteria for inhomogeneous spatial point processes.
  \emph{Australian & New Zealand Journal of Statistics}
  To appear. 
 
}
\examples{

  #model with one covariate
  fit1=kppm(bei~elev,data=bei.extra)
  ic1=ic.kppm(fit1,bei)
  

  #model with two covariates
  fit2=kppm(bei~elev+grad,data=bei.extra)
  ic2=ic.kppm(fit2,bei)
  
  #smallest cbic for fit1 but smallest cic for fit2
}
\author{
  Achmad Choiruddin, Jean-Francois Coeurjolly and Rasmus Waagepetersen.
}
\keyword{spatial}
\keyword{models}

