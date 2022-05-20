\name{spatstat.core-deprecated}
\alias{evalCovar}
\alias{which.max.im}
\title{Deprecated spatstat.core functions}
\description{
  Deprecated spatstat.core functions.
}
\usage{
evalCovar(model, covariate, \dots)
which.max.im(x)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat.core} package.
  

  \code{which.max.im(x)} is replaced by
  \code{\link{im.apply}(x, which.max)}.
  
  The internal function \code{evalCovar}
  is replaced by the internal function \code{spatialCovariateEvidence}.
  
}
\value{
\code{which.max.im} returns an integer.
}
\keyword{internal}
