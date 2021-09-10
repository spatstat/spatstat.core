\name{spatstat.core-deprecated}
\alias{which.max.im}
\title{Deprecated spatstat.core functions}
\description{
  Deprecated spatstat.core functions.
}
\usage{
which.max.im(x)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat.core} package.
  
  \code{as.psp.owin} has been replaced by \code{\link{edges}}.
  
  \code{which.max.im(x)} is replaced by
  \code{\link{im.apply}(x, which.max)}.

  \code{circumradius} is replaced by the more appropriately named
  \code{boundingradius}.

  \code{rjitterlpp} is replaced by \code{rjitter.lpp}, a method for the
  generic \code{rjitter}.
}
\value{
\code{which.max.im} returns an integer.
}
\keyword{internal}
