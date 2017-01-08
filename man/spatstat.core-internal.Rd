\name{spatstat.core-internal} 
\title{Internal spatstat.core functions}
\alias{existsSpatstatVariable}
\alias{getSpatstatVariable}
\alias{putSpatstatVariable}
\alias{warn.once}
\alias{waxlyrical}
\alias{ruletextline}
\alias{parbreak}
%% colourtables.R
\alias{print.colourmap}
\alias{print.lut}
\alias{print.summary.lut}
\alias{summary.lut}
\alias{colouroutputs}
\alias{colouroutputs<-}
%% colourtools.R
\alias{is.col.argname}
\alias{col.args.to.grey} 
\alias{rgbNA}
\alias{hsvNA}
\alias{do.call.plotfun}
\description{
  Internal spatstat.core functions.
}
\usage{
  existsSpatstatVariable(name)
  getSpatstatVariable(name)
  putSpatstatVariable(name, value)
  warn.once(key, \dots)
  waxlyrical(type, terse)
  ruletextline(ch, n, terse)
  parbreak(terse)
%% colourtables.R
  \method{print}{colourmap}(x, \dots)
  \method{print}{lut}(x, \dots)
  \method{print}{summary.lut}(x, \dots)
  \method{summary}{lut}(object, \dots)
  colouroutputs(x)
  colouroutputs(x) <- value
%% colourtools.R
  is.col.argname(x)
  col.args.to.grey(x, \dots)
  rgbNA(red, green, blue, alpha, maxColorValue)
  hsvNA(h, s, v, alpha)
  do.call.plotfun(fun, arglist, \dots)
}
\details{
  These internal \pkg{spatstat.core} functions are not usually called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.core} to the next.
}
\keyword{internal}
