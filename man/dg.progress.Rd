\name{dg.progress}
\alias{dg.progress}
\title{
   Progress Plot of Dao-Genton Test of Spatial Pattern
}
\description{
  Generates a progress plot (envelope representation) of the
  Dao-Genton test for a spatial point pattern.
}
\usage{
dg.progress(X, fun = Lest, \dots,
            exponent = 2, nsim = 19, nsimsub = nsim - 1,
            nrank = 1, alpha, leaveout=1, interpolate = FALSE, rmin=0,
            savefuns = FALSE, savepatterns = FALSE, verbose=TRUE)
}
\arguments{
  \item{X}{
    Either a point pattern (object of class \code{"ppp"}, \code{"lpp"}
    or other class), a fitted point process model (object of class \code{"ppm"},
    \code{"kppm"} or other class) or an envelope object (class
    \code{"envelope"}). 
  }
  \item{fun}{
    Function that computes the desired summary statistic
    for a point pattern. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{envelope}}.
    Useful arguments include \code{alternative} to
    specify one-sided or two-sided envelopes.
  }  
  \item{exponent}{
    Positive number. The exponent of the \eqn{L^p} distance.
    See Details.
  }
  \item{nsim}{
    Number of repetitions of the basic test.
  }
  \item{nsimsub}{
    Number of simulations in each basic test. There will be \code{nsim}
    repetitions of the basic test, each involving \code{nsimsub} simulated
    realisations, so there will be a total
    of \code{nsim * (nsimsub + 1)} simulations.
  }
  \item{nrank}{
    Integer. The rank of the critical value of the Monte Carlo test,
    amongst the \code{nsim} simulated values.
    A rank of 1 means that the minimum and maximum
    simulated values will become the critical values for the test.
  }
  \item{alpha}{
    Optional. The significance level of the test.
    Equivalent to \code{nrank/(nsim+1)} where \code{nsim} is the
    number of simulations.
  }
  \item{leaveout}{
    Optional integer 0, 1 or 2 indicating how to calculate the
    deviation between the observed summary function and the
    nominal reference value, when the reference value must be estimated
    by simulation. See Details.
  }
  \item{interpolate}{
    Logical value indicating how to compute the critical value.
    If \code{interpolate=FALSE} (the default), a standard Monte Carlo test
    is performed, and the critical value is the largest
    simulated value of the test statistic (if \code{nrank=1})
    or the \code{nrank}-th largest (if \code{nrank} is another number).
    If \code{interpolate=TRUE}, kernel density estimation
    is applied to the simulated values, and the critical value is
    the upper \code{alpha} quantile of this estimated distribution.
  }
  \item{rmin}{
    Optional. Left endpoint for the interval of \eqn{r} values
    on which the test statistic is calculated.
  }
  \item{savefuns}{
    Logical value indicating whether to save the simulated
    function values (from the first stage).
  }
  \item{savepatterns}{
    Logical value indicating whether to save the simulated
    point patterns (from the first stage). 
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports.
  }
}
\details{
  The Dao and Genton (2014) test for a spatial point pattern
  is described in \code{\link{dg.test}}.
  This test depends on the choice of an interval of
  distance values (the argument \code{rinterval}).
  A \emph{progress plot} or \emph{envelope representation}
  of the test (Baddeley et al, 2014, 2015; Baddeley, Rubak and Turner, 2015) is a plot of the
  test statistic (and the corresponding critical value) against the length of
  the interval \code{rinterval}.
  
  The command \code{dg.progress} effectively performs 
  \code{\link{dg.test}} on \code{X} using all possible intervals
  of the form \eqn{[0,R]}, and returns the resulting values of the test
  statistic, and the corresponding critical values of the test,
  as a function of \eqn{R}. 

  The result is an object of class \code{"fv"}
  that can be plotted to obtain the progress plot. The display shows
  the test statistic (solid black line) and the test
  acceptance region (grey shading).
  If \code{X} is an envelope object, then some of the data stored
  in \code{X} may be re-used:
  \itemize{
    \item
    If \code{X} is an envelope object containing simulated functions,
    and \code{fun=NULL}, then
    the code will re-use the simulated functions stored in \code{X}.
    \item
    If \code{X} is an envelope object containing
    simulated point patterns, 
    then \code{fun} will be applied to the stored point patterns
    to obtain the simulated functions.
    If \code{fun} is not specified, it defaults to \code{\link{Lest}}.
    \item
    Otherwise, new simulations will be performed,
    and \code{fun} defaults to  \code{\link{Lest}}.
  }
  If the argument \code{rmin} is given, it specifies the left endpoint
  of the interval defining the test statistic: the tests are
  performed using intervals \eqn{[r_{\mbox{\scriptsize min}},R]}{[rmin,R]}
  where \eqn{R \ge r_{\mbox{\scriptsize min}}}{R \ge rmin}.

  The argument \code{leaveout} specifies how to calculate the
  discrepancy between the summary function for the data and the
  nominal reference value, when the reference value must be estimated
  by simulation. The values \code{leaveout=0} and
  \code{leaveout=1} are both algebraically equivalent (Baddeley et al, 2014,
  Appendix) to computing the difference \code{observed - reference}
  where the \code{reference} is the mean of simulated values.
  The value \code{leaveout=2} gives the leave-two-out discrepancy
  proposed by Dao and Genton (2014).
}
\value{
  An object of class \code{"fv"} that can be plotted to
  obtain the progress plot. 
}
\references{
  Baddeley, A., Diggle, P., Hardegen, A., Lawrence, T.,
  Milne, R. and Nair, G. (2014)
  On tests of spatial pattern based on simulation envelopes.
  \emph{Ecological Monographs} \bold{84} (3) 477--489.

  Baddeley, A., Hardegen, A., Lawrence, L., 
  Milne, R.K., Nair, G.M. and Rakshit, S. (2015)
  Pushing the envelope: extensions of graphical
  Monte Carlo tests. Unpublished manuscript.

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press, Boca Raton, FL.

  Dao, N.A. and Genton, M. (2014)
  A Monte Carlo adjusted goodness-of-fit test for
  parametric models describing spatial point patterns.
  \emph{Journal of Graphical and Computational Statistics}
  \bold{23}, 497--517.
}
\author{
  Adrian Baddeley, Andrew Hardegen, Tom Lawrence,
  Robin Milne, Gopalan Nair and Suman Rakshit.
  Implemented by \spatstatAuthors.
}
\seealso{
  \code{\link{dg.test}},
  \code{\link{dclf.progress}}
}
\examples{
   ns <- if(interactive()) 19 else 5
   plot(dg.progress(cells, nsim=ns))
}
\keyword{spatial}
\keyword{htest}

