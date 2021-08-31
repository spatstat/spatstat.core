\name{collapse.fv}
\alias{collapse.fv}
\alias{collapse.anylist}
\title{
  Collapse Several Function Tables into One
}
\description{
  Combines several function tables (objects of class \code{"fv"})
  into a single function table, merging columns that are identical
  and relabelling columns that are different.
}
\usage{
\method{collapse}{fv}(object, \dots, same = NULL, different = NULL)

\method{collapse}{anylist}(object, \dots, same = NULL, different = NULL)
}
\arguments{
  \item{object}{
    An object of class \code{"fv"}, or a list of such objects.
  }
  \item{\dots}{
    Additional objects of class \code{"fv"}.
  }
  \item{same}{
    Character string or character vector specifying a column or columns
    of function values that are identical in different \code{"fv"}
    objects. These columns will be included only once in the result.
  }
  \item{different}{
    Character string or character vector specifying a column or columns
    of function values, that are different in different \code{"fv"} objects.
    Each of these columns of data will be included, with
    labels that distinguish them from each other.
  }
}
\details{
  This is a method for the generic function
  \code{\link[nlme]{collapse}}.
  
  It combines the data in several function tables 
  (objects of class \code{"fv"}, see \code{\link{fv.object}})
  to make a single function table.
  It is essentially a smart wrapper for
  \code{\link{cbind.fv}}.

  A typical application is to calculate the same summary statistic
  (such as the \eqn{K} function) for different point patterns,
  and then to use \code{collapse.fv} to combine the results into a
  single object that can easily be plotted. See the Examples.
  
  The arguments \code{object} and \code{\dots} should be function tables
  (objects of class \code{"fv"}, see \code{\link{fv.object}})
  that are compatible in the sense that they
  have the same values of the function argument.
  (This can be ensured by applying \code{\link{harmonise.fv}} to them.)

  The argument \code{same} identifies any columns that are present
  in some or all of the function tables, and 
  which are known to contain exactly the same values in each table
  that includes them. This column or columns will be
  included only once in the result.

  The argument \code{different} identifies any columns that are present
  in some or all of the function tables, and which may contain
  different numerical values in different tables.
  Each of these columns will be included, with labels
  to distinguish them. 

  Columns that are not named in \code{same} or \code{different} will not
  be included.

  The function argument is always included and does not need to be
  specified.
  
  The arguments \code{same} and \code{different} can be \code{NULL},
  or they can be character vectors containing the names
  of columns of \code{object}. The argument \code{different} can be
  one of the abbreviations recognised by \code{\link{fvnames}}.
}
\value{
  Object of class \code{"fv"}.
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{cbind.fv}}
}
\examples{
  # generate simulated data
  X <- replicate(3, rpoispp(100), simplify=FALSE)
  names(X) <- paste("Simulation", 1:3)
  # compute K function estimates
  Klist <- anylapply(X, Kest)
  # collapse
  K <- collapse(Klist, same="theo", different="iso")
  K
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}
