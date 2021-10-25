#'
#'   hotbox.R
#'
#'  Heat kernel for a one-dimensional rod
#'  and two-dimensional rectangle
#'
#'  Code from Greg McSwiggan and Adrian Baddeley

hotbox <- function(Xsource, Xquery, sigma,
                   ..., W=NULL, squared=FALSE, nmax=20) {
  #' heat kernel in a rectangle
  check.1.real(sigma)
  if(is.null(W)) {
    if(is.ppp(Xsource)) W <- Window(Xsource) else
    if(is.sob(Xquery)) W <- Window(Xquery) else
    stop("No window information is present")
  } else {
    stopifnot(is.owin(W))
    if(!is.sob(Xsource)) Xsource <- as.ppp(Xsource, W)
    if(!is.sob(Xquery)) Xquery <- as.ppp(Xquery, W)
  }
  if(!is.rectangle(W))
    stop("The window must be a rectangle")
  slen <- sidelengths(W)
  Xsource <- shift(Xsource, origin="bottomleft")
  Xquery <- shift(Xquery, origin="bottomleft")
  nsource <- npoints(Xsource)
  if(is.ppp(Xquery)) {
    nquery <- npoints(Xquery)
    answer <- numeric(nquery)
    for(i in seq_len(nsource)) {
      cx <- hotrod(slen[1], Xsource$x[i], Xquery$x, sigma,
                   ends="insulated", nmax=nmax)
      cy <- hotrod(slen[2], Xsource$y[i], Xquery$y, sigma,
                   ends="insulated", nmax=nmax)
      contrib <- cx * cy
      if(squared) contrib <- contrib^2
      answer <- answer + contrib
    }
  } else if(is.im(Xquery) || is.owin(Xquery)) {
    Xquery <- as.im(Xquery, ...)
    if(anyNA(Xquery)) stop("Image must be a full rectangle")
    ansmat <- matrix(0, nrow(Xquery), ncol(Xquery))
    xx <- Xquery$xcol
    yy <- Xquery$yrow
    for(i in seq_len(nsource)) {
      cx <- hotrod(slen[1], Xsource$x[i], xx, sigma,
                   ends="insulated", nmax=nmax)
      cy <- hotrod(slen[2], Xsource$y[i], yy, sigma,
                   ends="insulated", nmax=nmax)
      contrib <- outer(cy, cx, "*")
      if(squared) contrib <- contrib^2
      ansmat <- ansmat + contrib
    }
    answer <- Xquery
    answer[] <- ansmat
  } else stop("Unrecognised format for Xquery")
  return(answer)
}
