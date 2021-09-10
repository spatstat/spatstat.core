#
# randomseg.R
#
# $Revision: 1.16 $ $Date: 2021/09/10 08:10:54 $
#


rpoisline <- function(lambda, win=owin()) {
  win <- as.owin(win)
  # determine circumcircle
  xr <- win$xrange
  yr <- win$yrange
  xmid <- mean(xr)
  ymid <- mean(yr)
  width <- diff(xr)
  height <- diff(yr)
  rmax <- sqrt(width^2 + height^2)/2
  boundbox <- owin(xmid + c(-1,1) * rmax, ymid + c(-1,1) * rmax)
  # generate poisson lines through circumcircle
  n <- rpois(1, lambda * 2 * pi * rmax)
  if(n == 0) {
    X <- psp(numeric(0), numeric(0), numeric(0), numeric(0),
             marks=integer(0), 
             window=win)
    attr(X, "lines") <- infline(p=numeric(0), theta=numeric(0))
    attr(X, "linemap") <- integer(0)
    return(X)
  }
  theta <- runif(n, max= 2 * pi)
  p <- runif(n, max=rmax)
  # compute intersection points with circle
  q <- sqrt(rmax^2 - p^2)
  co <- cos(theta)
  si <- sin(theta)
  X <- psp(x0= xmid + p * co + q * si,
           y0= ymid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= ymid + p * si + q * co,
           marks = seq_len(n),
           window=boundbox, check=FALSE)
  # infinite lines
  L <- infline(p = p + xmid * co + ymid * si,
               theta = theta)
  # clip to window
  X <- X[win]
  # append info
  linemap <- as.integer(marks(X))
  X <- unmark(X)
  attr(X, "lines") <- L
  attr(X, "linemap") <- linemap
  return(X)
}

rjitter.psp <- function(X, radius, ..., clip=TRUE, nsim=1, drop=TRUE) {
  if(nsegments(X) == 0) {
    result <- rep(list(X), nsim)
    result <- simulationresult(result, nsim, drop)
    return(result)
  }
  Xfrom <- endpoints.psp(X, "first")
  Xto   <- endpoints.psp(X, "second")
  if(clip) 
    Window(Xfrom) <- Window(Xto) <- grow.rectangle(Frame(X), radius)
  result <- vector(mode="list", length=nsim)
  for(isim in seq_len(nsim)) {
    Xfrom <- rjitter(Xfrom, radius)
    Xto   <- rjitter(Xto, radius)
    Y <- as.psp(from=Xfrom, to=Xto)
    if(clip)
      Y <- Y[Window(X), clip=TRUE]
    result[[isim]] <- Y
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

