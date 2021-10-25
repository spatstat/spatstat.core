#'
#'   bw.CvLHeat
#'
#'  Cronie-van Lieshout bandwidth selection for Diffusion smoothing
#'
#'  Copyright (c) 2020 Adrian Baddeley, Tilman Davies and Suman Rakshit
#'  GNU Public Licence >= 2.0

bw.CvLHeat <- function(X, ..., srange=NULL, ns=16, sigma=NULL,
                     leaveoneout=TRUE, verbose=TRUE) {
  #' compute intensity estimates
  b <- HeatEstimates.ppp(X, ..., srange=srange, ns=ns, sigma=sigma,
                         leaveoneout=leaveoneout, verbose=verbose)
  lambda <- b$lambda
  h      <- b$h
  hname  <- b$hname
  #' compute Cronie-van Lieshout criterion
  AW <- area.owin(Window(X))
  CV <- (rowSums(1/lambda) - AW)^2
  iopt <- which.min(CV)
  result <- bw.optim(CV, h, iopt,
                     criterion="Cronie-van Lieshout criterion",
                     hname=hname)
  return(result)
}
