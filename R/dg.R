#
#     dg.S
#
#    $Revision: 1.25 $	$Date: 2022/05/21 08:53:38 $
#
#     Diggle-Gratton pair potential
#
#
DiggleGratton <- local({

  # .... auxiliary functions ......

  diggraterms <- function(X, Y, idX, idY, delta, rho, splitInf=FALSE) {
    stopifnot(is.numeric(delta))
    stopifnot(is.numeric(rho))
    stopifnot(delta < rho)
    ## sort in increasing order of x coordinate
    oX <- fave.order(X$x)
    oY <- fave.order(Y$x)
    Xsort <- X[oX]
    Ysort <- Y[oY]
    idXsort <- idX[oX]
    idYsort <- idY[oY]
    nX <- npoints(X)
    nY <- npoints(Y)
    ## call C routine
    if(!splitInf) {
      ## usual case: allow cif to be zero because of hard core
      out <- .C(SC_Ediggra,
                nnsource = as.integer(nX),
                xsource  = as.double(Xsort$x),
                ysource  = as.double(Xsort$y),
                idsource = as.integer(idXsort),
                nntarget = as.integer(nY),
                xtarget  = as.double(Ysort$x),
                ytarget  = as.double(Ysort$y),
                idtarget = as.integer(idYsort),
                ddelta   = as.double(delta),
                rrho     = as.double(rho),
                values   = as.double(double(nX)),
                PACKAGE="spatstat.core")
      answer <- integer(nX)
      answer[oX] <- out$values
    } else {
      ## split off the hard core terms and return them separately
      out <- .C(SC_ESdiggra,
                nnsource = as.integer(nX),
                xsource  = as.double(Xsort$x),
                ysource  = as.double(Xsort$y),
                idsource = as.integer(idXsort),
                nntarget = as.integer(nY),
                xtarget  = as.double(Ysort$x),
                ytarget  = as.double(Ysort$y),
                idtarget = as.integer(idYsort),
                ddelta   = as.double(delta),
                rrho     = as.double(rho),
                positive = as.double(double(nX)),
                hardcore = as.integer(integer(nX)),
                PACKAGE="spatstat.core")
      answer <- integer(nX)
      hardcore <- logical(nX)
      answer[oX] <- out$positive
      hardcore[oX] <- as.logical(out$hardcore)
      attr(answer, "hardcore") <- hardcore
    }
    return(answer)
  }

  # .......... template object ..........
  
  BlankDG <- 
  list(
         name     = "Diggle-Gratton process",
         creator  = "DiggleGratton",
         family    = "pairwise.family",  #evaluated later
         pot      = function(d, par) {
                       delta <- par$delta
                       rho <- par$rho
                       above <- (d > rho)
                       inrange <- (!above) & (d > delta)
                       h <- above + inrange * (d - delta)/(rho - delta)
                       return(log(h))
                    },
         par      = list(delta=NULL, rho=NULL),  # to be filled in later
         parnames = list("lower limit delta", "upper limit rho"),
         hasInf   = TRUE, 
         selfstart = function(X, self) {
           # self starter for DiggleGratton
           nX <- npoints(X)
           if(nX < 2) {
             # not enough points to make any decisions
             return(self)
           }
           md <- minnndist(X)
           if(!is.na(delta <- self$par$delta)) {
             # value fixed by user or previous invocation
             # check it
             if(md < delta)
               warning(paste("Hard core distance delta is too large;",
                             "some data points will have zero probability"))
             return(self)
           }
           if(md == 0) 
             warning(paste("Pattern contains duplicated points:",
                           "hard core distance delta must be zero"))
           # take hard core = minimum interpoint distance * n/(n+1)
           deltaX <- md * nX/(nX+1)
           DiggleGratton(delta=deltaX, rho=self$par$rho)
         },
         init = function(self) {
           delta <- self$par$delta
           rho   <- self$par$rho
           if(!is.numeric(rho) || length(rho) != 1L)
             stop("upper limit rho must be a single number")
           stopifnot(is.finite(rho))
           if(!is.na(delta)) {
             if(!is.numeric(delta) || length(delta) != 1L)
               stop("lower limit delta must be a single number")
             stopifnot(delta >= 0)
             stopifnot(rho > delta)
           } else stopifnot(rho >= 0)
         },
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           kappa <- as.numeric(coeffs[1L])
           return(list(param=list(kappa=kappa),
                       inames="exponent kappa",
                       printable=dround(kappa)))
         },
         valid = function(coeffs, self) {
           kappa <- as.numeric(coeffs[1L])
           return(is.finite(kappa) && (kappa >= 0))
         },
         project = function(coeffs, self) {
           kappa <- as.numeric(coeffs[1L])
           if(is.finite(kappa) && (kappa >= 0))
             return(NULL)
           return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           rho <- self$par$rho
           if(all(is.na(coeffs)))
             return(rho)
           kappa <- coeffs[1L]
           delta <- self$par$delta
           if(abs(kappa) <= epsilon)
             return(delta)
           else return(rho)
         },
         hardcore = function(self, coeffs=NA, epsilon=0, ...) {
           return(self$par$delta)
         },
       version=NULL, # evaluated later
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                      splitInf=FALSE, ...) {
         ## fast evaluator for DiggleGratton interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for DiggleGratton")
         delta <- potpars$delta
         rho   <- potpars$rho
         idX <- seq_len(npoints(X))
         idU <- rep.int(-1L, npoints(U))
         idU[EqualPairs[,2L]] <- EqualPairs[,1L]
         values <- diggraterms(U, X, idU, idX, delta, rho, splitInf)
         result <- log(pmax.int(0, values))
         result <- matrix(result, ncol=1L)
         if(!splitInf)
           attr(result, "-Inf") <- attr(values, "hardcore")
         return(result)
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         rho   <- self$par$rho
         delta <- self$par$delta
         width <- rho - delta
         kappa <- coeffs[1L]
         ans <- pi * (rho^2
                      - 2 * rho* width/(kappa + 1)
                      + 2 * width^2/((kappa + 1) * (kappa + 2)))
         return(ans)
       }
  )
  class(BlankDG) <- "interact"

  DiggleGratton <- function(delta=NA, rho) {
    instantiate.interact(BlankDG, list(delta=delta, rho=rho))
  }

  DiggleGratton <- intermaker(DiggleGratton, BlankDG)

  DiggleGratton
})
