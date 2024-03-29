#
# bw.optim.R
#
#  Class of optimised bandwidths
#  Plotting the object displays the optimisation criterion
#
#  $Revision: 1.33 $  $Date: 2022/01/04 05:30:06 $
#

bw.optim <- function(cv, h,
                     iopt=if(optimum == "min") which.min(cv) else which.max(cv),
                     ...,
                     cvname, hname,
                     criterion="cross-validation",
                     optimum = c("min", "max"), 
                     warnextreme=TRUE, hargnames=NULL,
                     unitname=NULL) {
  if(missing(cvname) || is.null(cvname)) cvname <- deparse(substitute(cv))
  if(missing(hname) || is.null(hname)) hname <- deparse(substitute(h))
  stopifnot(is.numeric(cv))
  stopifnot(is.numeric(h))
  stopifnot(length(h) == length(cv))
  optimum <- match.arg(optimum)
  result <- h[iopt]
  if(warnextreme) {
    optimised <- switch(optimum, min="minimised", max="maximised")
    if(is.infinite(result)) {
      warning(paste(criterion, "criterion was", optimised, "at",
                    hname, "=", as.numeric(result)),
              call.=FALSE)
    } else if((iopt == length(h) || iopt == 1)) {
      warning(paste(criterion, "criterion was", optimised, "at",
                    if(iopt == 1) "left-hand" else "right-hand",
                    "end of interval",
                    paste0(prange(signif(range(h[is.finite(h)]), 3)), ";"), 
                    "use", ngettext(length(hargnames), "argument", "arguments"),
                    paste(sQuote(hargnames), collapse=", "),
                    "to specify a wider interval for bandwidth", sQuote(hname)),
              call.=FALSE)
    }
  }
  attr(result, "cv") <- cv
  attr(result, "h") <- h
  attr(result, "iopt") <- iopt
  attr(result, "labels") <- list(hname=hname, cvname=cvname)
  attr(result, "info") <- list(...)
  attr(result, "criterion") <- criterion
  attr(result, "units") <- unitname
  class(result) <- "bw.optim"
  return(result)
}

print.bw.optim <- function(x, ...) {
  y <- as.numeric(x)
  names(y) <- attr(x, "labels")$hname
  print(y, ...)
  return(invisible(NULL))
}

as.data.frame.bw.optim <- function(x, ...) {
  h <- attr(x, "h")
  cv <- attr(x, "cv")
  df <- data.frame(h, cv)
  labels <- attr(x, "labels")
  colnames(df) <- labels[c("hname", "cvname")]
  info <- attr(x, "info")
  if(length(info) > 0) {
    lenfs <- lengths(info)
    if(any(ok <- (lenfs == nrow(df)))) {
      df <- cbind(df, as.data.frame(info[ok]))
    }
  }
  return(df)
}

as.fv.bw.optim <- function(x) {
  # convert to fv object
  df <- as.data.frame(x)
  dfnames <- colnames(df)
  hname <- dfnames[1L]
  cvname <- dfnames[2L]
  descrip <- c("smoothing parameter",
               paste(attr(x, "criterion"), "criterion"))
  if(ncol(df) > 2)
    descrip <- c(descrip, paste("Additional variable", sQuote(dfnames[-(1:2)])))
  labl <- c(hname, paste0(dfnames[-1L], paren(hname)))
  yexp <- substitute(CV(h), list(CV=as.name(cvname), h=as.name(hname)))
  xfv <- fv(df,
            argu=hname,
            ylab=yexp,
            valu=cvname,
            labl=labl,
            desc=descrip,
            fname=cvname,
            yexp=yexp)
  fvnames(xfv, ".") <- cvname
  unitname(xfv) <- unitname(x)
  return(xfv)
}

plot.bw.optim <- function(x, ...,
                          showopt=TRUE, optargs=list(lty=3, col="blue")) {
  xname <- short.deparse(substitute(x))
  # convert to fv object
  xfv <- as.fv(x)
  # plot cross-validation criterion
  out <- do.call(plot.fv,
                 resolve.defaults(list(x=quote(xfv)),
                                  list(...),
                                  list(main=xname)))
  # Turn off 'showopt' if the x-variable is not the bandwidth
  if(missing(showopt)) {
    argh <- list(...)
    isfmla <- unlist(lapply(argh, inherits, what="formula"))
    if(any(isfmla)) {
      fmla <- argh[[min(which(isfmla))]]
      xvar <- deparse(rhs.of.formula(fmla, tilde=FALSE))
      if(!(identical(xvar, fvnames(xfv, ".x")) || identical(xvar, ".x")))
        showopt <- FALSE
    }
  }
  # show optimal value?
  if(showopt) {
    hoptim <- as.numeric(x)
    if(spatstat.options('monochrome'))
      optargs <- col.args.to.grey(optargs)
    do.call(abline, append(list(v=hoptim), optargs))
  }
  if(is.null(out)) return(invisible(NULL))
  return(out)
}


