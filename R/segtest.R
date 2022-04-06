#'
#'      segtest.R
#'
#'   Monte Carlo test of segregation for multitype patterns
#'
#'    $Revision: 1.6 $ $Date: 2022/04/06 07:35:46 $
#'

segregation.test <- function(X, ...) {
  UseMethod("segregation.test")
}

segregation.test.ppp <- function(X, ..., nsim=19, permute=TRUE,
                                 verbose=TRUE, Xname) {
  if(missing(Xname))
    Xname <- short.deparse(substitute(X))
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  check.1.integer(nsim)
  stopifnot(nsim > 1)
  verboten <- c("at", "relative", "se", "leaveoneout",
                "casecontrol", "case", "control")
  if(any(nyet <- (verboten %in% names(list(...)))))
    stop(paste(ngettext(sum(nyet), "Argument", "Arguments"),
               commasep(sQuote(verboten[nyet])),
               "cannot be used"))
  lam <- intensity(X)
  pbar <- lam/sum(lam)
  np <- npoints(X)
  nt <- length(pbar)
  pbar <- matrix(pbar, byrow=TRUE, nrow=np, ncol=nt)
  if(verbose) cat("Computing observed value... ")
  phat <- relrisk(X, at="points", ..., casecontrol=FALSE)
  obs <- sum((phat-pbar)^2)
  if(verbose) {
    cat(paste("Done.\nComputing", nsim, "simulated values... "))
    pstate <- list()
  }
  sim <- numeric(nsim)
  for(i in 1:nsim) {
    Xsim <- rlabel(X, permute=permute)
    phatsim <- relrisk(Xsim, at="points", ..., casecontrol=FALSE)
    if(permute) pbarsim <- pbar else {
      lamsim <- intensity(Xsim)
      pbarsim <- lamsim/sum(lamsim)
      pbarsim <- matrix(pbarsim, byrow=TRUE, nrow=np, ncol=nt)
    }
    sim[i] <- sum((phatsim - pbarsim)^2)
    if(verbose) pstate <- progressreport(i, nsim, state=pstate)
  }
  if(verbose) cat("Done.\n")
  p.value <- (1+sum(sim >= obs))/(1+nsim)
  names(obs) <- "T"
  out <- list(statistic=obs,
              p.value=p.value,
              method="Monte Carlo test of spatial segregation of types",
              data.name=Xname)
  class(out) <- "htest"
  return(out)
}





