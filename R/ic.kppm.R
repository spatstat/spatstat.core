ic.kppm <- function(object){
  loglike<-object$loglike
  betahat<-object$betahat
  n <- object$n
  p <-sum(betahat!=0)
  if (is.poisson(object$obj.kppm))
    df<-p
  else {
    verifyclass(object$obj.kppm, "kppm")
    cov<-vcov.kppm(object$obj.kppm,what="internals")
    df<-p+sum(diag(as.matrix(cov$J.inv%*%cov$E))) #compute p_approx
  }

  cbic = -2*loglike+df*log(n)
  cic = -2*loglike+df*2
  #NB: cbic is BIC and cic is AIC in case of Poisson process
  
  return(list(loglike=loglike,cbic=cbic,cic=cic))
}
