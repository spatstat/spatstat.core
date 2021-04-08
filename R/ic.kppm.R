ic.kppm <- function(object,X){
  
  if (is.poisson(object)){
    loglike<-object$maxlogpl
    betahat<-object$coef
    n=X$n
  }
  else {
    loglike<-object$po$maxlogpl
    betahat<-object$po$coef
    n <- object$X$n
  } 
  p <-length(betahat)

  if (is.poisson(object))
    df<-p
  else {
    verifyclass(object, "kppm")
    cov<-vcov.kppm(object,what="internals")
    df<-p+sum(diag(as.matrix(cov$J.inv%*%cov$E))) #compute p_approx
  }

  cbic = -2*loglike+df*log(n)
  cic = -2*loglike+df*2
  #NB: cbic is BIC and cic is AIC in case of Poisson process
  
  return(list(loglike=loglike,cbic=cbic,cic=cic,df=df))
}
