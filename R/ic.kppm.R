ic.kppm <- function(object){
  loglike<-object$loglike
  betahat<-object$betahat
  n <- object$n;  W<- object$W; N<-object$N
  p <-sum(betahat!=0)
  if (is.poisson(object$obj.kppm))
    df<-p
  else {
    verifyclass(object$obj.kppm, "kppm")
    cov<-vcov.kppm(object$obj.kppm,what="internals")
    df<-p+sum(diag(as.matrix(cov$J.inv%*%cov$E))) #compute p_approx
  }

  bic = -2*loglike+df*log(n)
  aic = -2*loglike+df*2
  bic.W = -2*loglike+df*log(W)
  bic.N = -2*loglike+df*log(N)

  #NB: aic is cic (and bic is cbic) if not Poisson !

  return(list(loglike=loglike,bic=bic,aic=aic,bic.W=bic.W,bic.N=bic.N,df=df))
}
