#'
#'   intensity.interact.R
#'
#'    Intensity approximations for pairwise interaction models
#'    
#'    

intensity.interact <- function(X,beta=NULL,interaction.coef=NULL, method=c("Poisson","DPP"), ...) {
  method <- match.arg(method)
  if(!identical(X$family$name, "pairwise")) stop("Only pairwise interactions are considered")
  
  if(is.null(beta)||is.null(interaction.coef)){
    stop("The spatial trend beta and the interaction parameters interaction.coef must be supplied")}
  
  if(sum(is.na(X$par))!=0){stop("All irregular parameters of the interaction must be supplied")}
  
  if(method=="Poisson"){lambda <- PoisSaddlePairwise(beta, X,interaction.coef)}
  
  if(method=="DPP"){
    if(length(beta)!=1){stop("The DPP approximation is only implemented for a stationary pairwise interaction")}
    lambda <- DPPSaddlePairwise(beta,X,interaction.coef)
  }
  
  return(lambda)
  
}
