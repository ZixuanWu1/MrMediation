#' Produce a summary table of the output from gibbs_wrapper.
#' 
#' @param result The output of gibbs_wrapper or gibbs_wrapper_cor
#' @param pars A vector of parameters to be included in the summary table
#' @param K The number of mediation layers
#' @param T The warmup period. Default is 1000.
#' @param inv When inv is TRUE, the function shows summary statistics for B in
#' Gamma = B Gamma + alpha, otherwise it shows summary statistics for B in
#' Gamma = (I + B) alpha. Default is False
#' 
#' @return A table where each row corresponds to a parameter. It has columns:
#' \item{mean} The average of samples after warm-up
#' \item{variance} The variance of samples after warm-up
#' \item{sd} The standard deviation of samples after warm-up
#' \item{2.5\%} The 2.5\% quantile of samples after warm-up
#' \item{50\%} The 50% quantile of samples after warm-up
#' \item{97.5\%} The 97.5\% quantile of samples after warm-up
#' \item{ESS} The Effective sample sizes (combined from all the chains)
#' \item{R_hat} A measure of convergence. A value below 1.1 indicates sign
#' of convergence
#' 
#' @export
summary_gibbs <- function(result, pars, K, T = 1000, inv = FALSE){
  df = rep(0, 8)  
  
  #If inv is TRUE, then convert B into the other form.
  
  if(inv == TRUE & K >= 3){
    result = result_process(result, K)
  }
  for(par in pars){
    
    #Compute statistics for B
    if(par == "B"){
      for(i in 1:(K - 1)){
        for(j in (i + 1): K){
          es = estimates(result, "B", ind = c(i, j), T = T)
          rownames(es) = paste(par, "[", toString(i), ",", toString(j), "]", sep = "")
          df = rbind(df,es)
        }
        
      }
    }
    
    #Compute statistics for sigma
    else if(par == "sigma"){
      es = estimates(result, "sigma", 1, T = T)
      rownames(es) = "sigma"
      df = rbind(df,es)
    }
    #Compute statistics for sigma1, sigma0 and pi.
    else{
      for(i in 1:K){
        es = estimates(result, par, i, T = T)
        rownames(es) = rownames(es) = paste(par,"[", toString(i), "]", sep = "")
        df = rbind(df,es)
      }
    }
  }
  
  return(as.data.frame(df[-1, ]))
}

