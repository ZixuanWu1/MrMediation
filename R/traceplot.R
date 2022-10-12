#' Produce traceplot for the specified parameter
#' 
#' @param results The output of gibbs_wrapper.
#' @param pars A parameter to be ploted. Must be one of "B", "sigma", "sigma1", 
#' "sigma0" and "pi".
#' @param ind The index of the parameter. For "B" it is a two-vector. For "sigma1",
#' "sigma0" and "pi" it is a single number. For "sigma" ind is ignored.
#' @param chains The indices of chains to be included in the plot. Default is -1,
#' which uses all the chains.
#' @param ylim The range of y-axis in the final plot
#' @param T The number of warmup period.
#' 
#' @return A traceplot of a scalar parameter after warmup-period
#' 
#' @import ggplot2
#' @export
traceplot <- function(results, par, ind, chains = -1, ylim = NA, T = 1000){
  par(mar=c(5, 4, 4, 7), xpd=TRUE)
  len = tail(dim(results[[1]][[par]]), n = 1)
  
  #Consider all chains if not specified.
  if(chains[1] == -1){
    chains = c(1:length(results))
  }
  
  #Give each chain a name
  name = rep(0, length(chains))
  for(i in 1:length(chains)){
    name[i] = paste("chain", toString(chains[i]))
  }
  
  #Traceplot for B[,]
  if(par == "B"){
    plot((T + 1):len, results[[chains[1]]][[par]][ind[1], ind[2], ( (T+1 ):len)], type = "l",
         col = 1, ylab = paste(par, "[", toString(ind[1]), ",", toString(ind[2]), "]", 
                               sep = ""), xlab = "iterations", lwd=2, cex=1.2,)
    
    if(length(chains) > 1){
      for(i in 2:length(chains)){
        lines((T + 1):len, results[[chains[i]]][[par]][ind[1], ind[2], ( (T+1):len)], type = "l",
              col = i)
      }
      
    }
    
  }
  #Traceplot for sigma
  else if(par == "sigma"){
    plot((T + 1):len, results[[chains[1]]][[par]][(T+1):len], type = "l",
         col = 1, ylab = par, xlab = "iterations")
    if(length(chains) > 1){    for(i in 2:length(chains)){
      lines((T + 1):len, results[[chains[i]]][[par]][( (T+1):len)], type = "l",
            col = i)
    }}
    
    
  }
  #Traceplot for sigma1, sigma0 or pi
  else{
    if(is.na(ylim)){
      plot((T + 1):len, results[[chains[1]]][[par]][ind, ( (T+1):len)], type = "l",
           col = 1, ylab = paste(par, "[", toString(ind), "]", 
                                 sep = ""), xlab = "iterations")
    }
    else{
      plot((T + 1):len, results[[chains[1]]][[par]][ind, ( (T+1):len)], type = "l",
           col = 1, ylab = paste(par, "[", toString(ind), "]", 
                                 sep = ""), xlab = "iterations", ylim = ylim)
    }
    
    if(length(chains) > 1){    for(i in 2:length(chains)){
      lines((T + 1):len, results[[chains[i]]][[par]][ind, ( (T+1):len)], type = "l",
            col = i)
    }}
    
  }
  
  legend("bottomright", legend = name,inset = c(-0.3, 0), col =c(1: length(chains)),
         bty = "n", xpd=TRUE,  lty = 1)
}

