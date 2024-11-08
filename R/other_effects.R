#' Total effect
#' 
#' Estimate the total effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#' @param warump The length of warmup period. Default is 3000
#'  
#' @return A matrix of total effect and its quantiles.
#' 
#' @export

total_effect = function(results, K, warmup = 3000) {
  results = result_process(results, K)
  m = length(results)
  len = utils::tail(dim(results[[1]][["B"]]), n = 1)
  effects = matrix(0, (K - 1), m * (len - warmup))
  for(i in 1:m){
    for(j in 1:(len - warmup)){
      temp_B = results[[i]][["B"]][,,(j + warmup)]
      effects[,( (i - 1) * (len - warmup) + j )] = compute_total(temp_B)[1, dim(temp_B)[1]:2]
    }}
  df = matrix(nrow = (K - 1), ncol = 4)
  for(i in 1:(K - 1)){
    df[i, 1] = mean(effects[i,])
    df[i, 2:4] = quantile(effects[i,], c(0.025, 0.5, 0.975))
  }
  
  return(df)
  
  
}


#' Indirrect effedt
#' 
#' Estimate the indirrect effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#‘ @param path A vector representing the pathway of interests
#' @param quantiles A vector of quantiles of interests
#' @param warump The length of warmup period. Default is 3000
#' 
#' @return A matrix of indirect effect and its quantiles.
#' 
#' @export

indirect_effect = function(results, K, path = NULL, quantiles =  c(0.025, .5, 0.975),  warmup = 3000) {
  
  results = result_process(results, K)
  m = length(results)
  len = utils::tail(dim(results[[1]][["B"]]), n = 1)
  
  if(is.null(path)){
    effects = matrix(0, K, m * (len - warmup))
    for(i in 1:m){
      for(j in 1:(len - warmup)){
        temp_B = results[[i]][["B"]][,,(j + warmup)]
        effects[,( (i - 1) * (len - warmup) + j )] = (compute_total(temp_B) - temp_B)[1,]
      }
    }
    df = matrix(nrow = K, ncol = 4)
    for(i in 1:K){
      df[i, 1] = mean(effects[i,])
      df[i, 2:4] = quantile(effects[i,], quantiles)
      
    }
    df = as.data.frame(df)
    for(i in 1:K){
      row.names(df)[i] = paste("Exp", toString(i))
    }

    
  } else{
    effects = NULL
    for(i in 1:m){
      for(j in 1:(len - warmup)){
        temp_B = results[[i]][["B"]][,,(j + warmup)]
        cur = 1
        for (ele in 1:(length(path) - 1)){
          cur = cur * temp_B[path[ele + 1], path[ele]]
        }
        effects = c(effects, cur)
      }
    }
    df = c(mean(effects), quantile(effects, quantiles)  )
    df = as.data.frame(matrix(df, nrow= 1))

  }
  colnames(df)[1] = "Mean"
  for(j in 2:ncol(df)){
    colnames(df)[j] = paste(toString(quantiles[j - 1] * 100), "%", sep = "")
  }
  
  return(df)
  
  
}

#' @keywords internal
#'
compute_total <- function(B){
  K = length(B[,1])
  return( solve(diag(K) - B)  - diag(K) )
}
