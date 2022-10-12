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
  len = tail(dim(results[[1]][["B"]]), n = 1)
  effects = matrix(0, (K - 1), m * (len - warmup))
  for(i in 1:m){
    for(j in 1:(len - warmup)){
      temp_B = results[[i]][["B"]][,,(j + warmup)]
      for(u in 1:(K-1) )
        effects[u,( (i - 1) * (len - warmup) + j )] = compute_total(temp_B,u)
    }
  }
  df = matrix(nrow = (K - 1), ncol = 4)
  for(i in 1:(K - 1)){
    df[i, 1] = mean(effects[i,])
    df[i, 2:4] = quantile(effects[i,], c(0.025, 0.5, 0.975))
  }
  colnames(df) = c("mean", "2.5%", "50%", "97.5")
  row.names(df) = paste(rep("exp", (K - 1)), 1:(K - 1), sep = "")
  return(df)
  
  
  
}


#' Indirect effect
#' 
#' Estimate the indirect effect of an exposure on the outcome
#' 
#' @param results results from gibbs_wrapper
#' @param K number of phenotypes
#' @param warump The length of warm-up period. Default is 3000
#'  
#' @return A matrix of indirect effect and its quantiles.
#' 
#' @export

indirect_effect = function(results, K, warmup = 3000) {
  results = result_process(results, K)
  m = length(results)
  len = tail(dim(results[[1]][["B"]]), n = 1)
  effects = matrix(0, (K - 1), m * (len - warmup))
  for(i in 1:m){
    for(j in 1:(len - warmup)){
      temp_B = results[[i]][["B"]][,,(j + warmup)]
      for(u in 1:(K-1) )
        effects[u,( (i - 1) * (len - warmup) + j )] = compute_total(temp_B,u) - temp_B[1, (dim(temp_B)[1] + 1 - u)]
    }
  }
  df = matrix(nrow = (K - 1), ncol = 4)
  for(i in 1:(K - 1)){
    df[i, 1] = mean(effects[i,])
    df[i, 2:4] = quantile(effects[i,], c(0.025, 0.5, 0.975))
  }
  colnames(df) = c("mean", "2.5%", "50%", "97.5")
  row.names(df) = paste(rep("exp", (K - 1)), 1:(K - 1), sep = "")
  return(df)
  
  
}

#' @keywords internal
#'
compute_total <- function(B, k){
  p = (dim(B)[1] + 1 - k)
  B = B[1:p, 1:p]
  K = length(B[,1])
  if(K == 2){
    return(B[1,2])
  }
  else{
    return(B[1, K] + B[(K - 1), K] * compute_total(B[(1:(K - 1)), (1:(K-1))], 1))
  }
}
