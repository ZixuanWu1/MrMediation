#' @keywords internal
#' 
random_upper_tri <- function(K) {
  if (K==1) return(matrix(0))
  Beta <- matrix(0,K,K)
  Beta_upper_tri = runif(K*(K-1)/2,-1,1)
  n<-1;
  for (j in 2:K) {
    for (i in 1:(j-1)) {
      Beta[i,j] <- Beta_upper_tri[n];
      n <- n + 1;
    }
  }
  return(Beta);
}


#' @keywords internal
#' 
result_process <- function(results, K){
  m = length(results)
  len = tail(dim(results[[1]][["B"]]), n = 1)
  
  for(i in 1:m){
    for(j in 1:(len)){
      B_hat = matrix(0, nrow = K, ncol =K)
      for(u in 1:(K - 1)){
        for(r in (u + 1):K){
          B_hat[u, r] = results[[i]][["B"]][u, r, j]
        }
      }
      B_inv = convert(B_hat)
      for(u in 1:(K - 1)){
        for(r in (u + 1):K){
          results[[i]][["B"]][u, r, j]  = B_inv[u, r]
        }
      }
    }
  }
  return(results)
  
}

#' @keywords internal
#'

estimates <- function(results, par, ind, chains = -1, T = 1000){
  if(chains[1] == -1){
    m = length(results)
  }
  else{
    results_sub = list()
    m = length(chains)
    for(i in 1:m){
      results_sub[[i]] = results[[chains[i]]]
    }
    results = results_sub
    
  }
  len = tail(dim(results[[1]][[par]]), n = 1)
  means = rep(0, m)
  vars = rep(0, m)
  ESS = 0
  
  n = len - T
  x = 0
  if(par == "B"){
    
    for(i in 1:m){
      x = c(x, results[[i]][[par]][ind[1], ind[2], (T + 1):len])
      means[i] = mean(results[[i]][[par]][ind[1], ind[2], (T + 1):len])
      vars[i] =  var(results[[i]][[par]][ind[1], ind[2], (T + 1):len])
      ESS = ESS + effectiveSize(results[[i]][[par]][ind[1], ind[2], (T + 1):len])
    }
  }
  else if(par == "sigma"){
    for(i in 1:m){
      x = c(x, results[[i]][[par]][ (T + 1):len])
      means[i] = mean(results[[i]][[par]][(T + 1):len])
      vars[i] =  var(results[[i]][[par]][ (T + 1):len])
      ESS = ESS + effectiveSize(results[[i]][[par]][(T + 1):len])
    }
  }
  else{
    for(i in 1:m){
      x = c(x, results[[i]][[par]][ind, (T + 1):len])
      means[i] = mean(results[[i]][[par]][ind, (T + 1):len])
      vars[i] =  var(results[[i]][[par]][ind, (T + 1):len])
      ESS = ESS + effectiveSize(results[[i]][[par]][ind, (T + 1):len])
    }
  }
  x = x[-1]
  total_mean = sum(means)/m
  B_n = var(means)
  W = mean(vars)
  var = W * (n - 1)/n + B_n
  if(m == 1){
    var = W * (n - 1)/n
  }
  sd= sqrt(var)
  R_hat = sqrt(var/W)
  
  quant = quantile(x, c(0.025, 0.5, 0.975))
  output =matrix(c(total_mean,  var,sd, quant,  ESS, R_hat), nrow = 1)
  colnames(output) = c("mean", "var", "sd","2.5%" ,"50%" ,"97.5%" , "ESS", "Rhat")
  return(output)
}



#' @keywords internal
#'
convert <- function(B){
  K = dim(B)[1]
  return( -(solve(diag(K) +B) - diag(K) ))
}
