#' BayesMediation
#' 
#' Bayes approach for MR Mediation model.
#' 
#' @param Gamma_hat The estimated exposure and outcome effects gamma_hat
#' @param Sd_hat The standarded errors of gamma_hat
#' @param init Starting value for gibbs sampler. Either "Random" or "EM"
#' @param iter The number of iterations for gibbs sampler
#' @param warmup The length of warm-up periods
#' @param second Whether to run the second stage. Default is False.
#' @param inv When inv = False, we are estiming B in Gamma = (I + B) alpha;
#' inv = True, we are estimating B in Gamma = B Gamma + alpha.
#' @param cor The correlation matrix of noise.
#' @param Raw Whether to include the unprocessed raw outputs. Default
#' is F
#' @param total Whether to include the total effects. Default is F
#' @param indirect Whether to include the indirect effects. Default is F.
#'  
#' @return A list with elements
#' \item{summary_first} The summary table of first stage
#' \item{summary_second} The summary table of second stage (included if second  = T)
#' \item{total_effect_first} The total effect of each exposure on the outcome,
#' computed by first stage (included if total = T)
#' \item{total_effect_second} The total effect of each exposure on the outcome,
#' computed by second stage (included if total = T and second = T)
#' \item{indirect_effect_first} The indirect effect of each exposure on the outcome,
#' computed by first stage (included if indirect = T)
#' \item{indirect_effect_second} The indirect effect of each exposure on the outcome,
#' computed by second stage (included if indirect = T and  second = T)
#' \item{raw_first} The unprocessed output from the first stage (included if raw= T)
#' \item{raw_second} The unprocessed output from the second stage (included if raw = T and second = T)
#' 
#' @export
BayesMediation = function(Gamma_hat, Sd_hat, init = "Random", iter = 6000,
                     warmup = 3000, second = F,  inv = FALSE,
                     cor = NULL, Raw = T, total = F, indirect = F) {
  registerDoParallel(cores=4)
  P = length(Gamma_hat[1,])
  K = length(Gamma_hat[, 1])
  
  #Estimate pleiotropy for exposures and outcome
  exps = list()
  for(i in 1: (K-1)) {
    exps[[i]] = zero.centered.em(Gamma_hat[(i + 1), ], 
                                 Sd_hat[(i + 1),], show_outputs = FALSE)
  }
  out = zero.centered.em(Gamma_hat[(1), ], 
                         Sd_hat[(1),], show_outputs = FALSE)
  
  sigma1 = rep(0, K)
  sigma0 = rep(0, K)
  pi = rep(0, K)
  for(i in 1:(K - 1)){
    sigma1[i + 1] = exps[[i]]$S1
    sigma0[i + 1] = exps[[i]]$S2
    pi[i + 1] = exps[[i]]$Pi
  }
  
  sigma1[1] = out$S1 
  sigma0[1] = out$S2
  pi[1] = out$Pi 
  
  #First Stage
  result1 = gibbs_wrapper(iter, warmup, chains = 4,
                            Gamma_hat, Sd_hat, sigma = 0.1, 
                            sigma_1 = sigma1,
                            sigma_0 = sigma0,
                            p = pi,
                            parallel = TRUE, init = init, cor_mat = cor)
  
  
  
  summarys1  = summary_gibbs(result1, c("B", 
                                        "sigma", 
                                        'sigma1', "sigma0", 
                                        "p"), K, T = warmup, inv = inv)
  
  if(max(summarys1[, 8]) > 1.1){
    warning("The maximum R_hat is greater than 1.1. The gibbs sampler might fail
            to converge in the first stage")
  }
  # second stage  (if necessary)
  if(second == TRUE){
    B_hat_vec = summarys1[1:((K^2-K)/2), 1]
    B_hat = matrix(0, K, K)
    B_hat[lower.tri(B_hat, diag = FALSE)] = B_hat_vec
    B_hat = t(B_hat)
    if(inv == FALSE){
      B_hat = convert(B_hat)
    }
    
    #Note (I - B_hat) Gamma_hat = alpha + (I - B_hat) epsilon
    alpha_hat = (diag(K) - B_hat) %*% Gamma_hat
    cor_hat = matrix(nrow = K, ncol = P)
    for(c in 1:P){
      cor_hat[, c] = diag( (diag(K) - B_hat) %*% 
                             diag(Sd_hat[,c])^2 %*% t(diag(K) - B_hat) )
    }
    
    exps = list()
    for(i in 1: (K-1)) {
      exps[[i]] = zero.centered.em(alpha_hat[(i + 1), ], 
                                   sqrt( cor_hat[(i + 1),] ), show_outputs = FALSE)
    }
    out = zero.centered.em(alpha_hat[(1), ], 
                           sqrt( cor_hat[(1),] ), show_outputs = FALSE)
    
    sigma1 = rep(0, K)
    sigma0 = rep(0, K)
    pi = rep(0, K)
    for(i in 1:(K - 1)){
      sigma1[i + 1] = exps[[i]]$S1
      sigma0[i + 1] = exps[[i]]$S2
      pi[i + 1] = exps[[i]]$Pi
    }
    
    sigma1[1] = out$S1 
    sigma0[1] = out$S2
    pi[1] = out$Pi 
    
    result2 = gibbs_wrapper(iter, warmup, chains = 4,
                              Gamma_hat, Sd_hat, sigma = 0.1, 
                              sigma_1 = sigma1,
                              sigma_0 = sigma0,
                              p = pi,
                              parallel = TRUE, init = init, cor_mat = cor)
  
    
    summarys2  = summary_gibbs(result2, c("B", 
                                          "sigma", 
                                          'sigma1', "sigma0", 
                                          "p"), K, T = warmup, inv = inv)
    
    if(max(summarys2[, 8]) > 1.1){
      warning("The maximum R_hat is greater than 1.1. The gibbs sampler might fail
            to converge in the first stage")
    }
    
    ls = list(summary_first = summarys1, summary_second = summarys2)
    if(total){
      ls$total_effect_first = total_effect(result1, K, warmup)
      ls$total_effect_second = total_effect(result2, K, warmup)
    }
    if(indirect){
      ls$indirect_effect_first = indirect_effect(result1, K, warmup)
      ls$indirect_effect_second = indirect_effect(result2, K, warmup)
      
    }
    
    if(Raw == T){
      ls$raw_first = result1
      ls$raw_second = result2
      
    }
    return(ls)
    
  }
  ls = list(summary= summarys1)
  if(total){
    ls$total_effect = total_effect(result1, K, warmup)
  }
  if(indirect){
    ls$indirect_effect= indirect_effect(result1, K, warmup)
  }
  
  if(Raw == T){
    ls$raw = result1
  }
  return(ls)

  
}
