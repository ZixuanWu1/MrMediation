#' Gibbs sampler
#' 
#' Given the observed effect and standard error, run gibbs sampler to draw 
#' samples from the posterior density
#' 
#' @param N Number of total iterations
#' @param warmup Number of warm-up iterations.
#' @param chains Number of Chains.
#' @param Gamma_hat The observed effects
#' @param Sd_hat The standard deviation of Gamma_hat
#' @param sigma The prior mean for sigma
#' @param sigma_1 The prior mean for sigma_1
#' @param sigma_0 The prior mean for sigma_0
#' @param p The prior mean for p
#' @param parallel Indicating whether to parallelly run all the chains. 
#' Default is True.
#' @param init Method of Initialization. Must be one of "EM" or "Random". Default
#' is "Random"
#' @param ratio The ratio of initial sigma1 and sigma0 when the initialization method
#' is "Random". Default is 10.
#' @param cor_mat The correlation matrix for noise. Default is null (no correlation)
#' @param Lambda The block diagonal matrix such that each block is the covariance 
#' matrix of the noise for a SNP. Default is NULL (not provided)
#' @param Lambda_inv The block diagonal matrix such that each block is the inverse of 
#' covariance matrix of the noise for a SNP. Default is NULL (not provided).
#' 
#' 
#' @return A list of which the elements are outputs from each chain. Each of these
#' outputs is a list with attributes
#' \item{B}
#' \item{p}
#' \item{sigma1}
#' \item{sigma0}
#' \item{pi}
#' \item{sigma}
#' 
#' @import foreach
#' @import doParallel
#' @import coda
#' @import stats
#' @export

gibbs_wrapper = function(N, warmup, chains = 4,
                             Gamma_hat, Sd_hat, sigma, sigma_1, sigma_0, p,
                             parallel = TRUE, init = "Random", ratio = 10, cor_mat = NULL, Lambda = NULL,
                             Lambda_inv = NULL) {
  P = length(Gamma_hat[1,])
  K = length(Gamma_hat[, 1])
  
  #Set Hyperparameters
  alpha_B = 2
  beta_B = 0.5
  alpha_0 = rep(3, K)
  alpha_1 = rep(3, K)
  stopifnot(length(sigma_0) == K)
  stopifnot(length(sigma_1) == K)
  beta_0 = sigma_0^2 * 2
  beta_1 = sigma_1^2 * 2
  a = rep(2, K)
  stopifnot(length(p) == K)
  b = a/p - a
  
  A = matrix(nrow=K, ncol=P);
  Z = matrix(nrow=K, ncol = P);
  
  
  outputs = list()
  
  #If correlation is non-zero, compute the covariance and inverse covariance
  #matrices, if not provided.
  if(! is.null(cor_mat)){
    if(is.null(Lambda)){
      requireNamespace('magic')
      Lambda = diag( Sd_hat[,1] )%*% cor_mat %*% diag( Sd_hat[,1] )
      for(i in 2:length(Sd_hat[1, ])){
        Lambda = adiag(Lambda, diag( Sd_hat[,i] ) %*% cor_mat %*% diag( Sd_hat[,i]))
      }
      
    }
    if(is.null(Lambda_inv)){
      requireNamespace('magic')
      Lambda_inv = solve( diag(Sd_hat[,1] )%*% cor_mat %*% diag(Sd_hat[,1]) )
      for(i in 2:length(Sd_hat[1,])){
        Lambda_inv = adiag(Lambda_inv, solve( diag( Sd_hat[,i] )%*% cor_mat %*% 
                                                diag( Sd_hat[,i]) ))
      }
      
    }
  }

  
  each_chain_task = function(i) {
    
    #Set initial values for parameters.
    
    sigma = runif(1)
    if(init == "Random"){
      sigma1 = runif(K)
      sigma0 = sigma1 / ratio
      p = runif(K, min = 0, max = 0.1)
    }
    else{
      sigma1 = sigma_1
      sigma0 = sigma_0 
      p = p
    }
    B = random_upper_tri(K)
    
    if(is.null(cor_mat)){
      this_output = gibbs_sampler(Gamma_hat, Sd_hat, N, B, sigma, sigma1,
                                  sigma0, p, A, Z,
                                  alpha_B, beta_B, alpha_0, alpha_1, beta_0, beta_1, a, b)
    }
    else{
      this_output = gibbs_sampler_with_corr(Gamma_hat,
                                            Sd_hat, cor_mat, N, B,
                                            sigma, sigma1,
                                            sigma0, p, A, Z,
                                            alpha_B, beta_B, 
                                            alpha_0, alpha_1, beta_0,
                                            beta_1, a, b, Lambda, Lambda_inv)
    }

    this_output
    
  }
  
  
  print(Sys.time())
  if (parallel){
    outputs = foreach(i=1:chains) %dopar% each_chain_task(i)
  } else {
    outputs = lapply(1:1, each_chain_task)
  }
  print(Sys.time())
  
  return (outputs)
}

