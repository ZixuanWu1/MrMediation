#' Perform EM Algorithm on Pleiotropy
#'
#' Fit a two-component zero-mean Gaussian mixture on a sequence of observed 
#' values
#'
#' @param values A vector of observations.
#' @param std.error The standard error of noise of observations.
#' @param initial_slab_prob The initial value for the proportion of slab (the
#' Gaussian component with higher variance). Default is .05.
#' @param initial_sigma1_squared The initial value for the variance of the
#' slab component. Default is .1.
#' @param initial_sigma2_squared The initial value for the variance of the 
#' spike component. Default is .0001
#' @param iters The number of iterations of EM algorithm. Default is 100
#' @param optimizer The method used in Maximization Step. Must be one of "Brent",
#' "Nelder-Mead" and "L-BFGS-B". Default is "Brent". 
#' @param show_outputs Whether to included a plot of the fitted density and the
#' smoothed data density. Default is FALSE.
#'
#' @return A list of elements 
#' \item{S1} The standard deviation of the slab component
#' \item{S2} The standard deviation of the spike component
#' \item{Pi} The proportion of the slab component
#' 
#' @export
zero.centered.em <- function(values, std.errors,
                             initial_slab_prob = .05,
                             initial_sigma1_squared = .1,
                             initial_sigma2_squared = .0001,
                             iters = 100,
                             optimizer = "Brent",
                             show_outputs = FALSE) {
  
  #check that dimensions agree
  stopifnot(length(values) == length(std.errors))
  
  slab_prob <- initial_slab_prob;
  sigma1.sq <- initial_sigma1_squared;
  sigma2.sq <- initial_sigma2_squared;
  n <- length(values)
  
  estimates <- data.frame()
  
  
  se.to.weights <- function(se) {
    w <- 1/se
    w <- w/sum(w)
    return(w)
  }
  
  #this is main loop:
  for(iter in 1:iters) {
    # unnormalized membership probabilities:
    w1_un <- slab_prob * dnorm(values, 0, sqrt(sigma1.sq + std.errors^2))
    w2_un <- (1-slab_prob) * dnorm(values, 0, sqrt(sigma2.sq + std.errors^2))
    #normalize them:
    w1 <- w1_un / (w1_un + w2_un)
    w2 <- 1 - w1
    
    log_likelihood <- sum(log(w1_un + w2_un))
    
    
    #maximize slab_prob:
    slab_prob <- mean(w1)
    #use built in optimize function from stats library
    if(optimizer =="Nealder-Mead"){
      sigma1.sq <-inv.logit(optim(0.1, sigma_objective_n, 
                                  method = optimizer, weights = w1, 
                                  values = values, std.errors =  std.errors,
                                  control = list(warn.1d.NelderMead = FALSE))$par) 
      sigma2.sq <- inv.logit(optim(0.01,sigma_objective_n,
                                   method = optimizer, weights = w2, 
                                   values = values, std.errors =  std.errors,
                                   control = list(warn.1d.NelderMead = FALSE))$par)
      
    }
    
    else if(optimizer == "L-BFGS-B"){
      sigma1.sq <- optim(0.1, sigma_objective, sigma_gradient,
                         method = optimizer, weights = w1, 
                         values = values, std.errors =  std.errors, 
                         lower = 0, upper =sigma1.sq * 10 )$par
      sigma2.sq <- optim(0.01,sigma_objective, sigma_gradient,
                         method = optimizer, weights = w2, 
                         values = values, std.errors =  std.errors, 
                         lower = 0, upper =sigma2.sq * 10)$par
      
    }
    else{
      sigma1.sq <- optim(0.1, sigma_objective, 
                         method = optimizer, weights = w1, 
                         values = values, std.errors =  std.errors, 
                         lower = 0, upper =sigma1.sq * 10 )$par
      sigma2.sq <- optim(0.01,sigma_objective,
                         method = optimizer, weights = w2, 
                         values = values, std.errors =  std.errors, 
                         lower = 0, upper =sigma2.sq * 10)$par
      
    }
    
    
    #if( as.character(readline())=="c" ) {return()}
    estimates<-rbind(estimates,
                     cbind(iter,
                           c(slab_prob,
                             sqrt(sigma1.sq),
                             sqrt(sigma2.sq),
                             log_likelihood
                           ),
                           c("Pi","S1","S2","log_likelihood")))
  }
  #the rest is just formatting output/orintouts
  names(estimates) <- c("n","estimate","variable")
  estimates[,1]<-as.numeric(estimates[,1])
  estimates[,2]<-as.numeric(estimates[,2])
  if (show_outputs) {
    print(ggplot(estimates[estimates$variable!="log_likelihood",], aes(n, estimate, color=factor(variable))) +
            geom_line() + ggtitle("EM Results")
    )
    print(ggplot(estimates[estimates$variable=="log_likelihood",], aes(n, estimate, color=factor(variable))) +
            geom_line() + ggtitle("EM Log Likelihood")
    )
    graphics::curve(
      slab_prob*dnorm(x, 0, sqrt(sigma1.sq)) +
        (1-slab_prob)*dnorm(x, 0, sqrt(sigma2.sq)),
      min(values), max(values),
      xlab="beta", ylab="density"
    )
    title("Fitted Distribtuion vs Smoothed Data Density")
    lines(density(values, weights = se.to.weights(std.errors)), col="blue")
  }
  
  result =   list(
    S1 = utils::tail(estimates[estimates$variable=="S1", "estimate"] , 1),
    S2 = utils::tail(estimates[estimates$variable=="S2", "estimate"] , 1),
    Pi = utils::tail(estimates[estimates$variable=="Pi", "estimate"] , 1)
  )
  
  return(result)
}


#' The log likelihood function in the M-step
#' 
#' @keywords internal
#' 
sigma_objective <- function(sigma.sq, weights, values, std.errors) {
  n <- length(weights)
  return(sum(
    sapply(1:n, function(i) weights[i] * (log(sigma.sq + std.errors[i]^2) +
                                            values[i]^2 /(sigma.sq + std.errors[i]^2) )
    )))
}
