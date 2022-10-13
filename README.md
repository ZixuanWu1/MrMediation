# MrMediation

This is the R package for MrMediation, a Bayesian framework for Mendelian Randomization in the mediation setting.

## Installation

To run the package, you might need to first install RCpp from CRAN

```
install.packages("RCpp")
```

Then the package can be installed in R from this Github repository:

```
library(devtools)
install_github("ZixuanWu1/MrMediation")
```

## Basic Usage

Recall in the Mediation setting, we have $$\Gamma = B \cdot \Gamma + \alpha $$
where $\Gamma$ is the matrix of true SNP effects, $B$ is an upper-triangular matrix with zero diagonalsand $\alpha$ is the horizontal pleiotropy. One might equivalently write this as $$\Gamma = (I + \tilde{B}) \alpha$$

In order to run MrMediation, we will need at least two inputs, Gamma_hat and Sd_hat. The Gamma_hat matrix is an observation of $\Gamma$. Typically it is  a $K \times P$ matrix where $K$ is the number of phenotypes and $P$ is the number of the measured genetic variants. The Sd_hat matrix is the matrix of the standard deviation of noise and is of same dimension as Gamma_hat. 

To estimate the matrix $B$, one can use the following command

```
result = BayesMediation(Gamma_hat, Sd_hat, inv = TRUE)
```

Here ```inv = TRUE``` implies we are estimating B, as opposed to estimating $\tilde{B}$ when ```inv = FALSE ```.

One can also support the correlation among the measurement errors by setting ```cor``` to be the correlation matrix

```
result_cor = BayesMediation(Gamma_hat, Sd_hat, cor = cor_mat, inv = TRUE)
```

In addition, in order to compute the indirect effect of exposures on the outcome, one can simply use the following command

```
result = BayesMediation(Gamma_hat, Sd_hat, cor = cor_mat, indirect = T)
```

Sometimes we might see a warning that the algorithm might not have convergenced. This problem could be solved by using more iterations or changing the initialization method. One can also look at the traceplot of parameters for MCMC diagnosis. For instance, the following function gives the traceplot of the parameter of insterest of all the chains.

```
traceplot(result$raw, par = "B", ind = c(1,2))
```

(Here we are using $B[1,2]$ as an example.)
