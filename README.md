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

Recall in the Mediation setting, we have $$Gamma = B \Gamma + \alpha + \epsilon,$$
where $\Gamma$ is the matrix of true SNP effects, $B$ is an upper-triangular matrix with zero diagonals, $\alpha$ is the horizontal pleiotropy and $\epsilon$ is the noise. One might equivalently write this as $$\Gamma = (I + \tilde{B}) \alpha$$

In order to run MrMediation, we will need at least two inputs, Gamma_hat and Sd_hat. The Gamma_hat matrix is an observation of $\Gamma$. Typically it is  a $K \times P$ matrix where $K$ is the number of phenotypes and $P$ is the number of the measured genetic variants. The Sd_hat matrix is the matrix of the standard deviation of noise and is of same dimension as Gamma_hat. 

To estimate the matrix $B$, one can use the following command

```
BayesMediation(Gamma_hat, Sd_hat, inv = T)
```

Here ```inv``` is 

In order to compute the indirect effect of exposures on the outcome, one can simply use the following command

```
BayesMediation(Gamma_hat, Sd_hat, cor = cor_mat, indirect = T)
```




