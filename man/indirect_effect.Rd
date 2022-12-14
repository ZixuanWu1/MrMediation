% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_effects.R
\name{indirect_effect}
\alias{indirect_effect}
\title{Indirrect effedt}
\usage{
indirect_effect(results, K, warmup = 3000)
}
\arguments{
\item{results}{results from gibbs_wrapper}

\item{K}{number of phenotypes}

\item{warump}{The length of warmup period. Default is 3000}
}
\value{
A matrix of indirect effect and its quantiles.
}
\description{
Estimate the indirrect effect of an exposure on the outcome
}
