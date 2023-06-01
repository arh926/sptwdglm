
# sptwdglm: An R-package for performing Bayesian Variable selection using Spike and Slab priors for Double Generalized Tweedie Spatial Process Models

<!-- badges: start -->
![Maintainer](https://img.shields.io/badge/maintainer-arh926-blue)
<!-- badges: end -->

`sptwdglm` contains MCMC algorithms for fitting the following models:

Function   | Models
:---- | :-------------
`dglm.autograd.R`   | Double Generalzied Linear Model (DGLM) 
`ssdglm.autograd.R`  | Spike and Slab Priors for DGLMs 
`spdglm.autograd.R`    | Spatial DGLMs
`spssdglm.autograd.R`    | Spike and Slab Priors for Spatial DGLMs

 Variable selection is performed using the function `FDR.R` on the model coefficients. 
 
It supplements the paper titled, "Bayesian Variable Selection in Double Generalized Linear Tweedie Spatial Process Models", New England Journal of Statistics in Data Science: Special Issue on Modern Bayesian Methods with Applications in Data Science (Just accepted). 

All of the above MCMC samplers use a Metropolis Adjusted Langevin Algorithm (MALA, see Girolami and Calderhead, 2011 https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00765.x).

<p align="center">
  <img width="650" height="500" src="https://user-images.githubusercontent.com/73150479/234196908-ea672b3c-9ceb-4472-9d43-1e865d59738d.jpg">
<p>
Figure showing spatial patterns on the left column and logarithm of aggregated Tweedie response on the right column for 10,000 realizations across 100 locations.

## Installation

You can install the development version of `sptwdglm` like so:


``` r
# if you dont have devtools installed
# install.packages("devtools")
devtools::install_github("arh926/sptwdglm")
```

## Example

There are examples contained within every function. Please install the package to view them. 

``` r
require(sptwdglm)

# non-spatial
mc <- Function(response, mean covariates, dispersion covariates, mcmc parameters)
# spatial
mc <- Function(coordinates, response, mean covariates, dispersion covariates, mcmc parameters)

# Diagnostics
plot_mcmc(posterior samples)

# Variable selection through FDR for coefficients
FDR(mean coefficients)
FDR(dispersion coefficients)
```

