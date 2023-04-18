
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

All samplers use Metropolis Adjusted Langevin Algorithm (MALA, see Girolami and Calderhead, 2011 \url{https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00765.x}). Variable selection is performed using the function `FDR.R` on the model coefficients.

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

mc <- Function(coordinates, repsponse, mean covariates, dispersion covariates, mcmc parameters)

# Variable selection through FDR for coefficients
FDR(mean coefficients)
FDR(dispersion coefficients)
```

