
# sptwdglm: An R-package for performing Bayesian Variable selection using Spike and Slab priors for Double Generalized Tweedie Spatial Process Models

<!-- badges: start -->
<!-- badges: end -->

`sptdwglm` contains MCMC algorithms for fitting teh following models 

(a) Double Generalzied Linear Model (DGLM) 
(b) Spike and Slab Priors for DGLMs 
(c) Spatial DGLMs and 
(d) Spike and Slab for Spatial DGLMs

## Installation

You can install the development version of `sptwdglm` like so:

``` r
# if you dont have devtools installed
# install.packages("devtools")
devtools::install_github("arh926/sptwdglm")
```

## Example

There are examples contained within every function listed. I list the particular example pertaining to variable selection in DGLM

``` r
require(sptwdglm)
require(tweedie)
require(mvtnorm)
require(MASS)

#################
# Generate Data #
#################
N = 1e4
L = 1e2

coords = matrix(runif(2*L), nc=2)
par(mfcol=c(1,1))
# plot(coords)
sigma2.true = 1
phis.true = 3
Delta = as.matrix(dist(coords))
Sigma = sigma2.true*exp(-phis.true*Delta)
w.true = mvrnorm(1, rep(0, L), Sigma)

if(N > L) index = sample(1:L, N, replace = TRUE) else if(N == L) index = sample(1:L, N, replace = FALSE)

# Design matrices
z = x = cbind(1, rnorm(N), rnorm(N), rnorm(N), rnorm(N), rnorm(N), rnorm(N))
## scaling design matrices
x[,-1] = apply(x[,-1], 2, function(s) (s - mean(s))/sd(s))
z[,-1] = apply(z[,-1], 2, function(s) (s - mean(s))/sd(s))

p = ncol(x)
q = ncol(z)

# true coefficients
beta0 = 1
beta1 = 1.5
beta2 = 0.01
beta3 = 1.4
beta4 = 1.1
beta5 = 0.01
beta6 = 2.5
beta.true = c(beta0, beta1, beta2, beta3, beta4, beta5, beta6)
mu_sim.sp = exp(x %*% beta.true + w.true[index])

gamma0 = 1
gamma1 = 0.01
gamma2 = 1.5
gamma3 = 1.1
gamma4 = 0.01
gamma5 = -2.5
gamma6 = 0.01
gamma.true = c(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)
phi_sim = exp(z %*% gamma.true)
 
range(mu_sim.sp)
range(phi_sim)

xi.true = 1.5
 
y_sim = rtweedie(N, xi = xi.true, mu = mu_sim.sp, phi = phi_sim)
sum(y_sim == 0)/N # proportion of zeros
y = y_sim
x = x
z = z
# mcmc parameters
niter = 1e5
nburn = niter/2
report = 5e2

# hyperparameters
lower.xi = 1
upper.xi = 2

system.time(mc_ss <- spssdglm.autograd(coords = coords, y = y, x = x, z = z,
                                    niter = niter, nburn= nburn, report = report, thin = 20,
                                    index = index, lower.xi = lower.xi, upper.xi = upper.xi,
                                    verbose = TRUE))
# Check for convergence
plot_mcmc(samples = mc_ss$xi.mcmc, true = 1.5, col = "blue", cnames = "xi")
```

