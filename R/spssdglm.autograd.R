#' Tweedie Compound Poisson-Gamma (CP-g) Double Generalized Linear Spatial Process Models
#'
#' Bayesian Variable Selection for the Double Generalized Linear Model: \eqn{\log(\mu(s))=x(s)^T\beta+f(s)^T w(s)} and \eqn{\log(\phi)=z^T\gamma}. Parameters not listed below are optional.
#'
#' @param coords coordinates for observed process (ex. latitude-longitude)
#' @param y observed response
#' @param x covariates for the mean model
#' @param z covariates for the dispersion model
#' @param lower.xi lower bound for the index parameter, \eqn{\xi}
#' @param upper.xi upper bound for the index parameter, \eqn{\xi}
#' @keywords
#' @import stats tweedie coda Matrix
#' @importFrom mvtnorm dmvnorm
#' @export
#' @examples
#' require(tweedie)
#' require(mvtnorm)
#' require(coda)
#'
#' # Generate Data
#' N = 1e4
#' L = 1e2
#'
#' coords = matrix(runif(2*L), nc=2)
#' par(mfcol=c(1,1))
#' # plot(coords)
#' sigma2.true = 1
#' phis.true = 3
#' Delta = as.matrix(dist(coords))
#' Sigma = sigma2.true*exp(-phis.true*Delta)
#' w.true = crossprod(chol(Sigma), rnorm(L))
#'
#' if(N > L){
#' index = sample(1:L, N, replace = TRUE)
#' }else if(N == L){
#' index = sample(1:L, N, replace = FALSE)
#' }
#'
#' # Design matrices
#' z = x = cbind(1, rnorm(N), rnorm(N), rnorm(N), rnorm(N), rnorm(N), rnorm(N))
#' x[,-1] = apply(x[,-1], 2, function(s) (s - mean(s))/sd(s))
#' z[,-1] = apply(z[,-1], 2, function(s) (s - mean(s))/sd(s))
#'
#' p = ncol(x)
#' q = ncol(z)
#'
#' # covariates
#' beta0 = 1
#' beta1 = 1.5
#' beta2 = 1e-4
#' beta3 = 1.4
#' beta4 = 1.1
#' beta5 = 1e-4
#' beta6 = 2.5
#' beta.true = c(beta0, beta1, beta2, beta3, beta4, beta5, beta6)
#' mu_sim.sp = exp(x %*% beta.true + w.true[index])
#' gamma0 = 1
#' gamma1 = 1e-4
#' gamma2 = 1.5
#' gamma3 = 1.1
#' gamma4 = 1e-4
#' gamma5 = -2.5
#' gamma6 = 1e-4
#' gamma.true = c(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)
#' phi_sim = exp(z %*% gamma.true)
#'
#' range(mu_sim.sp)
#' range(phi_sim)
#'
#' xi.true = 1.5
#'
#' y_sim = rtweedie(N, xi = xi.true, mu = mu_sim.sp, phi = phi_sim)
#' sum(y_sim == 0)/N # proportion of zeros
#' y = y_sim
#' x = x
#' z = z
#' # mcmc parameters
#' niter = 1e4
#' nburn = niter/2
#' report = 1e2
#'
#' # hyperparameters
#' lower.xi = 1
#' upper.xi = 2
#' system.time(mc_ss <- spssdglm.autograd(coords = coords, y = y, x = x, z = z,
#'                                     niter = niter, nburn= nburn, report = report, thin = 50,
#'                                     index = index, lower.xi = lower.xi, upper.xi = upper.xi,
#'                                     verbose = TRUE))
#' # Check for convergence
#' plot_mcmc(samples = mc_ss$xi.mcmc, true = 1.5, col = "blue", cnames = "xi")
##############################################################
# Hierarchical Bayesian Tweedie Compound Poisson Gamma       #
# Double Generalized Linear Model with Spike and Slab Priors #
##############################################################
# to-dos: right now just fits with exponential covariance
# generalize this to Matern kernel
spssdglm.autograd <- function(coords = NULL,
                              y = NULL,
                              x = NULL,
                              z = NULL,
                              index = NULL,
                              beta.init = NULL,
                              w.init = NULL,
                              gamma.init = NULL,
                              xi.init = NULL,
                              alpha.beta.init = NULL,
                              alpha.gamma.init = NULL,
                              zeta.beta.init = NULL,
                              zeta.gamma.init = NULL,
                              sigma2.beta.init = NULL,
                              sigma2.gamma.init = NULL,
                              shape.sigma2.beta = NULL,
                              shape.sigma2.gamma = NULL,
                              rate.sigma2.beta = NULL,
                              rate.sigma2.gamma = NULL,
                              sigma2.init = NULL,
                              phis.init = NULL,
                              lower.xi = NULL,
                              upper.xi = NULL,
                              shape.sigma2 = NULL,
                              scale.sigma2 = NULL,
                              lower.phis = NULL,
                              upper.phis = NULL,
                              nu0 = 1e-4,
                              tau.beta = 1e-1,
                              tau.gamma = 1e-1,
                              tau.xi = 1e-1,
                              tau.phis = 1e-1,
                              niter = NULL,
                              nburn = NULL,
                              report = NULL,
                              thin = 1,
                              return.mcmc = TRUE,
                              verbose = FALSE,
                              track = FALSE,
                              digits = 3,
                              reg.factor = 0,
                              reg.factor.sp = 0){
  # distance matrix
  Delta = as.matrix(dist(coords))

  N = length(y)
  L = nrow(coords)
  p = ncol(x)
  q = ncol(z)

  flag.int = 0
  if(sum(x[,1] == 1) == N){
    flag.int = 1
    # print(" Intercept being estimated under linear constraints x . 1 = 0 and w . 1 = 0! ")
    # print(" Design matrices x and z are being scaled.. ")
    # x[,-1] = apply(x[,-1], 2, function(s) (s - mean(s))/sd(s))
    # z[,-1] = apply(z[,-1], 2, function(s) (s - mean(s))/sd(s))
  }else{
    print(" Fitting mean model without intercept.. ")
  }

  # to-do:: construct index for replication automatically
  if(N > L){
    if(is.null(index)) stop(" Error: index missing! ")
  }else if(N == L){
    if(is.null(index)) index = 1:L
  }else{
    stop(" Error: Computationally singular system -- number of observations (N) less than number of locations (L) -- ")
  }

  # some housekeeping
  if(is.null(x) | is.null(z)) stop(" Error: either mean (x) or dispersion (z) (or both) design matrices missing! ")
  if(is.null(lower.xi) | is.null(upper.xi)) stop(" Error: specify prior for index parameter in Tweedie! ")
  if(is.null(niter)){
    warning(" Warning: chain length not specified, setting defaults to length = 1e4, burnin = 5e3, report = 1e2. ")
    niter = 1e4
    nburn = niter/2
    report = 1e2
  }
  if(is.null(nburn) & !is.null(niter)){
    warning(" Warning: burn-in not specified, setting default to niter/2. ")
    nburn = niter/2
  }
  if(is.null(report)){
    warning(" Warning: batch length not specified, setting default to 100. ")
    report = 1e2
  }

  if(is.null(shape.sigma2)) shape.sigma2 = 2
  if(is.null(scale.sigma2)) scale.sigma2 = 1
  if(is.null(lower.phis)) lower.phis = 0
  if(is.null(upper.phis)) upper.phis = 60

  # acceptance
  accept.beta = accept.gamma = accept.xi = accept.phis = 0
  accept.mat = c()

  # storage
  res_beta = res_zeta.beta = res_sigma2.beta = matrix(0, nrow = niter, ncol = p)
  res_w = matrix(0, nrow = niter, ncol = L)
  res_gamma = res_zeta.gamma = res_sigma2.gamma = matrix(0, nrow = niter, ncol = q)
  res_xi = res_sigma2 = res_phis = res_alpha.beta = res_alpha.gamma = rep(0, niter)

  # default hyperparameters
  shape.sigma2.beta = ifelse(is.null(shape.sigma2.beta), 2, shape.sigma2.beta)
  shape.sigma2.gamma = ifelse(is.null(shape.sigma2.gamma), 2, shape.sigma2.gamma)
  rate.sigma2.beta = ifelse(is.null(rate.sigma2.beta), 1, rate.sigma2.beta)
  rate.sigma2.gamma = ifelse(is.null(rate.sigma2.gamma), 1, rate.sigma2.gamma)

  if(is.null(beta.init) | is.null(w.init)){
    betaw.init = solve(rbind(cbind(crossprod(x, x),
                                   t(apply(x, 2, function(s) aggregate(s, list(index), sum)[, 2]))),
                             cbind(apply(x, 2, function(s) aggregate(s, list(index), sum)[, 2]),
                                   diag(as.vector(table(index))))) + 1e-3 * diag(p + L)) %*%
      matrix(c(crossprod(x,log(y + 1e-1)), aggregate(log(y + 1e-1), list(index), sum)[, 2]), ncol = 1)
    beta = betaw.init[1 : p]
    w = betaw.init[-(1 : p)]
  }else{
    beta = beta.init
    w = w.init
  } # starting at MLE
  if(is.null(gamma.init)){
    gamma = solve(crossprod(z, z)) %*% crossprod(z, log(y + 1e-3))
  }else{
    gamma = gamma.init
  } # starting at MLE
  if(is.null(sigma2.init)) sigma2 = 5 else sigma2 = sigma2.init
  if(is.null(phis.init)) phis = 1 else phis = phis.init
  xi = ifelse(is.null(xi.init), lower.xi + (upper.xi - lower.xi)/2, xi.init)
  alpha.beta = ifelse(is.null(alpha.beta.init), 0.5, alpha.beta.init)
  alpha.gamma = ifelse(is.null(alpha.gamma.init), 0.5, alpha.gamma.init)
  if(is.null(zeta.beta.init)) zeta.beta = rep(1, p) else zeta.beta = zeta.beta.init
  if(is.null(zeta.gamma.init)) zeta.gamma = rep(1, q) else zeta.gamma = zeta.gamma.init
  if(is.null(sigma2.beta.init)) sigma2.beta = rep(1, p) else sigma2.beta = sigma2.beta.init
  if(is.null(sigma2.gamma.init)) sigma2.gamma = rep(1, q) else sigma2.gamma = sigma2.gamma.init

  # target density
  trgt_dfn = rep(0, niter)

  # to-do:: generalize to Matern family
  R = exp(- phis * Delta) # exponential
  chol.R = chol(R)
  R.inv = chol2inv(chol.R)
  E = diag(chol.R)

  ztz = crossprod(z, z)
  xb = as.vector(x %*% beta + w[index])
  zg = as.vector(z %*% gamma)

  A.gamma.mat = ztz + (reg.factor + 1) * diag(q)
  A.gamma.chol = chol(A.gamma.mat)
  A.gamma = as.matrix(chol2inv(A.gamma.chol))
  # browser()
  for(i in 1:niter){
    #####################
    # update beta and w #
    #####################
    # ** intercept (beta0) is estimated **
    if(flag.int == 1){
      t1 = exp((2 - xi) * xb - zg)
      t2 = y * exp((1 - xi) * xb - zg)

      A.a = x[,-1] * as.vector(t1)
      A.beta.mat.11 = crossprod(A.a, x[, -1]) + reg.factor * diag((p - 1)) + diag(1/(sigma2.beta[-1] * zeta.beta[-1]))
      A.beta.mat.21 = apply(A.a, 2, function(s) aggregate(s, list(index), sum)[, 2])
      A.beta.mat.22 = diag(aggregate(as.vector(t1), list(index), sum)[, 2]) + reg.factor.sp * diag(L) + R.inv/sigma2
      A.beta.mat = rbind(cbind(A.beta.mat.11, t(A.beta.mat.21)), cbind(A.beta.mat.21, A.beta.mat.22))
      A.beta.chol = chol(A.beta.mat)
      A.beta = as.matrix(chol2inv(A.beta.chol))
      nabla.beta = - matrix(c(beta[-1]/(sigma2.beta[-1] * zeta.beta[-1]),
                              crossprod(R.inv, w)/sigma2),
                            ncol = 1) -
        matrix(c(crossprod(x[, -1], (t1 - t2)),
                 aggregate((t1 - t2), list(index), sum)[, 2]),
               ncol = 1)

      beta.draw =  matrix(c(beta[-1], w), ncol = 1) +
        tau.beta * A.beta %*% nabla.beta/2 +
        chol(tau.beta * A.beta) %*% rnorm((p + L - 1))
      # check this
      beta.n0.draw = beta.draw[1 : (p - 1)]
      w.draw = beta.draw[-(1 : (p - 1))] - mean(beta.draw[-(1 : (p - 1))][index])
      xb.draw = as.vector((x[, -1] %*% beta.n0.draw) + w.draw[index])


      t1.draw = exp((2 - xi) * xb.draw - zg)
      t2.draw = y * exp((1 - xi) * xb.draw - zg)

      # update beta0
      A.beta0 = 1/(sum((2 - xi) * t1.draw * exp((2 - xi) * beta[1]) - (1 - xi) * t2.draw * exp((1 - xi) * beta[1])) + 1/(sigma2.beta[1] * zeta.beta[1]))
      nabla.beta0 = - beta[1]/(sigma2.beta[1] * zeta.beta[1]) - sum(t1.draw * exp((2 - xi) * beta[1]) - t2.draw * exp((1 - xi) * beta[1]))

      beta0.draw = beta[1] + tau.beta * A.beta0 * nabla.beta0/2 + sqrt(tau.beta * A.beta0) * rnorm(1)
      xb.draw = xb.draw + beta0.draw

      nabla.beta0.draw = - beta0.draw/(sigma2.beta[1] * zeta.beta[1]) - sum(t1.draw * exp((2 - xi) * beta0.draw) - t2.draw * exp((1 - xi) * beta0.draw))

      t1.draw = t1.draw * exp((2 - xi) * beta0.draw)
      t2.draw = t2.draw * exp((1 - xi) * beta0.draw)
      nabla.beta.draw = - matrix(c(beta.n0.draw/(sigma2.beta[-1] * zeta.beta[-1]),
                                   crossprod(R.inv, w.draw)/sigma2),
                                 ncol = 1) -
        matrix(c(crossprod(x[, -1], (t1.draw - t2.draw)),
                 aggregate(t1.draw - t2.draw, list(index), sum)[, 2]),
               ncol = 1)

      ra1 = sum(dev.tw(y = y, mu = exp(xb.draw), phi = exp(zg), xi = xi)) - sum(dev.tw(y = y, mu = exp(xb), phi = exp(zg), xi = xi))
      ra2 = 0.5 * (sum(c(beta0.draw, beta.n0.draw)^2/(sigma2.beta * zeta.beta)) - sum(beta^2/(sigma2.beta * zeta.beta))) + 0.5 * (crossprod(t(crossprod(w.draw, R.inv)), w.draw) - crossprod(t(crossprod(w, R.inv)), w))/sigma2
      ra = ra1 - ra2
      if(is.nan(ra)) ra = -Inf # error checks

      rb1.1 = matrix(c(beta.n0.draw, w.draw), ncol= 1) - matrix(c(beta[-1], w), ncol = 1) - tau.beta * A.beta %*% nabla.beta/2
      rb1 = -1/(2 * tau.beta) * crossprod(t(crossprod(rb1.1, A.beta.mat)), (rb1.1))
      rb2.1 = matrix(c(beta[-1], w), ncol = 1) - matrix(c(beta.n0.draw, w.draw), ncol= 1) - tau.beta * A.beta %*% nabla.beta.draw/2
      rb2 = -1/(2 * tau.beta) * crossprod(t(crossprod(rb2.1, A.beta.mat)), (rb2.1))
      rb = rb2 - rb1
      if(is.nan(rb)) rb = -Inf # error checks

      rc1.1 = beta0.draw - beta[1] - tau.beta * A.beta0 * nabla.beta0/2
      rc1 = - rc1.1^2/(2 * tau.beta * A.beta0)
      rc2.1 = beta[1] - beta0.draw - tau.beta * A.beta0 * nabla.beta0.draw/2
      rc2 = - rc2.1^2/(2 * tau.beta * A.beta0)
      rc = rc2 - rc1
      if(is.nan(rc)) rc = -Inf # error checks

      accept.prob = min(ra + rb + rc, 0)
      if(log(runif(1)) < accept.prob){
        res_beta[i,] = beta = c(beta0.draw,  as.vector(beta.n0.draw))
        res_w[i,] = w = as.vector(w.draw)
        xb = xb.draw
        accept.beta = accept.beta + 1
      }else{
        res_beta[i,] = beta
        res_w[i,] = w
      }
    }else{
      t1 = exp((2 - xi) * xb - zg)
      t2 = y * exp((1 - xi) * xb - zg)

      A.a = x * as.vector(t1)
      A.beta.mat.11 = crossprod(A.a, x) + reg.factor * diag(p) + diag(1/(sigma2.beta * zeta.beta))
      A.beta.mat.21 = apply(A.a, 2, function(s) aggregate(s, list(index), sum)[, 2])
      A.beta.mat.22 = diag(aggregate(as.vector(t1), list(index), sum)[, 2]) + reg.factor.sp * diag(L) + R.inv/sigma2
      A.beta.mat = rbind(cbind(A.beta.mat.11, t(A.beta.mat.21)), cbind(A.beta.mat.21, A.beta.mat.22))
      A.beta.chol = chol(A.beta.mat)
      A.beta = as.matrix(chol2inv(A.beta.chol))
      nabla.beta = - matrix(c(beta/(sigma2.beta * zeta.beta),
                              crossprod(R.inv, w)/sigma2),
                            ncol = 1) -
        matrix(c(crossprod(x, (t1 - t2))[,1],
                 aggregate((t1 - t2), list(index), sum)[, 2]),
               ncol = 1)

      betaw.draw =  matrix(c(beta, w), ncol = 1) +
        tau.beta * A.beta %*% nabla.beta/2 +
        chol(tau.beta * A.beta) %*% rnorm((p + L))
      beta.draw = betaw.draw[1 : p]
      w.draw = betaw.draw[-(1 : p)]
      xb.draw = as.vector((x %*% beta.draw) + (w.draw)[index])

      t1.draw = exp((2 - xi) * xb.draw - zg)
      t2.draw = y * exp((1 - xi) * xb.draw - zg)
      nabla.beta.draw = - matrix(c(beta.draw/(sigma2.beta * zeta.beta),
                                   crossprod(R.inv, w.draw)/sigma2),
                                 ncol = 1) -
        matrix(c(crossprod(x, (t1.draw - t2.draw))[,1],
                 aggregate(t1.draw - t2.draw, list(index), sum)[, 2]),
               ncol = 1)

      ra1 = sum(dev.tw(y = y, mu = exp(xb.draw), phi = exp(zg), xi = xi)) - sum(dev.tw(y = y, mu = exp(xb), phi = exp(zg), xi = xi))
      ra2 = 0.5 * (sum(beta.draw^2/(sigma2.beta * zeta.beta)) - sum(beta^2/(sigma2.beta * zeta.beta))) + 0.5 * (crossprod(t(crossprod(w.draw, R.inv)), w.draw) - crossprod(t(crossprod(w, R.inv)), w))/sigma2
      ra = ra1 - ra2
      if(is.nan(ra)) ra = -Inf # error checks

      rb1.1 = betaw.draw - matrix(c(beta, w), ncol = 1) - tau.beta * A.beta %*% nabla.beta/2
      rb1 = -1/(2 * tau.beta) * crossprod(t(crossprod(rb1.1, A.beta.mat)), (rb1.1))
      rb2.1 = matrix(c(beta, w), ncol = 1) - betaw.draw - tau.beta * A.beta %*% nabla.beta.draw/2
      rb2 = -1/(2 * tau.beta) * crossprod(t(crossprod(rb2.1, A.beta.mat)), (rb2.1))
      rb = rb2 - rb1
      if(is.nan(rb)) rb = -Inf # error checks

      accept.prob = min(as.numeric(ra + rb), 0)
      if(log(runif(1)) < accept.prob){
        res_beta[i,] = beta = as.vector(beta.draw)
        res_w[i,] = w = as.vector(w.draw)
        xb = xb.draw
        accept.beta = accept.beta + 1
      }else{
        res_beta[i,] = beta
        res_w[i,] = w
      }
    }

    ##############################
    # update spike-slab for beta #
    ##############################
    # update zeta.beta
    tmp = exp(- 0.5 * beta^2/sigma2.beta)
    t1 = (1 - alpha.beta) * 1/sqrt(nu0) * tmp^(1/nu0)
    t2 = alpha.beta * tmp
    zeta.pr = t2/(t1 + t2)
    zeta.pr[is.na(zeta.pr)] = 0.5
    res_zeta.beta[i,] = zeta.beta = rbinom(p, 1, prob = zeta.pr)
    zeta.beta[zeta.beta == 0] = nu0

    # update sigma2.beta
    res_sigma2.beta[i,] = sigma2.beta = 1/rgamma(p,
                                                 shape = shape.sigma2.beta + 0.5,
                                                 rate = rate.sigma2.beta + 0.5 * beta^2/zeta.beta)

    # update alpha.beta
    res_alpha.beta[i] = alpha.beta = rbeta(1, 1 + sum(zeta.beta == 1), 1 + sum(zeta.beta == nu0))

    ################
    # update gamma #
    ################
    dldphi = delphi(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    nabla.gamma = - gamma/(sigma2.gamma * zeta.gamma) + crossprod(z, exp(zg) * dldphi)

    gamma.draw =  gamma + tau.gamma * A.gamma %*% nabla.gamma/2 + chol(tau.gamma * A.gamma) %*% rnorm(q)
    zg.draw = as.vector(z %*% gamma.draw)

    dldphi.draw = delphi(y = y, xi = xi, mu = exp(xb), phi = exp(zg.draw))
    nabla.gamma.draw = - gamma.draw/(sigma2.gamma * zeta.gamma) +  crossprod(z, exp(zg) * dldphi.draw)

    ra1.tmp1 = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg.draw))
    ra1.tmp1[ra1.tmp1 == 0] = 1e-300
    ra1.tmp2 = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    ra1.tmp2[ra1.tmp2 == 0] = 1e-300
    ra1 = sum(log(ra1.tmp1/ra1.tmp2))
    ra2 = 0.5 * (sum(gamma.draw^2/(sigma2.gamma * zeta.gamma)) - sum(gamma^2/(sigma2.gamma * zeta.gamma)))
    ra = ra1 - ra2
    if(is.nan(ra)) ra = -Inf # error checks

    rb1.1 = gamma.draw - gamma - tau.gamma * A.gamma %*% nabla.gamma/2
    rb1 = -1/(2 * tau.gamma) * crossprod(t(crossprod(rb1.1, A.gamma.mat)), (rb1.1))
    rb2.1 = gamma - gamma.draw - tau.gamma * A.gamma %*% nabla.gamma.draw/2
    rb2 = -1/(2 * tau.gamma) * crossprod(t(crossprod(rb2.1, A.gamma.mat)), (rb2.1))
    rb = rb2 - rb1
    if(is.nan(rb)) rb = -Inf # error checks

    accept.prob = min(ra + rb, 0)
    if(log(runif(1)) < accept.prob){
      res_gamma[i,] = gamma = as.vector(gamma.draw)
      zg = zg.draw
      accept.gamma = accept.gamma + 1
    }else{
      res_gamma[i,] = gamma
    }

    ###############################
    # update spike-slab for gamma #
    ###############################
    # update zeta.gamma
    tmp = exp(- 0.5 * gamma^2/sigma2.gamma)
    t1 = (1 - alpha.gamma) * 1/sqrt(nu0) * tmp^(1/nu0)
    t2 = alpha.gamma * tmp
    zeta.pr = t2/(t1 + t2)
    zeta.pr[is.na(zeta.pr)] = 0.5
    res_zeta.gamma[i,] = zeta.gamma = rbinom(q, 1, prob = zeta.pr)
    zeta.gamma[zeta.gamma == 0] = nu0

    # update sigma2.gamma
    res_sigma2.gamma[i,] = sigma2.gamma = 1/rgamma(q,
                                                   shape = shape.sigma2.gamma + 0.5,
                                                   rate = rate.sigma2.gamma + 0.5 * gamma^2/zeta.gamma)

    # update alpha.gamma
    res_alpha.gamma[i] = alpha.gamma = rbeta(1, 1 + sum(zeta.gamma == 1), 1 + sum(zeta.gamma == nu0))

    #############
    # update xi #
    #############
    xi.draw = tau.xi * rnorm(1) + xi

    if(xi.draw < lower.xi | xi.draw > upper.xi){
      accept.prob = -Inf
    }else{
      ra1.tmp1 = dtweedie(y = y, xi = xi.draw, mu = exp(xb), phi = exp(zg))
      ra1.tmp1[ra1.tmp1 == 0] = 1e-300
      ra1.tmp2 = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
      ra1.tmp2[ra1.tmp2 == 0] = 1e-300
      ra = sum(log(ra1.tmp1/ra1.tmp2))
      accept.prob = min(ra, 0)
    }
    if(log(runif(1)) < accept.prob){
      res_xi[i] = xi = xi.draw
      accept.xi = accept.xi + 1
    }else{
      res_xi[i] = xi
    }
    #################
    # update sigma2 #
    #################
    res_sigma2[i] = sigma2 = 1/rgamma(1, shape = shape.sigma2 + L/2,
                                      rate = 1/scale.sigma2 +
                                        crossprod(t(crossprod(w, R.inv)), w)/2)

    ###############
    # update phis #
    ###############
    phis.draw = tau.phis * rnorm(1) + phis
    if(phis.draw < lower.phis | phis.draw > upper.phis){
      accept.prob = -Inf
    }else{
      R.draw_s = exp(-phis.draw*Delta)
      chol.R.draw_s = chol(R.draw_s)
      R.draw.inv_s = chol2inv(chol.R.draw_s)
      E.draw_s = diag(chol.R.draw_s)

      ra = sum(log((E/E.draw_s)))
      rb = - 0.5 * crossprod(t(crossprod(w, (R.draw.inv_s - R.inv))), w)/sigma2
      accept.prob = min(ra + rb, 0)
    }

    if(log(runif(1)) < accept.prob){
      res_phis[i] = phis = phis.draw
      R = R.draw_s
      R.inv = R.draw.inv_s
      E = E.draw_s
      accept.phis = accept.phis + 1
    }else{
      res_phis[i] = phis
    }

    # target density evaluation
    lik.val = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    lik.val[lik.val == 0] = 1e-300
    trgt_dfn[i] = sum(log(lik.val)) +
      mvtnorm::dmvnorm(as.vector(beta), rep(0, p), diag((sigma2.beta * zeta.beta)), log = T) +
      mvtnorm::dmvnorm(as.vector(w), rep(0, L), sigma2 * R, log = T) +
      mvtnorm::dmvnorm(as.vector(gamma), rep(0, q), diag((sigma2.gamma * zeta.gamma)), log = T) +
      dgamma(1/sigma2, shape = shape.sigma2, scale = scale.sigma2) + 0 + 0
    if(track){
      if(i %% report == 0) cat(round(trgt_dfn[i], 2), "\n") else cat(round(trgt_dfn[i], 2), "\t")
    }
    #############################
    # proposal variance scaling #
    #############################
    if(i %% report == 0){
      accept = c(accept.beta, accept.gamma, accept.xi, accept.phis)
      accept.p = accept/report
      accept.mat = rbind(accept.mat, accept.p)

      if(verbose){
        if(i>nburn){
          cat("Iteration::", i, "Acceptance:", accept.p,"\n",
              "=============================================","\n",
              "Beta::", round(apply(res_beta[nburn:i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(tau.beta, 2),"\n",
              "Gamma::", round(apply(res_gamma[nburn:i,], 2, median), digits = (digits - 1)), ", tuning.gamma = ", round(tau.gamma, 2),"\n",
              "Sigma::", round(median(res_sigma2[nburn:i]), digits = (digits - 1)), "\n",
              "Phis::", round(median(res_phis[nburn:i]), digits = (digits - 1)), ", tuning.phis = ", round(tau.phis, 2),"\n",
              "Xi::", round(median(res_xi[nburn:i]), digits = (digits - 1)), ", tuning.xi = ", round(tau.xi, 2),"\n",
              "=============================================","\n")
        }else{
          if(track){
            # cat("Iteration::", i, "Acceptance:", accept.p, "\n")
            cat("Iteration::", i, "Acceptance:", accept.p,"\n",
                "=============================================","\n",
                "Beta::", round(apply(res_beta[(i-report+1):i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(tau.beta, 2),"\n",
                "Gamma::", round(apply(res_gamma[(i-report+1):i,], 2, median), digits = (digits - 1)), ", tuning.gamma = ", round(tau.gamma, 2),"\n",
                "Sigma::", round(median(res_sigma2[(i-report+1):i]), digits = (digits - 1)), "\n",
                "Phis::", round(median(res_phis[(i-report+1):i]), digits = (digits - 1)), ", tuning.phis = ", round(tau.phis, 2),"\n",
                "Xi::", round(median(res_xi[(i-report+1):i]), digits = (digits - 1)), ", tuning.xi = ", round(tau.xi, 2),"\n",
                "=============================================","\n")
          }else{
            cat("Iteration::", i, "Acceptance:", accept.p, "\n")
          }
        }
      }

      # MALA: optimal scales (58%)
      step = sapply(accept.p[1:2], function(x){
        out = 1
        x = max(0.25, min(x, 0.75))
        if(x > 0.65) out = x/0.65
        else if(x < 0.45) out = x/0.45
        out
      })
      tau.beta = tau.beta * step[1]
      tau.gamma = tau.gamma * step[2]

      # RW: optimal scales (33%)
      step = sapply(accept.p[3:4], function(x){
        out = 1
        x = max(0.17, min(x, 0.75))
        if(x > 0.5) out = x/0.5
        else if(x < 0.25) out = x/0.25
        out
      })
      tau.xi = tau.xi * step[1]
      tau.phis = tau.phis * step[2]

      accept.beta = accept.gamma = accept.xi = accept.phis = 0
    }
  }

  #############
  # Inference #
  #############
  sample.win = seq((nburn + 1), niter, by = thin)
  beta.est = cbind.data.frame(MAP = round(apply(res_beta[sample.win,], 2, function(s){ den = density(s); ind = which.max(den$y); den$x[ind]}), digits = digits),
                              median = round(apply(res_beta[sample.win,], 2, median), digits = digits),
                              mean = round(apply(res_beta[sample.win,], 2, mean), digits = digits),
                              sd = round(apply(res_beta[sample.win,], 2, sd), digits = digits),
                              round(coda::HPDinterval(coda::as.mcmc(res_beta[sample.win,])), digits = digits))
  gamma.est = cbind.data.frame(MAP = round(apply(res_gamma[sample.win,], 2, function(s){ den = density(s); ind = which.max(den$y); den$x[ind]}), digits = digits),
                               median = round(apply(res_gamma[sample.win,], 2, median), digits = digits),
                               mean = round(apply(res_gamma[sample.win,], 2, mean), digits = digits),
                               sd = round(apply(res_gamma[sample.win,], 2, sd), digits = digits),
                               round(coda::HPDinterval(coda::as.mcmc(res_gamma[sample.win,])), digits = digits))
  xi.est = cbind.data.frame(MAP = round(apply(matrix(res_xi[sample.win], ncol = 1), 2, function(s){ den = density(s); ind = which.max(den$y); den$x[ind]}), digits = digits),
                            median = round(median(res_xi[sample.win]), digits = digits),
                            mean = round(mean(res_xi[sample.win]), digits = digits),
                            sd = round(sd(res_xi[sample.win]), digits = digits),
                            round(coda::HPDinterval(coda::as.mcmc(res_xi[sample.win])), digits = digits))
  sigma2.est = cbind.data.frame(median = round(median(res_sigma2[sample.win]), digits = digits),
                                mean = round(mean(res_sigma2[sample.win]), digits = digits),
                                sd = round(sd(res_sigma2[sample.win]), digits = digits),
                                round(coda::HPDinterval(coda::as.mcmc(res_sigma2[sample.win])), digits = digits))
  phis.est = cbind.data.frame(median = round(median(res_phis[sample.win]), digits = digits),
                              mean = round(mean(res_phis[sample.win]), digits = digits),
                              sd = round(sd(res_phis[sample.win]), digits = digits),
                              round(coda::HPDinterval(coda::as.mcmc(res_phis[sample.win])), digits = digits))
  w.est = cbind.data.frame(median = round(apply(res_w[sample.win,], 2, median), digits = digits),
                           mean = round(apply(res_w[sample.win,], 2, mean), digits = digits),
                           sd = round(apply(res_w[sample.win,], 2, sd), digits = digits),
                           round(coda::HPDinterval(coda::as.mcmc(res_w[sample.win,])), digits = digits))

  model.summary.fixed = rbind.data.frame(beta.est, gamma.est, xi.est)
  model.summary.fixed$sig = apply(model.summary.fixed, 1, function(x){
    if(x[5] > 0 | x[6] < 0) return("*")
    else if(x[5] < 0 & x[6] > 0) return("")
    else if(x[5] == 0 | x[6] == 0) return("")
  })
  if(flag.int == 1){
    rownames(model.summary.fixed) = c(paste("beta", 0 : (p - 1), sep=""),
                                      paste("gamma", 0 : (q - 1), sep=""),
                                      "xi")
  }else{
    rownames(model.summary.fixed) = c(paste("beta", 1 : p, sep=""),
                                      paste("gamma", 0 : (q - 1), sep=""),
                                      "xi")
  }
  model.summary.latent = rbind.data.frame(sigma2.est, phis.est, w.est)
  model.summary.latent$sig = apply(model.summary.latent, 1, function(x){
    if(x[4] > 0 | x[5] < 0) return("*")
    else if(x[4] < 0 & x[5] > 0) return("")
    else if(x[4] == 0 | x[5] == 0) return("")
  })
  rownames(model.summary.latent) = c("sigma2", "phis", paste("loc", 1:L, sep=""))

  if(flag.int == 1){
    A.beta.mat = rbind(c(1/A.beta0, rep(0, p + L - 1)), cbind(0, A.beta.mat))
  }

  if(return.mcmc){
    if(flag.int == 1) colnames(res_beta) = paste("beta", 0 : (p - 1), sep="") else colnames(res_beta) = paste("beta", 1 : p, sep="")
    colnames(res_w) = paste("loc", 1:L, sep="")
    colnames(A.beta.mat) = rownames(A.beta.mat) = c(colnames(res_beta),colnames(res_w))
    colnames(A.gamma.mat) = rownames(A.gamma.mat) = colnames(res_gamma) = colnames(res_zeta.gamma) = paste("gamma", 0 : (q - 1), sep="")
    res_xi = matrix(res_xi, ncol = 1)
    colnames(res_xi) = "xi"
    res_alpha.beta = matrix(res_alpha.beta, ncol = 1)
    colnames(res_alpha.beta) = "alpha.beta"
    res_alpha.gamma = matrix(res_alpha.gamma, ncol = 1)
    colnames(res_alpha.gamma) = "alpha.gamma"
    res_sigma2 = matrix(res_sigma2, ncol = 1)
    colnames(res_sigma2) = "sigma2"
    res_phis = matrix(res_phis, ncol = 1)
    colnames(res_phis) = "phis"
    return(list(model.fixed = model.summary.fixed,
                model.latent = model.summary.latent,
                mean.cov.model = A.beta.mat,
                disp.cov.model = A.gamma.mat,
                target_fn = as.ts(trgt_dfn),
                niter = niter,
                nburn = nburn,
                report = report,
                beta.mcmc = coda::as.mcmc(res_beta[sample.win,]),
                gamma.mcmc = coda::as.mcmc(res_gamma[sample.win,]),
                xi.mcmc = coda::as.mcmc(res_xi[sample.win, "xi"]),
                zeta.beta.mcmc = coda::as.mcmc(res_zeta.beta[sample.win,]),
                zeta.gamma.mcmc = coda::as.mcmc(res_zeta.gamma[sample.win,]),
                sigma2.beta.mcmc = coda::as.mcmc(res_sigma2.beta[sample.win,]),
                sigma2.gamma.mcmc = coda::as.mcmc(res_sigma2.gamma[sample.win,]),
                alpha.beta.mcmc = coda::as.mcmc(res_alpha.beta[sample.win,"alpha.beta"]),
                alpha.gamma.mcmc = coda::as.mcmc(res_alpha.gamma[sample.win,"alpha.gamma"]),
                sigma2.mcmc = coda::as.mcmc(res_sigma2[sample.win, "sigma2"]),
                phis.mcmc = coda::as.mcmc(res_phis[sample.win, "phis"]),
                w.mcmc = coda::as.mcmc(res_w[sample.win,])))
  }else{
    if(flag.int == 1) colnames(res_beta) = paste("beta", 0 : (p - 1), sep="") else colnames(res_beta) = paste("beta", 1 : p, sep="")
    colnames(res_w) = paste("loc", 1:L, sep="")
    colnames(A.beta.mat) = rownames(A.beta.mat) = c(colnames(res_beta)[-1],colnames(res_w))
    colnames(A.gamma.mat) = rownames(A.gamma.mat) = paste("gamma", 0 : (q - 1), sep="")
    return(list(model.fixed = model.summary.fixed,
                model.latent = model.summary.latent,
                mean.cov.model = A.beta.mat,
                disp.cov.model = A.gamma.mat,
                target_fn = as.ts(trgt_dfn)))
  }
}
