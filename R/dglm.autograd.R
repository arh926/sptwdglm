#' Tweedie Compound Poisson-Gamma (CP-g) Double Generalized Linear Models
#'
#' Fits a Double Generalized Linear Model: \eqn{\log(\mu)=x^T\beta} and \eqn{\log(\phi)=z^T\gamma}. Parameters not listed below are optional.
#'
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
#' \dontrun{
#' set.seed(2022)
#' require(tweedie)
#' require(Matrix)
#' # require(dglm) # for fitting traditional dglms
#' # require(statmod)
#'
#' par(mfcol=c(1,1))
#' # Generate Data
#' N = 1e3
#' x = z = cbind(1, rnorm(N), rnorm(N), rnorm(N)) # Matrix()
#' p = ncol(x)
#' q = ncol(z)
#' # Covariates
#' beta0 = 1
#' beta1 = 1.5
#' beta2 = 1.1
#' beta3 = 1.4
#' beta.true = c(beta0, beta1, beta2, beta3)
#' mu_sim = as.vector(exp(x %*% beta.true))
#' gamma0 = 1
#' gamma1 = 0.5
#' gamma2 = 0.1
#' gamma3 = 1.1
#' gamma.true = c(gamma0, gamma1, gamma2, gamma3)
#' phi_sim =  as.vector(exp(z %*% gamma.true))
#'
#' # weights = rgamma(N, shape = 1, scale = 1) # exposure
#'
#' xi.true = 1.5
#'
#' y_sim = rtweedie(N, xi = xi.true,
#'  mu = mu_sim,
#'   phi = phi_sim); range(y_sim)
#' # y_sim = y_sim/weights; range(y_sim)
#' sum(y_sim == 0)/N # proportion of zeros
#'
#' # Traditional DGLM
#' # mdglm = try(dglm(y_sim~as.matrix(x[,-1]),
#' # ~as.matrix(z[,-1]),
#' # family=tweedie(link.power=0, var.power=1.5)))
#' # if(!("try-error" %in% class(mdglm))){
#'  #  mdglm.mean = mdglm$coefficients # mean model
#'   # mdglm.disp = mdglm$dispersion.fit$coefficients # dispersion model
#' # }
#' # Bayesian DGLM
#' y = y_sim
#' x = x
#' z = z
#' # MCMC parameters
#' niter = 1e4
#' nburn = niter/2
#' report = 1e2
#'
#' # Hyperparameters
#' prec.beta = prec.gamma = 1e-6
#' lower.xi = 1
#' upper.xi = 2
#' system.time(mc <- dglm.autograd(y=y,
#'  x=x,
#'   z=z,
#'    lower.xi = lower.xi,
#'     upper.xi = upper.xi,
#'      niter = niter,
#'       verbose=TRUE,
#'        thin = 20))
#' # Model summary
#' cbind(mc$model, true=c(beta.true,gamma.true,xi.true))
#' # --mean model variance covariance matrix
#' mc$mean.cov.model
#' # --dispersion model variance covariance matrix
#' mc$disp.cov.model
#'
#' # Check convergence
#' plot_mcmc(mc$xi.mcmc, true =  xi.true, cnames= "xi")
#' }
####################################
#  Tweedie Compound Poisson Gamma  #
#  Double Generalized Linear Model #
####################################

dglm.autograd <- function(y = NULL,
                          x = NULL,
                          z = NULL,
                          beta.init = NULL,
                          gamma.init = NULL,
                          xi.init = NULL,
                          prec.beta = NULL,
                          prec.gamma = NULL,
                          lower.xi = NULL,
                          upper.xi = NULL,
                          tau.beta = 1e-1,
                          tau.gamma = 1e-1,
                          tau.xi = 1e-1,
                          niter = NULL,
                          nburn = NULL,
                          report = NULL,
                          thin = 1,
                          return.mcmc = TRUE,
                          verbose = FALSE,
                          track = FALSE,
                          digits = 3,
                          reg.factor = 2e1){

  if(is.null(y)) stop(" Error: response missing! ")
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

  p = ncol(x)
  q = ncol(z)

  if(is.null(prec.beta)) prec.beta = 1e-6
  if(is.null(prec.gamma)) prec.gamma = 1e-6

  # acceptance
  accept.beta = accept.gamma = accept.xi = 0
  accept.mat = c()

  # storage
  res_beta = matrix(0, nrow = niter, ncol = p)
  res_gamma = matrix(0, nrow = niter, ncol = q)
  res_xi = rep(0, niter)

  if(is.null(beta.init)){
    beta = solve(crossprod(x, x)) %*% crossprod(x, log(y + 1e-3))
  }else{
    beta = beta.init
  } # starting at MLE
  if(is.null(gamma.init)){
    gamma = solve(crossprod(z, z)) %*% crossprod(z, log(y + 1e-3))
  }else{
    gamma = gamma.init
  } # starting at MLE
  xi = ifelse(is.null(xi.init), lower.xi + (upper.xi - lower.xi)/2, xi.init)

  # target density
  trgt_dfn = rep(0, niter)

  ztz = crossprod(z, z)
  A.gamma.mat = ztz + (reg.factor + prec.gamma) * diag(q)
  A.gamma.chol = chol(A.gamma.mat)
  A.gamma = chol2inv(A.gamma.chol)

  xb = as.vector(x %*% beta)
  zg = as.vector(z %*% gamma)

  for(i in 1:niter){
    ###############
    # update beta #
    ###############
    t1 = exp((2 - xi) * xb - zg)
    t2 = y * exp((1 - xi) * xb - zg)

    A.beta.mat = crossprod(x * as.vector(t1), x) + (reg.factor + prec.beta) * diag(p) # forceSymmetric()
    A.beta.chol =  chol(A.beta.mat)
    A.beta = chol2inv(A.beta.chol)
    nabla.beta = - prec.beta * beta - crossprod(x, (t1 - t2))

    beta.draw =  beta + tau.beta * A.beta %*% nabla.beta/2 + chol(tau.beta * A.beta) %*% rnorm(p)
    xb.draw = as.vector(x %*% beta.draw)

    t1.draw = exp((2 - xi) * xb.draw - zg)
    t2.draw = y * exp((1 - xi) * xb.draw - zg)
    nabla.beta.draw = - prec.beta * beta.draw - crossprod(x, (t1.draw - t2.draw))

    ra1.tmp1 = dtweedie(y = y, xi = xi, mu = exp(xb.draw), phi = exp(zg))
    ra1.tmp1[ra1.tmp1 == 0] = 1e-300
    ra1.tmp2 = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    ra1.tmp2[ra1.tmp2 == 0] = 1e-300
    ra1 = sum(log(ra1.tmp1/ra1.tmp2))
    ra2 = 0.5 * prec.beta * (sum(beta.draw^2) - sum(beta^2))
    ra = ra1 - ra2

    rb1.1 = beta.draw - beta - tau.beta * A.beta %*% nabla.beta/2
    rb1 = -1/(2 * tau.beta) * crossprod(t(crossprod(rb1.1, A.beta.mat)), (rb1.1))
    rb2.1 = beta - beta.draw - tau.beta * A.beta %*% nabla.beta.draw/2
    rb2 = -1/(2 * tau.beta) * crossprod(t(crossprod(rb2.1, A.beta.mat)), (rb2.1))
    rb = rb2 - rb1

    accept.prob = min(as.numeric(ra + rb), 0)
    if(log(runif(1)) < accept.prob){
      res_beta[i,] = beta = as.vector(beta.draw)
      xb = xb.draw
      accept.beta = accept.beta + 1
    }else{
      res_beta[i,] = beta
    }

    ################
    # update gamma #
    ################
    dldphi = delphi(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    nabla.gamma = - prec.gamma * gamma + crossprod(z, exp(zg) * dldphi)

    gamma.draw =  gamma + tau.gamma * A.gamma %*% nabla.gamma/2 + chol(tau.gamma * A.gamma) %*% rnorm(q)
    zg.draw = as.vector(z %*% gamma.draw)

    dldphi.draw = delphi(y = y, xi = xi, mu = exp(xb), phi = exp(zg.draw))
    nabla.gamma.draw = - prec.gamma * gamma.draw +  crossprod(z, exp(zg) * dldphi.draw)

    ra1.tmp1 = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg.draw))
    ra1.tmp1[ra1.tmp1 == 0] = 1e-300
    ra1.tmp2 = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    ra1.tmp2[ra1.tmp2 == 0] = 1e-300
    ra1 = sum(log(ra1.tmp1/ra1.tmp2))
    ra2 = 0.5 * prec.gamma * (sum(gamma.draw^2) - sum(gamma^2))
    ra = ra1 - ra2
    rb1.1 = gamma.draw - gamma - tau.gamma * A.gamma %*% nabla.gamma/2
    rb1 = -1/(2 * tau.gamma) * crossprod(t(crossprod(rb1.1, A.gamma.mat)), (rb1.1))
    rb2.1 = gamma - gamma.draw - tau.gamma * A.gamma %*% nabla.gamma.draw/2
    rb2 = -1/(2 * tau.gamma) * crossprod(t(crossprod(rb2.1, A.gamma.mat)), (rb2.1))
    rb = rb2 - rb1

    accept.prob = min(as.numeric(ra + rb), 0)
    if(log(runif(1)) < accept.prob){
      res_gamma[i,] = gamma = as.vector(gamma.draw)
      zg = zg.draw
      accept.gamma = accept.gamma + 1
    }else{
      res_gamma[i,] = gamma
    }
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

    # target density evaluation
    lik.val = dtweedie(y = y, xi = xi, mu = exp(xb), phi = exp(zg))
    lik.val[lik.val == 0] = 1e-300
    trgt_dfn[i] = sum(log(lik.val)) +
      dmvnorm(as.vector(beta), rep(0, p), diag(p)/prec.beta, log = TRUE) +
      dmvnorm(as.vector(gamma), rep(0, q), diag(q)/prec.beta, log = TRUE) + 0
    if(track){
      if(i %% report == 0) cat(round(trgt_dfn[i], 2), "\n") else cat(round(trgt_dfn[i], 2), "\t")
    }
    #############################
    # proposal variance scaling #
    #############################
    if(i %% report == 0){
      accept = c(accept.beta, accept.gamma, accept.xi)
      accept.p = accept/report
      accept.mat = rbind(accept.mat, accept.p)

      if(verbose){
        if(i>nburn){
          cat("Iteration::", i, "Acceptance:", accept.p,"\n",
              "=============================================","\n",
              "Beta::", round(apply(res_beta[nburn:i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(tau.beta, 2),"\n",
              "Gamma::", round(apply(res_gamma[nburn:i,], 2, median), digits = (digits - 1)), ", tuning.gamma = ", round(tau.gamma, 2),"\n",
              "Xi::", round(median(res_xi[nburn:i]), digits = (digits - 1)), ", tuning.xi = ", round(tau.xi, 2),"\n",
              "=============================================","\n")
        }else{
          if(track){
            cat("Iteration::", i, "Acceptance:", accept.p,"\n",
                "=============================================","\n",
                "Beta::", round(apply(res_beta[(i-report+1):i,], 2, median), digits = (digits - 1)), ", tuning.beta = ", round(tau.beta, 2),"\n",
                "Gamma::", round(apply(res_gamma[(i-report+1):i,], 2, median), digits = (digits - 1)), ", tuning.gamma = ", round(tau.gamma, 2),"\n",
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
      accept.xi = max(0.17, min(accept.p[3], 0.75))
      if(accept.xi > 0.5) tau.xi = tau.xi * accept.xi/0.5
      else if(accept.xi < 0.25) tau.xi = tau.xi * accept.xi/0.25

      accept.beta = accept.gamma = accept.xi = 0
      gc()
    }
  }
  #############
  # Inference #
  #############
  sample.win = seq((nburn + 1), niter, by = thin)
  beta.est = cbind.data.frame(median = round(apply(res_beta[sample.win,], 2, median), digits = digits),
                              mean = round(apply(res_beta[sample.win,], 2, mean), digits = digits),
                              sd = round(apply(res_beta[sample.win,], 2, sd), digits = digits),
                              round(HPDinterval(as.mcmc(res_beta[sample.win,])), digits = digits))
  gamma.est = cbind.data.frame(median = round(apply(res_gamma[sample.win,], 2, median), digits = digits),
                               mean = round(apply(res_gamma[sample.win,], 2, mean), digits = digits),
                               sd = round(apply(res_gamma[sample.win,], 2, sd), digits = digits),
                               round(HPDinterval(as.mcmc(res_gamma[sample.win,])), digits = digits))
  xi.est = cbind.data.frame(median = round(median(res_xi[sample.win]), digits = digits),
                            mean = round(mean(res_xi[sample.win]), digits = digits),
                            sd = round(sd(res_xi[sample.win]), digits = digits),
                            round(HPDinterval(as.mcmc(res_xi[sample.win])), digits = digits))
  model.summary = rbind.data.frame(beta.est, gamma.est, xi.est)
  model.summary$sig = apply(model.summary, 1, function(x){
    if(x[4] > 0 | x[5] < 0) return("*")
    if(x[4] < 0 & x[5] > 0) return("")
  })
  rownames(model.summary) = c(paste("beta", 0 : (p - 1), sep=""),
                              paste("gamma", 0 : (q - 1), sep=""),
                              "xi")
  if(return.mcmc){
    colnames(A.beta.mat) = rownames(A.beta.mat) = colnames(res_beta) = paste("beta", 0 : (p - 1), sep="")
    colnames(A.gamma.mat) = rownames(A.gamma.mat) = colnames(res_gamma) = paste("gamma", 0 : (q - 1), sep="")
    res_xi = matrix(res_xi, ncol = 1)
    colnames(res_xi) = "xi"
    return(list(model = model.summary,
                mean.cov.model = A.beta.mat,
                disp.cov.model = A.gamma.mat,
                target_fn = as.ts(trgt_dfn),
                niter = niter,
                nburn = nburn,
                report = report,
                beta.mcmc = as.mcmc(res_beta[sample.win,]),
                gamma.mcmc = as.mcmc(res_gamma[sample.win,]),
                xi.mcmc = as.mcmc(res_xi[sample.win,"xi"])))
  }else{
    colnames(A.beta.mat) = rownames(A.beta.mat) =  paste("beta", 0 : (p - 1), sep="")
    colnames(A.gamma.mat) = rownames(A.gamma.mat) = paste("gamma", 0 : (q - 1), sep="")
    return(list(model = model.summary,
                mean.cov.model = A.beta.mat,
                disp.cov.model = A.gamma.mat,
                target_fn = as.ts(trgt_dfn)))
  }
}
