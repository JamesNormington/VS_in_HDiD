################################################################
##### Code to produce Tables 3-6 in "Bayesian variable selection 
##### in hierarchical difference-in-differences models"
##### Submitted to Biostatistics, September 2019
################################################################
rm(list = ls())

################################
### MAIN SIMULATION FUNCTION ###
################################
### Inputs ###
# NSIM: Number of simulations.
# nMCMC: Number of Gibbs draws using each variable selection method.
# BI: Number of initial Gibbs draws to discard as burn-in.
# J: Number of groups.
# beta.tilde.true: True values of beta.tilde (coefficients for baseline model) for data generation.
# beta.true: True values of beta (coefficients for change model, including Delta) for data generation.
# alpha.true: True values of alpha (coefficients for exposure model) for data generation.
# string: Character string to include in file name.
# seed: Set random seed for reproducibility.
### Outputs ###
# Write bias, MSE, coverage, and inclusion probabilities to separate text files.
run_sim = function(NSIM=5000, nMCMC=10000, BI=1000, J=50, beta.tilde.true, beta.true, alpha.true, string, seed) {
  setwd("/home/normington/Paper 2/PosteriorSimResults")
  require(mvtnorm); require(truncnorm)
  
  # Number of clinics and predictors.
  n.pre = n.post = rep(10, J)
  K = length(beta.true) - 1 # no variable selection on treatment
  
  # Initialize vectors and matrices to store MCMC draws.
  beta.tilde = rbind(beta.tilde.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
  beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
  mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
  mu.diff= rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
  sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
  gamma.tilde = gamma = c(5, rep(NA, nMCMC - 1))
  tau2.tilde = tau2 = gamma2.tilde = gamma2 = c(1, rep(NA, nMCMC-1))
  w.tilde =  rbind(rep(1, length(beta.tilde.true) - 1), matrix(NA, nrow = nMCMC-1, ncol = length(beta.tilde.true) - 1))
  w =  rbind(rep(1, length(beta.true) - 1), matrix(NA, nrow = nMCMC-1, ncol = length(beta.true) - 1))
  p.tilde = p = c(0.5, rep(NA, nMCMC - 1))
  a.tilde = a = rep(0.01, length(beta.true) - 1) 
  sigma2.tilde.b = sigma2.b = rep(NA, J)
  
  # Create lists to store individual-level outcome data.
  Y.pre = Y.post = list()
  for(j in 1:J) {
    Y.pre[[j]] = rep(NA, n.pre[j])
    Y.post[[j]] = rep(NA, n.post[j])
  }
  
  # Define post burn-in MCMC iterations.
  postBI = (BI+1):nMCMC
  
  # Set random seed.
  set.seed(seed)
  
  # Main simulation loop.
  for(nsim in 1:NSIM) {
    
    ##################
    ### Generate data.
    ##################
    X1 = rnorm(J); X2 = rnorm(J); X3 = rnorm(J); X4 = rnorm(J)
    X5 = rnorm(J); X6 = rnorm(J); X7 = rnorm(J); X8 = rnorm(J)
    X0 = cbind(X1, X2, X3, X4, X5, X6, X7, X8)
    Trt =  rnorm(J, X0%*%alpha.true)
    X = cbind(Trt, X1, X2, X3, X4, X5, X6, X7, X8)
    gram = t(X)%*%X
    gram0 = t(X0)%*%X0
    mu = rnorm(J, X%*%beta.tilde.true)
    mu.diff = rnorm(J, X%*%beta.true)
    for(j in 1:J) {
      Y.pre[[j]] = rnorm(n.pre[j], mu[j])
      Y.post[[j]] = rnorm(n.post[j], mu[j] + mu.diff[j])
    }
    
    Y.pre.means = unlist(lapply(Y.pre, mean))
    Y.post.means = unlist(lapply(Y.post, mean))
    
    ############
    ### Separate
    ############
    # Initialize vectors and matrices to store MCMC draws.
    if(TRUE) {
      beta.tilde = rbind(beta.tilde.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.tilde.true)))
      beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
      mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      mu.diff = rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
      tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
      
      # Separate indicators for baseline and change models.
      w.tilde = rbind(rep(1, length(beta.tilde.true) - 1), matrix(NA, nrow = nMCMC-1, ncol = K))
      w = rbind(rep(1, length(beta.true) - 1), matrix(NA, nrow = nMCMC-1, ncol = K))
      
      a.tilde = a = rep(0.1, length(beta.tilde.true) - 2)
      gamma.tilde = gamma = rbind(rep(1/25, K), matrix(NA, nrow = nMCMC - 1, ncol = K))
      sigma2.tilde.b = sigma2.b = rep(NA, J)
      n.vars = rep(NA, nMCMC-1)
      n.vars.tilde = rep(NA, nMCMC-1)
    }
    
    # Gibbs sampler.
    for(t in 2:nMCMC) {
      mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
      mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
      mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
      
      for(j in 1:J) {
        sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
        sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
      }
      
      sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
      sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
      
      tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X%*%beta.tilde[t-1,])^2))
      tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
      
      # Use the multi-variate normal representation of the spike-and-slab prior.
      a.tilde = ifelse(w.tilde[t-1,], 1/sqrt(gamma.tilde[t-1,]), 0.1)
      a = ifelse(w[t-1,], 1/sqrt(gamma[t-1,]), 0.1)
      
      # N(0, 100^2) prior on intercept.
      D.tilde = diag(c(100, a.tilde))
      D = diag(c(100, a))
      
      DD.tilde.inv = solve(D.tilde%*%D.tilde)
      DD.inv = solve(D%*%D)     
      
      V.tilde = solve(gram/tau2.tilde[t] + DD.tilde.inv)
      V = solve(gram/tau2[t] + DD.inv)
      
      beta.tilde[t,] = rmvnorm(1, (1/tau2.tilde[t])*V.tilde%*%t(X)%*%mu[t,], V.tilde)
      beta[t,] = rmvnorm(1, (1/tau2[t])*V%*%t(X)%*%mu.diff[t,], V)
      
      # Bayes Factors: P(beta | slab) / P(beta | spike)
      BF.tilde = exp(dnorm(beta.tilde[t,-1], sd = 0.1, log = TRUE) - dnorm(beta.tilde[t,-1], sd = 1/sqrt(gamma.tilde[t-1,]), log = TRUE))
      BF = exp(dnorm(beta[t,-1], sd = 0.1, log = TRUE) - dnorm(beta[t,-1], sd = 1/sqrt(gamma[t-1,]), log = TRUE))
      
      # P(w = 1 | beta) = P(beta | w = 1)P(w = 1) / [P(beta | w = 1)P(w = 1) + P(beta | w = 0)P(w = 0)] = 1 / (1 + BF)
      selp.tilde = 1 / (1 + BF.tilde)
      selp = 1 / (1 + BF)
      for(k in 1:K) {
        w.tilde[t,k] = sample(0:1, 1, prob = c(1-selp.tilde[k], selp.tilde[k]))
        w[t,k] = sample(0:1, 1, prob = c(1-selp[k], selp[k]))
      }
      
      # Number of variables selected in draw t.
      n.vars.tilde[t-1] = sum(w.tilde[t,])
      n.vars[t-1] = sum(w[t,])
      
      gamma.tilde[t,] = rgamma(K, shape = 2.5 + 0.5*w.tilde[t,], rate = 2.5*5^2 + 0.5*w.tilde[t,]*beta.tilde[t,-1]^2)
      gamma[t,] = rgamma(K, shape = 2.5 + 0.5*w[t,], rate = 2.5*5^2 + 0.5*w[t,]*beta[t,-1]^2)
      
    }
    
    ############################################
    ## Write iteration "n" results to text files
    ############################################
    
    # Compute posterior mean of Delta
    betaSeparate = colMeans(beta[postBI,])[1]
    
    # Compute inclusion probabilities for each covariate
    inclusionSeparate = colMeans(w[postBI,])
    
    # Write Bias for Delta 
    write.table(t(betaSeparate - beta.true[1]), file = paste("biasBetaSeparate", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(t(inclusionSeparate), file = paste("inclusionSeparate", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write MSE for Delta
    write.table(t((betaSeparate - beta.true[1])^2), file = paste("MSE_Separate", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write Coverage for Delta
    CI = quantile(beta[postBI,2], probs = c(.025, .975))
    write.table(ifelse(CI[1] < beta.true[1] & beta.true[1] < CI[2], 1, 0), file = paste("Coverage_Separate", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write number of covariates selected
    write.table(mean(n.vars[(BI+1):(nMCMC-1)]), file = paste("nVarsSeparate", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(mean(n.vars.tilde[(BI+1):(nMCMC-1)]), file = paste("nVarsTildeSeparate", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    ##########
    ### Shared
    ##########
    # Initialize.
    if(TRUE) {
      beta.tilde = rbind(beta.tilde.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.tilde.true)))
      beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
      mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      mu.diff = rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
      tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
      
      # Baseline and change mdoels share inclusion indicator.
      w = rbind(rep(1, length(beta.true) - 1), matrix(NA, nrow = nMCMC-1, ncol = K))
      
      a = rep(0.1, length(beta.true) - 2)
      gamma = rbind(rep(1/25, K), matrix(NA, nrow = nMCMC - 1, ncol = K))
      sigma2.tilde.b = sigma2.b = rep(NA, J)
      n.vars = rep(NA, nMCMC-1)
    }
    
    # Gibbs sampler.
    for(t in 2:nMCMC) {
      mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
      mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
      mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
      
      for(j in 1:J) {
        sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
        sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
      }
      
      sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
      sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
      
      tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X%*%beta.tilde[t-1,])^2))
      tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
      
      # Use the multi-variate normal representation of the spike-and-slab prior.
      a = ifelse(w[t-1,], 1/sqrt(gamma[t-1,]), 0.1)
      D = diag(c(100, a))
      DD.inv = solve(D%*%D)     
      
      V.tilde = solve(gram/tau2.tilde[t] + DD.inv)
      V = solve(gram/tau2[t] + DD.inv)
      
      beta.tilde[t,] = rmvnorm(1, (1/tau2.tilde[t])*V.tilde%*%t(X)%*%mu[t,], V.tilde)
      beta[t,] = rmvnorm(1, (1/tau2[t])*V%*%t(X)%*%mu.diff[t,], V)
      
      # Bayes Factors: P(beta | slab) / P(beta | spike)
      BF = NULL
      for(k in 1:K) {
        BF[k] = exp(dmvnorm(c(beta.tilde[t,k+1],beta[t,k+1]), sigma = diag(0.1^2, 2), log = TRUE) - dmvnorm(c(beta.tilde[t,k+1],beta[t,k+1]), sigma = diag(1/sqrt(gamma[t-1,k]),2), log = TRUE))
      }
      
      # P(w = 1 | beta) = P(beta | w = 1)P(w = 1) / [P(beta | w = 1)P(w = 1) + P(beta | w = 0)P(w = 0)] = 1 / (1 + BF)
      selp = 1 / (1 + BF)
      
      for(k in 1:K) {
        w[t,k] = sample(0:1, 1, prob = c(1-selp[k], selp[k]))
      }
      
      # Number of variables selected in draw t.
      n.vars[t-1] = sum(w[t,])
      
      gamma[t,] = rgamma(K, shape = 2.5 + 0.5*w[t,], rate = 2.5*5^2 + 0.5*w[t,]*(beta.tilde[t,-1]^2 + beta[t,-1]^2))
    }
    
    ############################################
    ## Write iteration "n" results to text files
    ############################################
    
    # Compute posterior mean of Delta
    betaShared = colMeans(beta[postBI,])[1]
    
    # Compute inclusion probabilities for each covariate
    inclusionShared = colMeans(w[postBI,])
    
    # Write Bias for Delta and inclusion
    write.table(t(betaShared - beta.true[1]), file = paste("biasBetaShared", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(t(inclusionShared), file = paste("inclusionShared", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write MSE
    write.table(t((betaShared - beta.true[1])^2), file = paste("MSE_Shared", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write Coverage for Delta
    CI = quantile(beta[postBI,2], probs = c(.025, .975))
    write.table(ifelse(CI[1] < beta.true[1] & beta.true[1] < CI[2], 1, 0), file = paste("Coverage_Shared", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write number of covariates selected
    write.table(mean(n.vars[(BI+1):(nMCMC-1)]), file = paste("nVarsShared", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    #####################
    ### Sufficient method
    #####################
    # Initialize.
    if(TRUE) {
      beta.tilde = rbind(beta.tilde.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.tilde.true)))
      beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
      mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      mu.diff = rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.alpha = c(1, rep(NA, nMCMC-1))
      tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
      
      w.tilde = rbind(beta.tilde.true[-1], matrix(NA, nrow = nMCMC-1, ncol = K))
      w = rbind(beta.true[-1], matrix(NA, nrow = nMCMC-1, ncol = K))
      w.c = rbind(alpha.true[-1], matrix(NA, nrow = nMCMC-1, ncol = K))
      
      a.tilde = a = rep(0.1, length(beta.true) - 2)
      gamma = gamma.tilde = rbind(rep(1/25, K), matrix(NA, nrow = nMCMC - 1, ncol = K))
      gamma.c = rbind(rep(1/25, K), matrix(NA, nrow = nMCMC - 1, ncol = K))
      sigma2.tilde.b = sigma2.b = rep(NA, J)

      alpha = rbind(alpha.true, matrix(NA, nrow = nMCMC - 1, ncol = length(alpha.true)))
      n.vars = rep(NA, nMCMC-1)
      n.vars.tilde = rep(NA, nMCMC-1)
    }
    
    # Gibbs sampler.
    for(t in 2:nMCMC) {
      mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
      mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
      mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
      
      for(j in 1:J) {
        sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
        sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
      }
      
      sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
      sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
      
      tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X%*%beta.tilde[t-1,])^2))
      tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
      
      # Use multi-variate normal representation of the spike-and-slab prior.
      ## Change model.
      a = ifelse(w[t-1,], 1/sqrt(gamma[t-1,]), 0.1)
      D = diag(c(100, a))
      DD.inv = solve(D%*%D)
      V = solve(gram/tau2[t] + DD.inv)
      beta[t,] = rmvnorm(1, (1/tau2[t])*V%*%t(X)%*%mu.diff[t,], V)
      
      # Baseline model.
      a.tilde = ifelse(w.tilde[t-1,], 1/sqrt(gamma.tilde[t-1,]), 0.1)
      D.tilde = diag(c(100, a.tilde))
      DD.tilde.inv = solve(D.tilde%*%D.tilde)
      V.tilde = solve(gram/tau2.tilde[t] + DD.tilde.inv)
      beta.tilde[t,] = rmvnorm(1, (1/tau2.tilde[t])*V.tilde%*%t(X)%*%mu[t,], V.tilde)
      
      # Exposure model.
      a.c = ifelse(w.c[t-1,], 1/sqrt(gamma.c[t-1,]), 0.1)
      D.c = diag(a.c)
      DcRcDc.inv = solve(D.c%*%D.c)     
      V.alpha = solve(gram0/sigma2.alpha[t-1] + DcRcDc.inv)
      alpha[t,] = rmvnorm(1, (1/sigma2.alpha[t-1])*V.alpha%*%t(X0)%*%Trt, V.alpha)
      sigma2.alpha[t] = 1 / rgamma(1, J/2, 0.5*sum((Trt - X0%*%alpha[t,])^2))
      
      # Bayes Factors: P(beta | slab) / P(beta | spike)
      # P(w = 1 | beta) = P(beta | w = 1)P(w = 1) / [P(beta | w = 1)P(w = 1) + P(beta | w = 0)P(w = 0)] = 1 / (1 + BF)
      BF.c = exp(dnorm(alpha[t,-1], sd = 0.1, log = TRUE) - dnorm(alpha[t,-1], sd = sqrt(1/gamma.c[t-1,]), log = TRUE))
      p.c = 1 / (1 + BF.c)
      BF.tilde = exp(dnorm(beta.tilde[t,-1], sd = 0.1, log = TRUE) - dnorm(beta.tilde[t,-1], sd = sqrt(1/gamma.tilde[t-1,]), log = TRUE))
      p.tilde = 1 / (1 + BF.tilde)
      BF = exp(dnorm(beta[t,-1], sd = 0.1, log = TRUE) - dnorm(beta[t,-1], sd = 1/sqrt(gamma[t-1,]), log = TRUE))
      p = 1 / (1 + BF)
      
      ## To avoid model feedback, draw w.c without conditioning on the baseline or change models,
      ## and draw w without conditioning on the baseline model.
      for(k in 1:K) {
        w.c[t,k] = sample(0:1, 1, prob = c(1-p.c[k], p.c[k]))
        w[t,k] = w.c[t,k]*sample(0:1, 1, prob = c(1-p[k], p[k]))
        w.tilde[t,k] = w[t,k]*sample(0:1, 1, prob = c(1-p.tilde[k], p.tilde[k]))
      }
    
      # Number of covariates included in each model.
      n.vars.tilde[t-1] = sum(w.tilde[t,])
      n.vars[t-1] = sum(w[t,])
      
      gamma.c[t,] = rgamma(K, shape = 2.5 + 0.5*w.c[t,], rate = 2.5*5^2 + 0.5*w.c[t,]*alpha[t,-1]^2)
      gamma.tilde[t,] = rgamma(K, shape = 2.5 + 0.5*w.tilde[t,], rate = 2.5*5^2 + 0.5*w.tilde[t,]*beta.tilde[t,-1]^2)
      gamma[t,] = rgamma(K, shape = 2.5 + 0.5*w[t,], rate = 2.5*5^2 + 0.5*w[t,]*beta[t,-1]^2)
      
    } 
    
    # Compute posterior mean of Delta
    betaSufficient = colMeans(beta[postBI,])[1]
    
    # Compute inclusion probabilities for each covariate
    inclusionSufficient = colMeans(w[postBI,])
    
    # Write Bias and inclusion
    write.table(t(betaSufficient - beta.true[1]), file = paste("biasBetaSufficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(t(inclusionSufficient), file = paste("inclusionSufficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write MSE
    write.table(t((betaSufficient - beta.true[1])^2), file = paste("MSE_Sufficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write Coverage for Delta
    CI = quantile(beta[postBI,2], probs = c(.025, .975))
    write.table(ifelse(CI[1] < beta.true[1] & beta.true[1] < CI[2], 1, 0), file = paste("Coverage_Sufficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write number of covariates selected
    write.table(mean(n.vars[(BI+1):(nMCMC-1)]), file = paste("nVarsSufficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(mean(n.vars.tilde[(BI+1):(nMCMC-1)]), file = paste("nVarsTildeSufficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    ###################
    ## Efficient method
    ###################
    if(TRUE) {
      beta.tilde = rbind(beta.tilde.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.tilde.true)))
      beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
      mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      mu.diff = rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.alpha = c(1, rep(NA, nMCMC-1))
      tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
      w.tilde = rbind(beta.tilde.true[-1], matrix(NA, nrow = nMCMC-1, ncol = K))
      w = rbind(beta.true[-1], matrix(NA, nrow = nMCMC-1, ncol = K))
      a.tilde = a = rep(0.1, length(beta.true) - 2)
      gamma = gamma.tilde = rbind(rep(1/25, K), matrix(NA, nrow = nMCMC - 1, ncol = K))
      sigma2.tilde.b = sigma2.b = rep(NA, J)
      n.vars = rep(NA, nMCMC-1)
      n.vars.tilde = rep(NA, nMCMC-1)
    }
    
    # Gibbs sampler.
    for(t in 2:nMCMC) {
      mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
      mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
      mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
      
      for(j in 1:J) {
        sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
        sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
      }
      
      sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
      sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
      
      tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X%*%beta.tilde[t-1,])^2))
      tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
      
      # Use multi-variate normal representation of the spike-and-slab prior.
      a = ifelse(w[t-1,], 1/sqrt(gamma[t-1,]), 0.1)
      D = diag(c(100, a))
      DD.inv = solve(D%*%D)
      V = solve(gram/tau2[t] + DD.inv)
      beta[t,] = rmvnorm(1, (1/tau2[t])*V%*%t(X)%*%mu.diff[t,], V)
      
      a.tilde = ifelse(w.tilde[t-1,], 1/sqrt(gamma.tilde[t-1,]), 0.1)
      D.tilde = diag(c(100, a.tilde))
      DD.tilde.inv = solve(D.tilde%*%D.tilde)
      V.tilde = solve(gram/tau2.tilde[t] + DD.tilde.inv)
      beta.tilde[t,] = rmvnorm(1, (1/tau2.tilde[t])*V.tilde%*%t(X)%*%mu[t,], V.tilde)
      
      # Bayes Factors: P(beta | slab) / P(beta | spike)
      BF.tilde = exp(dnorm(beta.tilde[t,-1], sd = 0.1, log = TRUE) - dnorm(beta.tilde[t,-1], sd = sqrt(1/gamma.tilde[t-1,]), log = TRUE))
      # P(w = 1 | beta) = P(beta | w = 1)P(w = 1) / [P(beta | w = 1)P(w = 1) + P(beta | w = 0)P(w = 0)] = 1 / (1 + BF)
      p.tilde = 1 / (1 + BF.tilde)
      BF = exp(dnorm(beta[t,-1], sd = 0.1, log = TRUE) - dnorm(beta[t,-1], sd = 1/sqrt(gamma[t-1,]), log = TRUE))
      p = 1 / (1 + BF)
      
      for(k in 1:K) {
        w[t,k] = sample(0:1, 1, prob = c(1-p[k], p[k]))
        w.tilde[t,k] = w[t,k]*sample(0:1, 1, prob = c(1-p.tilde[k], p.tilde[k]))
      }
      
      # Number of covariates included in each model
      n.vars.tilde[t-1] = sum(w.tilde[t,])
      n.vars[t-1] = sum(w[t,])
      
      gamma.tilde[t,] = rgamma(K, shape = 2.5 + 0.5*w.tilde[t,], rate = 2.5*5^2 + 0.5*w.tilde[t,]*beta.tilde[t,-1]^2)
      gamma[t,] = rgamma(K, shape = 2.5 + 0.5*w[t,], rate = 2.5*5^2 + 0.5*w[t,]*beta[t,-1]^2)
      
    } 
    
    ############################################
    ## Write iteration "n" results to text files
    ############################################
    
    # Compute posterior mean of Delta
    betaEfficient = colMeans(beta[postBI,])[1]
    
    # Compute inclusion probabilities for each covariate
    inclusionEfficient = colMeans(w[postBI,])
    
    # Write Bias and inclusion
    write.table(t(betaEfficient - beta.true[1]), file = paste("biasBetaEfficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(t(inclusionEfficient), file = paste("inclusionEfficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write MSE
    write.table(t((betaEfficient - beta.true[1])^2), file = paste("MSE_Efficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write Coverage for Delta
    CI = quantile(beta[postBI,2], probs = c(.025, .975))
    write.table(ifelse(CI[1] < beta.true[1] & beta.true[1] < CI[2], 1, 0), file = paste("Coverage_Efficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
    # Write number of covariates selected
    write.table(mean(n.vars[(BI+1):(nMCMC-1)]), file = paste("nVarsEfficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(mean(n.vars.tilde[(BI+1):(nMCMC-1)]), file = paste("nVarsTildeEfficient", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    
  }
}

alpha.true = c(rep(1, 4), rep(0, 4))
beta.tilde.true = c(0, 1, 1, 0, 0, 1, 1, 0, 0)
beta.true = c(1, 1, 0, 1, 0, 1, 0, 1, 0)

setwd("/home/normington/Paper 2/PosteriorSimResults")
string0 = "Test_Run"

# Use multiple cores for faster computation.
n.cores = 14

cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=1:n.cores) %dopar% {
  run_sim(NSIM = ceiling(5000/n.cores), J=50, beta.tilde.true = beta.tilde.true, beta.true = beta.true, alpha.true = alpha.true, string = string0, seed = 2019)
}
stopCluster(cl) 

cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=1:n.cores) %dopar% {
  run_sim(NSIM = ceiling(5000/n.cores), J=100, beta.tilde.true = beta.tilde.true, beta.true = beta.true, alpha.true = alpha.true, string = string0, seed = 2019)
}
stopCluster(cl) 

for(J in c(50, 100)) {
  string = paste0(string0, J)
  
  # Bias, beta 
  biasBetaSeparate = colMeans(read.table(paste0("biasBetaSeparate", string, ".txt")))
  biasBetaShared = colMeans(read.table(paste0("biasBetaShared", string, ".txt")))
  biasBetaSufficient = colMeans(read.table(paste0("biasBetaSufficient", string, ".txt")))
  biasBetaEfficient = colMeans(read.table(paste0("biasBetaEfficient", string, ".txt")))
  
  # MSE, beta
  MSEBetaSeparate = colMeans(read.table(paste0("MSE_Separate", string, ".txt")))
  MSEBetaShared = colMeans(read.table(paste0("MSE_Shared", string, ".txt")))
  MSEBetaSufficient = colMeans(read.table(paste0("MSE_Sufficient", string, ".txt")))
  MSEBetaEfficient = colMeans(read.table(paste0("MSE_Efficient", string, ".txt")))
  
  # Inclusion probabilities, mu_diff
  inclusionSeparate = colMeans(read.table(paste0("inclusionSeparate", string, ".txt")))
  inclusionShared = colMeans(read.table(paste0("inclusionShared", string, ".txt")))
  inclusionSufficient = colMeans(read.table(paste0("inclusionSufficient", string, ".txt")))
  inclusionEfficient = colMeans(read.table(paste0("inclusionEfficient", string, ".txt")))
  
  # Coverage probabilities, beta
  CoverageSeparate = mean(read.table(paste0("Coverage_Separate", string, ".txt"))[,1])
  CoverageShared = mean(read.table(paste0("Coverage_Shared", string, ".txt"))[,1])
  CoverageSufficient = mean(read.table(paste0("Coverage_Sufficient", string, ".txt"))[,1])
  CoverageEfficient = mean(read.table(paste0("Coverage_Efficient", string, ".txt"))[,1])
  
  # Number of predictors, mudiff
  n.varsSeparate = round(mean(read.table(paste0("nVarsSeparate", string, ".txt"))[,1]), 2)
  n.varsShared = round(mean(read.table(paste0("nVarsShared", string, ".txt"))[,1]), 2)
  n.varsSufficient = round(mean(read.table(paste0("nVarsSufficient", string, ".txt"))[,1]), 2)
  n.varsEfficient = round(mean(read.table(paste0("nVarsEfficient", string, ".txt"))[,1]), 2)
  
  # Number of predictors, mu
  n.VarsTildeSeparate = round(mean(read.table(paste0("nVarsTildeSeparate", string, ".txt"))[,1]), 2)
  n.VarsTildeShared = n.varsShared
  n.VarsTildeSufficient = round(mean(read.table(paste0("nVarsTildeSufficient", string, ".txt"))[,1]), 2)
  n.VarsTildeEfficient = round(mean(read.table(paste0("nVarsTildeEfficient", string, ".txt"))[,1]), 2)
  
  # Write out results
  round(c(biasBetaSeparate, biasBetaShared, biasBetaSufficient, biasBetaEfficient), 4)
  round(c(MSEBetaSeparate, MSEBetaShared, MSEBetaSufficient, MSEBetaEfficient), 4)
  round(rbind(inclusionSeparate, inclusionShared, inclusionSufficient, inclusionEfficient), 4)
  round(c(CoverageSeparate, CoverageShared, CoverageSufficient, CoverageEfficient), 4)
  
  print(c(n.varsSeparate, n.varsShared, n.varsSufficient, n.varsEfficient))
  print(c(n.VarsTildeSeparate, n.VarsTildeShared, n.VarsTildeSufficient, n.VarsTildeEfficient))
}