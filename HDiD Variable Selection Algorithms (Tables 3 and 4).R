################################################################
##### Code to produce Tables 4 and 5 in "Bayesian variable selection 
##### in hierarchical difference-in-differences models"
##### Submitted to Biometrics, October 2019
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
    if(0) {
      ##############
      ### Full Model
      ##############
      # Initialize vectors and matrices to store MCMC draws.
      if(TRUE) {
        beta.tilde = rbind(beta.tilde.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.tilde.true)))
        beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
        mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
        mu.diff = rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
        sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
        tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
        sigma2.tilde.b = sigma2.b = rep(NA, J)
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
        
        V.tilde = solve(gram/tau2.tilde[t])
        V = solve(gram/tau2[t])
        
        beta.tilde[t,] = rmvnorm(1, (1/tau2.tilde[t])*V.tilde%*%t(X)%*%mu[t,], V.tilde)
        beta[t,] = rmvnorm(1, (1/tau2[t])*V%*%t(X)%*%mu.diff[t,], V)
        
      }
      
      ############################################
      ## Write iteration "n" results to text files
      ############################################
      
      # Compute posterior mean of Delta
      betaFull = colMeans(beta[postBI,])[1]
      
      # Compute inclusion probabilities for each covariate
      inclusionFull = colMeans(w[postBI,])
      
      # Write Bias for Delta 
      write.table(t(betaFull - beta.true[1]), file = paste("biasBetaFull", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(t(inclusionFull), file = paste("inclusionFull", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      
      # Write MSE for Delta
      write.table(t((betaFull - beta.true[1])^2), file = paste("MSE_Full", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      
      # Write Coverage for Delta
      CI = quantile(beta[postBI,2], probs = c(.025, .975))
      write.table(ifelse(CI[1] < beta.true[1] & beta.true[1] < CI[2], 1, 0), file = paste("Coverage_Full", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      
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
    }
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
      BF.c = exp(dnorm(alpha[t,], sd = 0.1, log = TRUE) - dnorm(alpha[t,], sd = sqrt(1/gamma.c[t-1,]), log = TRUE))
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
    
    if(0) {
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
      
      ##############
      ### Null Model
      ##############
      # Initialize vectors and matrices to store MCMC draws.
      if(TRUE) {
        mu0 = c(0, rep(NA, nMCMC-1))
        Delta = c(0, rep(NA, nMCMC-1))
        mu = rbind(rep(0, J), matrix(NA, nrow = nMCMC-1, ncol = J))
        mu.diff = rbind(rep(0, J), matrix(NA, nrow = nMCMC-1, ncol = J))
        sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
        tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
        sigma2.tilde.b = sigma2.b = rep(NA, J)
        gramT.inv = solve(t(Trt)%*%Trt)
      }
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*mu0[t-1]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*Trt*Delta[t-1]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - mu0[t-1])^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - Trt*Delta[t-1])^2))
        
        mu0[t] = rnorm(1, mean(mu[t,]), sqrt(tau2.tilde[t] / J))
        Delta[t] = rnorm(1, as.numeric(gramT.inv%*%t(Trt)%*%mu.diff[t,]), sqrt(as.numeric(tau2[t]*gramT.inv)))
        
      }
      
      ############################################
      ## Write iteration "n" results to text files
      ############################################
      
      # Compute posterior mean of Delta
      betaNull = mean(Delta[postBI])
      
      # Compute inclusion probabilities for each covariate
      inclusionNull = colMeans(w[postBI,])
      
      # Write Bias for Delta 
      write.table(t(betaNull - beta.true[1]), file = paste("biasBetaNull", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      write.table(t(inclusionNull), file = paste("inclusionNull", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      
      # Write MSE for Delta
      write.table(t((betaNull - beta.true[1])^2), file = paste("MSE_Null", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
      
      # Write Coverage for Delta
      CI = quantile(Delta[postBI], probs = c(.025, .975))
      write.table(ifelse(CI[1] < beta.true[1] & beta.true[1] < CI[2], 1, 0), file = paste("Coverage_Null", string, J, ".txt", sep = ""), append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
}

alpha.true = c(rep(1, 4), rep(0, 4))
beta.tilde.true = c(0, 1, 1, 0, 0, 1, 1, 0, 0)
beta.true = c(1, 1, 0, 1, 0, 1, 0, 1, 0)

setwd("/home/normington/Paper 2/PosteriorSimResults")
string0 = "Test_Run"

# Use multiple cores for faster computation.
n.cores = 25

# Random seed for each core.
seeds = c(2019, 1992, 2018, 37, 777, 337, 554, 
          654, 2014, 89, 84, 2010, 66, 321,
          83, 97, 12345, 54321, 1776, 1692, 2008, 
          3023, 5, 65, 23)

cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=1:n.cores) %dopar% {
  run_sim(NSIM = ceiling(2200/n.cores), J=50, beta.tilde.true = beta.tilde.true, beta.true = beta.true, alpha.true = alpha.true, string = string0, 
          seed = seeds[i])
}
stopCluster(cl) 

cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=1:n.cores) %dopar% {
  run_sim(NSIM = ceiling(2200/n.cores), J=100, beta.tilde.true = beta.tilde.true, beta.true = beta.true, alpha.true = alpha.true, string = string0, 
          seed = seeds[i])
}
stopCluster(cl) 

z.star = qnorm(0.975)
for(J in c(50, 100)) {
  string = paste0(string0, J)
  
  print(J)
  
  biasFull = read.table(paste0("biasBetaFull", string, ".txt"))[,1]
  biasSeparate = read.table(paste0("biasBetaSeparate", string, ".txt"))[,1]
  biasShared = read.table(paste0("biasBetaShared", string, ".txt"))[,1]
  biasSufficient = read.table(paste0("biasBetaSufficient", string, ".txt"))[,1]
  biasEfficient = read.table(paste0("biasBetaEfficient", string, ".txt"))[,1]
  biasNull = read.table(paste0("biasBetaNull", string, ".txt"))[,1]
  
  # Bias, beta 
  biasBetaFull = round(mean(biasFull), 3)
  biasBetaSeparate = round(mean(biasSeparate), 3)
  biasBetaShared = round(mean(biasShared), 3)
  biasBetaSufficient = round(mean(biasSufficient), 3)
  biasBetaEfficient = round(mean(biasEfficient), 3)
  biasBetaNull = round(mean(biasNull), 3)
  
  # MoE Bias
  biasMoEFull = round(z.star*sd(biasFull) / sqrt(length(biasFull)), 3)
  biasMoESeparate = round(z.star*sd(biasSeparate) / sqrt(length(biasSeparate)), 3)
  biasMoEShared = round(z.star*sd(biasShared) / sqrt(length(biasShared)), 3)
  biasMoESufficient = round(z.star*sd(biasSufficient) / sqrt(length(biasSufficient)), 3)
  biasMoEEfficient = round(z.star*sd(biasEfficient) / sqrt(length(biasEfficient)), 3)
  biasMoENull = round(z.star*sd(biasNull) / sqrt(length(biasNull)), 3)
  
  MSEFull = read.table(paste0("MSE_Full", string, ".txt"))[,1]
  MSESeparate = read.table(paste0("MSE_Separate", string, ".txt"))[,1]
  MSEShared = read.table(paste0("MSE_Shared", string, ".txt"))[,1]
  MSESufficient = read.table(paste0("MSE_Sufficient", string, ".txt"))[,1]
  MSEEfficient = read.table(paste0("MSE_Efficient", string, ".txt"))[,1]
  MSENull = read.table(paste0("MSE_Null", string, ".txt"))[,1]
  
  # MSE, beta
  MSEBetaFull = round(mean(MSEFull), 3)
  MSEBetaSeparate = round(mean(MSESeparate), 3)
  MSEBetaShared = round(mean(MSEShared), 3)
  MSEBetaSufficient = round(mean(MSESufficient), 3)
  MSEBetaEfficient = round(mean(MSEEfficient), 3)
  MSEBetaNull = round(mean(MSENull), 3)
  
  # MoE MSE
  MSEMoEFull = round(z.star*sd(MSEFull) / sqrt(length(MSEFull)), 3)
  MSEMoESeparate = round(z.star*sd(MSESeparate) / sqrt(length(MSESeparate)), 3)
  MSEMoEShared = round(z.star*sd(MSEShared) / sqrt(length(MSEShared)), 3)
  MSEMoESufficient = round(z.star*sd(MSESufficient) / sqrt(length(MSESufficient)), 3)
  MSEMoEEfficient = round(z.star*sd(MSEEfficient) / sqrt(length(MSEEfficient)), 3)
  MSEMoENull = round(z.star*sd(MSENull) / sqrt(length(MSENull)), 3)
  
  # Coverage probabilities, beta
  CoverageFull = mean(read.table(paste0("Coverage_Full", string, ".txt"))[,1])
  CoverageSeparate = mean(read.table(paste0("Coverage_Separate", string, ".txt"))[,1])
  CoverageShared = mean(read.table(paste0("Coverage_Shared", string, ".txt"))[,1])
  CoverageSufficient = mean(read.table(paste0("Coverage_Sufficient", string, ".txt"))[,1])
  CoverageEfficient = mean(read.table(paste0("Coverage_Efficient", string, ".txt"))[,1])
  CoverageNull = mean(read.table(paste0("Coverage_Null", string, ".txt"))[,1])
  
  # Number of predictors, mudiff
  varsSeparate = read.table(paste0("nVarsSeparate", string, ".txt"))[,1]
  varsShared = read.table(paste0("nVarsShared", string, ".txt"))[,1]
  varsSufficient = read.table(paste0("nVarsSufficient", string, ".txt"))[,1]
  varsEfficient = read.table(paste0("nVarsEfficient", string, ".txt"))[,1]
  
  n.varsSeparate = round(mean(varsSeparate), 2)
  n.varsShared = round(mean(varsShared), 2)
  n.varsSufficient = round(mean(varsSufficient), 2)
  n.varsEfficient = round(mean(varsEfficient), 2)
  
  # MoE, number of predictors mudiff
  n.varsSeparateMoE = round(z.star * sd(varsSeparate) / sqrt(length(varsSeparate)), 2)
  n.varsSharedMoE = round(z.star * sd(varsShared) / sqrt(length(varsShared)), 2)
  n.varsSufficientMoE = round(z.star * sd(varsSufficient) / sqrt(length(varsSufficient)), 2)
  n.varsEfficientMoE = round(z.star * sd(varsEfficient) / sqrt(length(varsEfficient)), 2)
  
  # Number of predictors, mu
  varsTildeSeparate = read.table(paste0("nVarsTildeSeparate", string, ".txt"))[,1]
  varsTildeShared = varsShared
  varsTildeSufficient = read.table(paste0("nVarsTildeSufficient", string, ".txt"))[,1]
  varsTildeEfficient = read.table(paste0("nVarsTildeEfficient", string, ".txt"))[,1]
  
  n.varsTildeSeparate = round(mean(varsTildeSeparate), 2)
  n.varsTildeShared = round(mean(varsTildeShared), 2)
  n.varsTildeSufficient = round(mean(varsTildeSufficient), 2)
  n.varsTildeEfficient = round(mean(varsTildeEfficient), 2)
  
  # MoE, number of predictors mudiff
  n.varsTildeSeparateMoE = round(z.star * sd(varsTildeSeparate) / sqrt(length(varsTildeSeparate)), 2)
  n.varsTildeSharedMoE = round(z.star * sd(varsTildeShared) / sqrt(length(varsTildeShared)), 2)
  n.varsTildeSufficientMoE = round(z.star * sd(varsTildeSufficient) / sqrt(length(varsTildeSufficient)), 2)
  n.varsTildeEfficientMoE = round(z.star * sd(varsTildeEfficient) / sqrt(length(varsTildeEfficient)), 2)
  
  # Inclusion probabilities, mu_diff
  inclusionSeparate = round(colMeans(read.table(paste0("inclusionSeparate", string, ".txt"))), 3)
  inclusionShared = round(colMeans(read.table(paste0("inclusionShared", string, ".txt"))), 3)
  inclusionSufficient = round(colMeans(read.table(paste0("inclusionSufficient", string, ".txt"))), 3)
  inclusionEfficient = round(colMeans(read.table(paste0("inclusionEfficient", string, ".txt"))), 3)
  
  ## Table 3
  print(round(c(biasBetaFull, biasBetaSeparate, biasBetaShared, biasBetaSufficient, biasBetaEfficient, biasBetaNull), 4))
  print(round(c(biasMoEFull, biasMoESeparate, biasMoEShared, biasMoESufficient, biasMoEEfficient, biasMoENull), 4))
  
  print(round(c(MSEBetaFull, MSEBetaSeparate, MSEBetaShared, MSEBetaSufficient, MSEBetaEfficient, MSEBetaNull), 4))
  print(round(c(MSEMoEFull, MSEMoESeparate, MSEMoEShared, MSEMoESufficient, MSEMoEEfficient, MSEMoENull), 4))
  
  print(round(c(CoverageFull, CoverageSeparate, CoverageShared, CoverageSufficient, CoverageEfficient, CoverageNull), 4))
  
  
  print(c(n.varsSeparate, n.varsShared, n.varsSufficient, n.varsEfficient))
  print(c(n.varsSeparateMoE, n.varsSharedMoE, n.varsSufficientMoE, n.varsEfficientMoE))
  
  print(c(n.varsTildeSeparate, n.varsTildeShared, n.varsTildeSufficient, n.varsTildeEfficient))
  print(c(n.varsTildeSeparateMoE, n.varsTildeSharedMoE, n.varsTildeSufficientMoE, n.varsTildeEfficientMoE))
  
  
  ## Table 4
  print(round(rbind(inclusionSeparate, inclusionShared, inclusionSufficient, inclusionEfficient), 4))
}
