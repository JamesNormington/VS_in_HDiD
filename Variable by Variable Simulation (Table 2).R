################################################################
##### Code to produce Table 2 in "Bayesian variable selection 
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

run_sim = function(NSIM=5000, nMCMC=10000, BI=1000, J=50, beta.tilde.true, beta.true, alpha.true, string, case, scen, label, seed) {
  setwd("/home/normington/Paper 2/PosteriorSimResults")
  require(mvtnorm)
  
  # Number of clinics and predictors.
  n.pre = n.post = rep(10, J)
  K = length(beta.true) - 1
  
  # Initialize vectors and matrices to store MCMC draws.
  beta =  rbind(beta.true, matrix(0, nrow = nMCMC-1, ncol = length(beta.true)))
  mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
  mu.diff= rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
  sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
  gamma.tilde = gamma = c(5, rep(NA, nMCMC - 1))
  tau2.tilde = tau2 = gamma2.tilde = gamma2 = c(1, rep(NA, nMCMC-1))
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
  
  for(nsim in 1:NSIM) {
    
    ##################
    ### Generate data.
    ##################
    X.case = rnorm(J)
    Trt =  rnorm(J, alpha.true[case]*X.case)
    mu = rnorm(J, beta.tilde.true[case+2]*X.case)
    mu.diff = rnorm(J, beta.true[2]*Trt + beta.true[case+2]*X.case)
    for(j in 1:J) {
      Y.pre[[j]] = rnorm(n.pre[j], mu[j])
      Y.post[[j]] = rnorm(n.post[j], mu[j] + mu.diff[j])
    }
    
    Y.pre.means = unlist(lapply(Y.pre, mean))
    Y.post.means = unlist(lapply(Y.post, mean))
    
    # Initialize for Gibbs sampling in each scenario.
    if(TRUE) {
      mu = rbind(rep(beta.tilde.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      mu.diff = rbind(rep(beta.true[1], J), matrix(NA, nrow = nMCMC-1, ncol = J))
      sigma2.tilde = sigma2 = rbind(rep(1, J), matrix(NA, nrow = nMCMC-1, ncol = J))
      tau2.tilde = tau2 = c(1, rep(NA, nMCMC - 1))
      sigma2.tilde.b = sigma2.b = rep(NA, J)
    }
    
    ## Scenario 1: no adjustment for mu or mu^diff
    if(scen == 1) {
      beta.tilde = c(0, rep(NA, nMCMC-1))
      beta = rbind(c(0, 0), matrix(NA, nMCMC-1, 2))
      X.tilde = rep(1, J)
      X = cbind(1, Trt)
      gram = t(X)%*%X
      gram.inv = solve(gram)
      gram.tilde.inv = matrix(1/J)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*beta.tilde[t-1]*X.tilde) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - beta.tilde[t-1]*X.tilde)^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
        beta.tilde[t] = rnorm(1, as.numeric(gram.tilde.inv%*%t(X.tilde)%*%mu[t,]), as.numeric(tau2.tilde[t]*gram.tilde.inv))
      }
    }
    
    ## Scenario 2: no adjustment for mu^diff
    if(scen == 2) {
      beta.tilde = rbind(c(0, 0), matrix(NA, nMCMC-1, 2))
      beta = rbind(c(0, 0), matrix(NA, nMCMC-1, 2))
      X.tilde = cbind(1, X.case)
      X = cbind(1, Trt)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X.tilde%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X.tilde%*%beta.tilde[t-1,])^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
        beta.tilde[t,] = rmvnorm(1, gram.tilde.inv%*%t(X.tilde)%*%mu[t,], tau2.tilde[t]*gram.tilde.inv)
      }
    }
    
    ## Scenario 3: no adjustment for mu
    if(scen == 3) {
      beta.tilde = c(0, rep(NA, nMCMC-1))
      beta = rbind(c(0, 0, 0), matrix(NA, nMCMC-1, 3))
      X.tilde = rep(1, J)
      X = cbind(1, Trt, X.case)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*beta.tilde[t-1]*X.tilde) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - beta.tilde[t-1]*X.tilde)^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta.tilde[t] = rnorm(1, as.numeric(gram.tilde.inv%*%t(X.tilde)%*%mu[t,]), as.numeric(tau2.tilde[t]*gram.tilde.inv))
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
      }
    }
    
    ## Scenario 4: adjustment for both
    if(scen == 4) {
      beta.tilde = rbind(c(0, 0), matrix(NA, nMCMC-1, 2))
      beta = rbind(c(0, 0, 0), matrix(NA, nMCMC-1, 3))
      X.tilde = cbind(1, X.case)
      X = cbind(1, Trt, X.case)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X.tilde%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X.tilde%*%beta.tilde[t-1,])^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
        beta.tilde[t,] = rmvnorm(1, gram.tilde.inv%*%t(X.tilde)%*%mu[t,], tau2.tilde[t]*gram.tilde.inv)
      }
    }
    
    ## Scenario 5: no adjustment for mu or mu^diff, adjust mu for T
    if(scen == 5) {
      beta.tilde = c(0, rep(NA, nMCMC-1))
      beta = rbind(c(0, 0), matrix(NA, nMCMC-1, 2))
      X.tilde = cbind(1, Trt)
      X = cbind(1, Trt)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*beta.tilde[t-1]*X.tilde) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - beta.tilde[t-1]*X.tilde)^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
        beta.tilde[t] = rnorm(1, as.numeric(gram.tilde.inv%*%t(X.tilde)%*%mu[t,]), as.numeric(tau2.tilde[t]*gram.tilde.inv))
      }
    }
    
    ## Scenario 6: no adjustment for mu^diff, adjust mu for T
    if(scen == 6) {
      beta.tilde = rbind(c(0, 0, 0), matrix(NA, nMCMC-1, 3))
      beta = rbind(c(0, 0), matrix(NA, nMCMC-1, 2))
      X.tilde = cbind(1, Trt, X.case)
      X = cbind(1, Trt)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X.tilde%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X.tilde%*%beta.tilde[t-1,])^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
        beta.tilde[t,] = rmvnorm(1, gram.tilde.inv%*%t(X.tilde)%*%mu[t,], tau2.tilde[t]*gram.tilde.inv)
      }
    }
    
    ## Scenario 7: no adjustment for mu, adjust mu for T
    if(scen == 7) {
      beta.tilde = c(0, rep(NA, nMCMC-1))
      beta = rbind(c(0, 0, 0), matrix(NA, nMCMC-1, 3))
      X.tilde = cbind(1, Trt)
      X = cbind(1, Trt, X.case)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*beta.tilde[t-1]*X.tilde) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - beta.tilde[t-1]*X.tilde)^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))

        beta.tilde[t] = rnorm(1, as.numeric(gram.tilde.inv%*%t(X.tilde)%*%mu[t,]), as.numeric(tau2.tilde[t]*gram.tilde.inv))
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
      }
    }
    
    ## Scenario 8: adjustment for both, adjust mu for T
    if(scen == 8) {
      beta.tilde = rbind(c(0, 0, 0), matrix(NA, nMCMC-1, 3))
      beta = rbind(c(0, 0, 0), matrix(NA, nMCMC-1, 3))
      X.tilde = cbind(1, Trt, X.case)
      X = cbind(1, Trt, X.case)
      gram.inv = solve(t(X)%*%X)
      gram.tilde.inv = solve(t(X.tilde)%*%X.tilde)
      
      # Gibbs sampler.
      for(t in 2:nMCMC) {
        mu.mean = (n.pre*sigma2[t-1,]*tau2.tilde[t-1]*Y.pre.means + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1]*(Y.post.means - mu.diff[t-1,]) + sigma2.tilde[t-1,]*sigma2[t-1,]*X.tilde%*%beta.tilde[t-1,]) / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])
        mu[t,] = rnorm(J, mu.mean, sqrt(sigma2.tilde[t-1,]*sigma2[t-1,]*tau2.tilde[t-1] / (n.pre*sigma2[t-1,]*tau2.tilde[t-1] + n.post*sigma2.tilde[t-1,]*tau2.tilde[t-1] + sigma2.tilde[t-1,]*sigma2[t-1,])))
        mu.diff[t,] = rnorm(J, (n.post*tau2[t-1]*(Y.post.means - mu[t,]) + sigma2[t-1,]*X%*%beta[t-1,]) / (n.post*tau2[t-1] + sigma2[t-1,]), sqrt(sigma2[t-1,]*tau2[t-1] / (n.post*tau2[t-1] + sigma2[t-1,])))
        
        for(j in 1:J) {
          sigma2.tilde.b[j] = sum((Y.pre[[j]] - mu[t,j])^2)
          sigma2.b[j] = sum((Y.post[[j]] - mu[t,j] - mu.diff[t,j])^2)
        }
        
        sigma2.tilde[t,] = 1 / rgamma(J, n.pre/2, 0.5*sigma2.tilde.b)
        sigma2[t,] = 1 / rgamma(J, n.post/2, 0.5*sigma2.b)
        
        tau2.tilde[t] = 1 / rgamma(1, J/2, 0.5*sum((mu[t,] - X.tilde%*%beta.tilde[t-1,])^2))
        tau2[t] = 1 / rgamma(1, J/2, 0.5*sum((mu.diff[t,] - X%*%beta[t-1,])^2))
        
        beta[t,] = rmvnorm(1, gram.inv%*%t(X)%*%mu.diff[t,], tau2[t]*gram.inv)
        beta.tilde[t,] = rmvnorm(1, gram.tilde.inv%*%t(X.tilde)%*%mu[t,], tau2.tilde[t]*gram.tilde.inv)
      }
    }
    
    # Compute posterior means of beta.
    betaHat = colMeans(beta[(BI+1):nMCMC,])
    
    label2 = paste0(label, "Scen", scen)
    
    # Write Bias, MSE, and Coverage
    write.table(t(betaHat[2] - beta.true[2]), file = paste0("bias", label2, string, ".txt"), append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(t((betaHat[2] - beta.true[2])^2), file = paste0("MSE", label2, string, ".txt"), append = TRUE, col.names = FALSE, row.names = FALSE)
    CI = quantile(beta[(BI+1):nMCMC,2], probs = c(.025, .975))
    write.table(ifelse(CI[1] < beta.true[2] & beta.true[2] < CI[2], 1, 0), file = paste0("Coverage", label2, string, ".txt"), append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}

# True values for parameters; each include an intercept = 0.
alpha.true = c(0, rep(1, 4), rep(0, 4))
beta.tilde.true = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0)
beta.true = c(0, 1, 1, 0, 1, 0, 1, 0, 1, 0)

# Combine cases in a matrix.
cases_matrix = cbind(rep(1:8, each = 4), rep(5:8, 8))
label_vec = rep(paste0("X", 1:8), each = 4)

string0 = "Test_Sep9_v2"
n.cores = 16
seeds = c(2019, 1992, 2018, 37, 777, 337, 554, 
          654, 2014, 89, 84, 87, 63, 222, 
          23, 32)

cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=1:n.cores) %dopar% {
  run_sim(NSIM = 10, nMCMC = 5000, BI = 500, J=50, beta.tilde.true = beta.tilde.true, 
               beta.true = beta.true, alpha.true = alpha.true, string = string0,
               case = cases_matrix[i,1], scen = cases_matrix[i,2], label = label_vec[i], 
               seed = seeds[i])
}
stopCluster(cl) 

cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=(n.cores+1):(2*n.cores)) %dopar% {
  run_sim(NSIM = 5000, nMCMC = 10000, BI = 1000, J=50, beta.tilde.true = beta.tilde.true, 
          beta.true = beta.true, alpha.true = alpha.true, string = string0,
          case = cases_matrix[i,1], scen = cases_matrix[i,2], label = label_vec[i],
          seed = seeds[i]-n.cores)
}
stopCluster(cl) 

n.cores = 10
cl = makeCluster(n.cores)
registerDoParallel(cl)
foreach(i=23:32) %dopar% {
  run_sim(NSIM = 5000, nMCMC = 10000, BI = 1000, J=50, beta.tilde.true = beta.tilde.true, 
          beta.true = beta.true, alpha.true = alpha.true, string = string0,
          case = cases_matrix[i,1], scen = cases_matrix[i,2], label = label_vec[i],
          seed = seeds[i]-2*n.cores)
}
stopCluster(cl) 

setwd("/home/normington/Paper 2/PosteriorSimResults")
z.star = qnorm(.975)
for(lab in unique(label_vec)) {
  for(scen in 1:4) {
  print(lab)
  print(paste0("Scen ", scen))
  tempBias = read.table(paste0("bias", lab, "Scen", scen, string0, ".txt"))[,1]
  print(paste0("Bias: ", round(mean(tempBias), 3)))
  print(paste0("Bias MoE: ", round(z.star*sd(tempBias)/sqrt(length(tempBias)), 3)))
  
  tempMSE = read.table(paste0("MSE", lab, "Scen", scen, string0, ".txt"))[,1]
  print(paste0("MSE: ", round(mean(tempMSE), 3)))
  print(paste0("MSE MoE: ", round(z.star*sd(tempMSE)/sqrt(length(tempMSE)), 3)))
  
  cov = colMeans(read.table(paste0("Coverage", lab, "Scen", scen, string0, ".txt")))
  print(paste0("Coverage: ", round(cov, 3)))
  }
}
