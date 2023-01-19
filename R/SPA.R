# Super Population Approach (SPA) for regular state space used for fitting Parula data
spaR <- function(data, inits=NULL, tune=NULL, prior, control = list( n.iters=5000,
                                                                     n.chain=3, n.burn=5000, n.update=2, monitor=FALSE, tol=1e-15))
{
  # extract data
  # observed data n_j summed over times
  n <- data$n
  # total sampling times
  K <- data$K
  # trap locations 
  X <- data$X
  J <- dim(X)[1]
  # x-limit
  xlim <- data$xlim
  # y-limit
  ylim <- data$ylim
  # upper limit of population
  M <- data$M
  
  #extract tuning paramaters
  tune.sigma <- tune$sigma
  tune.lam0 <- tune$lam0
  tune.S <- tune$S
  
  #extract control parameters
  n.iters <- control$n.iters # number of iterations
  n.burn <- control$n.burn # number of burn-in periods
  n.chain <- control$n.chain # number of chains
  monitor <- control$monitor # FALSE in default; set TRUE to see trace of updating parameters per 100 iterations
  tol <- control$tol # tolerance set for efficient calculation
  
  
  # check if prior on sigma is specified otherwise assumed to be uniform
  I1 <- I2 <- 1
  if(is.null(prior$sigma)){
    I1 <- 0
    prior$sigma <- c(1,1)
  }
  
  result <- mcmc.list()
  niters <-  n.iters + n.burn
  for (j in 1:n.chain) {
    # set initial values for model parameters
    sigma <- inits$sigma
    lam0 <- inits$lam0
    S <- cbind(runif(M, xlim[1], xlim[2]),
               runif(M, ylim[1], ylim[2]))
    w <- inits$z
    psi <- inits$psi
    # calculate initial euclidean distance matrix
    D <- eudistC(S, X)
    # calculate initial detection matrix lambda_ij
    lam <- K*lamC(D, sigma)*w 
    # calculate initial vector Lambda_j by summing over individual i 
    lam.j <- colsumC(lam*w)
    # calculate initial log likelihood
    llcur <- sum(dpoisC(n, lam0*lam.j, TRUE))
    
    # matrix to hold samples
    out <- matrix(NA, nrow = niters, ncol = 4)
    accept.rate <- matrix( NA, nrow = niters, ncol = 2)
    colnames(out) <- c("sigma", "lam0", "N", "psi")
    colnames(accept.rate) <- c("sigma", "lam0")
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    acc.sigma <- 0
    acc.lam0 <- 0
    
    for(iter in 2:niters) {
      
      if(monitor == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("Acc. S =", Sups/M,"\n")
          cat("Acc. W =", Wups/M,"\n")
        }
      }
      # update sigma
      sigma.cand <- sigma + tune.sigma*rnorm(1)
      
      if(sigma.cand > 0){
        ldsigma <- dgamma(sigma, prior$sigma[1], prior$sigma[2], log = TRUE) * I1
        ldsigma.cand <- dgamma(sigma.cand, prior$sigma[1], prior$sigma[2], log = TRUE) * I1
        
        lam.cand <- K*lamC(D, sigma.cand)
        lam.j.cand <- colsumC(lam.cand*w)
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE)) 
        
        if(runif(1) < exp(llcand - llcur + ldsigma.cand - ldsigma)){
          llcur <- llcand
          lam <- lam.cand
          lam.j <- lam.j.cand
          sigma <- sigma.cand
          acc.sigma <- acc.sigma + 1
        }
      }
      # update sigma
      sigma.cand <- sigma + tune.sigma*rnorm(1)
      
      # move to the next step if the proposed sigma <= 0
      if (sigma.cand > 0) {
        # compute the prior distributions for sigma assuming gamma prior
        ldsigma <- dgamma(sigma, prior$sigma[1], prior$sigma[2], log = TRUE)*I1
        ldsigma.cand <- dgamma(sigma.cand, prior$sigma[1], prior$sigma[2], log = TRUE)*I1
        #ge the lambda_ij matrix for a new sigma
        lam.cand <- K*lamC(D, sigma.cand)
        # update Lambda_j by summing the lambda_ij candidate over individuals
        lam.j.cand <- colsumC(lam.cand*w)
        # compute the log-likelihood for a proposed value
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE)) 
        
        # rejection criteria in Metropolis Hastings
        if(runif(1) < exp(llcand - llcur)){
          llcur <- llcand
          lam <- lam.cand
          lam.j <- lam.j.cand
          sigma <- sigma.cand
          acc.sigma <- acc.sigma + 1
        }
      }
      
      # update lam0 assuming uniform prior
      # propose a new value
      lam0.cand <- lam0 + tune.lam0*rnorm(1)
      # move to the next step if the proposed lambda_0 <= 0
      if(lam0.cand > 0){
        # compute the log-likelihood for a proposed value
        llcand <- sum(dpoisC(n, lam0.cand*lam.j, TRUE)) 
        if(runif(1) < exp(llcand  - llcur)){
          llcur <- llcand
          lam0 <- lam0.cand
          acc.lam0 <- acc.lam0 + 1
        }
      }
      
      
      # update w
      Wups <- 0
      for(i in 1:M) {
        w.cand <- w
        w.cand[i] <- ifelse(w[i]==0, 1, 0)
        # we implement the efficient likelihood calculation
        # substracting and adding the contribution of the corresponding updated indvidual to the current Lambda_j
        lam.j.cand <- lam.j + lam[i,]*(w.cand[i] - w[i]) + tol
        # compute the log-likelihood for the proposed values
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
        # compute the prior distrbutions for psi
        prior.psi <- dbinom(w[i], 1, psi, log=TRUE)
        prior.psi.cand <- dbinom(w.cand[i], 1, psi, log=TRUE)
        if(runif(1) < exp(llcand + prior.psi.cand - llcur - prior.psi)) {
          w <- w.cand
          lam.j <- lam.j.cand
          llcur <- llcand
          Wups <- Wups+1
        }
      }
      
      # update psi
      psi <- rbeta(1, 1+sum(w), 1+M-sum(w))
      
      # Update S
      Sups <- 0
      # generate proposed activity centres from bivariate normal distrbution
      S.cand <- cbind(rnorm(M, S[,1], tune.S), rnorm(M, S[,2], tune.S))
      # check if the proposed activity centres are within the state space
      inbox <- S.cand[ ,1] >= xlim[1] & S.cand[ ,1] <= xlim[2] &
        S.cand[ ,2] >= ylim[1] & S.cand[ ,2] <= ylim[2]
      # compute the proposed euclidean distance
      D.cand <- eudistC(S.cand, X)
      # get the lambda_ij matrix for proposed eculidean distance
      lam.cand <- K*lamC(D.cand, sigma)
      for(i in which(inbox)) {
        # implemenet the efficient trick and the rest is the same as above
        lam.j.cand <- lam.j + (lam.cand[i,] - lam[i,])*w[i] + tol
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
        if(runif(1) < exp(llcand - llcur)) {
          llcur <- llcand
          lam.j <- lam.j.cand
          lam[i,] <- lam.cand[i,]
          S[i,] <- S.cand[i,]
          D[i,] <- D.cand[i,]
          Sups <- Sups+1
        }
      }
      
      setTxtProgressBar(pb,iter)    
      out[iter,] <- c(sigma, lam0, sum(w), psi)
      accept.rate[iter,] <- c(acc.sigma/iter, acc.lam0/iter)
    }
    close(pb)
    out.mcmc <- mcmc(out[-(1:n.burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}


# Super Population Approach (SPA) for irregular state space for fitting barking deer
spaIR <- function(data, inits=NULL, tune=NULL, prior, control = list( n.iters=5000,
                  n.chain=3, n.burn=5000, n.update=2, monitor=FALSE, tol=1e-15))
{
  # extract data
  # observed data n_j summed over times
  n <- data$n
  # total sampling times
  K <- data$K
  # trap locations 
  X <- data$X
  J <- dim(X)[1]
  # x-limit
  xlim <- data$xlim
  # y-limit
  ylim <- data$ylim
  # observation window 
  wind <- data$wind
  # upper limit of population
  M <- data$M
  
  #extract tuning paramaters
  tune.sigma <- tune$sigma
  tune.lam0 <- tune$lam0
  tune.S <- tune$S
  
  #extract control parameters
  n.iters <- control$n.iters # number of iterations
  n.burn <- control$n.burn # number of burn-in periods
  n.chain <- control$n.chain # number of chains
  monitor <- control$monitor # FALSE in default; set TRUE to see trace of updating parameters per 100 iterations
  tol <- control$tol # tolerance set for efficient calculation
  

  # check if prior on lambda_0 is specified otherwise assumed to be uniform
  if(is.null(prior$lam0)){
    I2 <- 0
    prior$lam0 <- c(0, 1)
  }
  
  result <- mcmc.list()
  niters <-  n.iters + n.burn
  for (j in 1:n.chain) {
    # set initial values for model parameters
    sigma <- inits$sigma
    lam0 <- inits$lam0
    S <- cbind(runif(M, xlim[1], xlim[2]),
               runif(M, ylim[1], ylim[2]))
    # check if initial values for activity centres lie within the state space
    inbox <- inside.owin(S[,1], S[,2], wind)
    for (i in 1:M) {
      if(!inbox[i]){
        while (!inbox[i]) {
          S[i,1:2] <- c(runif(1, xlim[1], xlim[2]),
                        runif(1, ylim[1], ylim[2]))
          inbox[i] <- inside.owin(S[i,1], S[i,2], wind)
        }
      }
    }
    w <- inits$z
    psi <- inits$psi
    # calculate initial euclidean distance matrix
    D <- eudistC(S, X)
    # calculate initial detection matrix lambda_ij
    lam <- K*lamC(D, sigma)*w 
    # calculate initial vector Lambda_j by summing over individual i 
    lam.j <- colsumC(lam*w)
    # calculate initial log likelihood
    llcur <- sum(dpoisC(n, lam0*lam.j, TRUE))
    
    # matrix to hold samples
    out <- matrix(NA, nrow = niters, ncol = 4)
    accept.rate <- matrix( NA, nrow = niters, ncol = 2)
    colnames(out) <- c("sigma", "lam0", "N", "psi")
    colnames(accept.rate) <- c("sigma", "lam0")
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    acc.sigma <- 0
    acc.lam0 <- 0
    
    for(iter in 2:niters) {
      
      if(monitor == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("Acc. S =", Sups/M,"\n")
          cat("Acc. W =", Wups/M,"\n")
        }
      }
      
      # update sigma
      sigma.cand <- sigma + tune.sigma*rnorm(1)
      
      # move to the next step if the proposed sigma <= 0
      if (sigma.cand > 0) {
        #ge the lambda_ij matrix for a new sigma
        lam.cand <- K*lamC(D, sigma.cand)
        # update Lambda_j by summing the lambda_ij candidate over individuals
        lam.j.cand <- colsumC(lam.cand*w)
        # compute the log-likelihood for a proposed value
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE)) 
        
        # rejection criteria in Metropolis Hastings
        if(runif(1) < exp(llcand - llcur)){
          llcur <- llcand
          lam <- lam.cand
          lam.j <- lam.j.cand
          sigma <- sigma.cand
          acc.sigma <- acc.sigma + 1
        }
      }
      
      # update lam0
      # we assume a normal prior on log lambda_0 as used in the paper
      
      # propose a new value in log scale
      log.lam0.cand <- log(lam0) + tune.lam0*rnorm(1) 
      lam0.cand <- exp(log.lam0.cand) 
      # compute the log-likelihood for a proposed value
      llcand <- sum(dpoisC(n, lam0.cand*lam.j, TRUE)) 
      # compute the prior distributions
      lllam0 <- dnorm(log(lam0), prior$lam0[1], prior$lam0[2], log = TRUE)*I2 
      lllam0.cand <- dnorm(log.lam0.cand, prior$lam0[1], prior$lam0[2], log = TRUE)*I2
      
      # rejection criteria in Metropolis Hastings
      if(runif(1) < exp(llcand - llcur + lllam0.cand - lllam0) ){ 
        llcur <- llcand
        lam0 <- lam0.cand
        acc.lam0 <- acc.lam0 + 1
      }
      
      # update w
      Wups <- 0
      for(i in 1:M) {
        w.cand <- w
        w.cand[i] <- ifelse(w[i]==0, 1, 0)
        # we implement the efficient likelihood calculation
        # substracting and adding the contribution of the corresponding updated indvidual to the current Lambda_j
        lam.j.cand <- lam.j + lam[i,]*(w.cand[i] - w[i]) + tol
        # compute the log-likelihood for the proposed values
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
        # compute the prior distrbutions for psi
        prior.psi <- dbinom(w[i], 1, psi, log=TRUE)
        prior.psi.cand <- dbinom(w.cand[i], 1, psi, log=TRUE)
        if(runif(1) < exp(llcand + prior.psi.cand - llcur - prior.psi)) {
          w <- w.cand
          lam.j <- lam.j.cand
          llcur <- llcand
          Wups <- Wups+1
        }
      }
      
      # update psi
      psi <- rbeta(1, 1+sum(w), 1+M-sum(w))
      
      # Update S
      Sups <- 0
      # generate proposed activity centres from bivariate normal distrbution
      S.cand <- cbind(rnorm(M, S[,1], tune.S), rnorm(M, S[,2], tune.S))
      # check if the proposed activity centres are within the state space
      inbox <- inside.owin(S.cand[,1], S.cand[,2], wind)
      # compute the proposed euclidean distance
      D.cand <- eudistC(S.cand, X)
      # get the lambda_ij matrix for proposed eculidean distance
      lam.cand <- K*lamC(D.cand, sigma)
      
      for(i in which(inbox)) {
        # implemenet the efficient trick and the rest is the same as above
        lam.j.cand <- lam.j + (lam.cand[i,] - lam[i,])*w[i] + tol
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
        if(runif(1)< exp(llcand - llcur)) {
          llcur <- llcand
          lam.j <- lam.j.cand
          lam[i,] <- lam.cand[i,]
          S[i,] <- S.cand[i,]
          D[i,] <- D.cand[i,]
          Sups <- Sups+1
        }
      }
      
      # setTxtProgressBar(pb,iter)    
      out[iter,] <- c(sigma, lam0, sum(w), psi)
      accept.rate[iter,] <- c(acc.sigma/iter, acc.lam0/iter)
    }
    # close(pb)
    out.mcmc <- mcmc(out[-(1:n.burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}