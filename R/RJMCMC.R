# RJMCMC algorithm using a fixed and a stochastic proposal distribution used 
# for fitting Parula data for regular state space
rjmcmcR <- function(data, inits=NULL, tune=NULL, prior, control = list( n.iters=5000,
                    n.chain=3, n.burn=5000, n.update=2, monitor=FALSE, tol=1e-15 ), 
                    stochastic = TRUE)
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
  # create a set of step size for RJMCMC
  delta <- seq(2, tune$N, 1)
  sx <- c(delta, -delta)
  
  #extract control parameters
  n.iters <- control$n.iters # number of iterations
  n.burn <- control$n.burn # number of burn-in periods
  n.chain <- control$n.chain # number of chains
  n.update <- control$n.update # number of multiple updates on N
  monitor <- control$monitor # FALSE in default; set TRUE to see trace of updating parameters per 100 iterations
  tol <- control$tol # tolerance set for efficient calculation
  
  # check if prior on all parameters are specified otherwise assumed to be uniform
  I1 <- I2 <- I3 <- 1
  if(is.null(prior$sigma)){
    I1 <- 0
    prior$sigma <- c(1,1)
  }
  if(is.null(prior$N)){
    I2 <- 0
    prior$N <- c(1,0.001)
  }
  if(is.null(prior$lam0)){
    I3 <- 0
    prior$lam0 <- c(0, 1)
  }
  # choose the proposal distribution for omitting the activity centres
  # stochastic = TRUE is default for a stochastic removal
  if(stochastic == TRUE){
    omit.s <- function(N, N.cand){
      return(sample(1:N, (N-N.cand), replace = FALSE))
    }
  } else {
    omit.s <- function(N, N.cand){
      return((N.cand+1):N)
    }
  }
  
  result <- mcmc.list()
  niters <-  n.iters + n.burn
  for (j in 1:n.chain) {
    # set initial values for model parameters
    sigma <- inits$sigma
    lam0 <- inits$lam0
    N <- N.cand <- inits$N
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]))
    # calculate initial euclidean distance matrix
    D <- eudistC(S, X)
    # calculate initial detection matrix lambda_ij
    lam <-  K*lamC(D, sigma) 
    # calculate initial log likelihood
    llcur <- sum(dpoisC(n, lam0*colsumC(lam), TRUE))
    
    # matrix to hold samples
    out <- matrix( NA, nrow = niters, ncol = 3 )
    accept.rate <- matrix( NA, nrow = niters, ncol = 3 )
    colnames(out) <- c("sigma", "lam0", "N")
    colnames(accept.rate) <- c("sigma", "lam0", "N")
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    acc.sigma <- 0
    acc.lam0 <- 0
    acc.N <- 0
    
    for(iter in 2:niters) {
      
      if(monitor == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("Acc. S =", Sups/N,"\n")
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
        lam.j.cand <- colsumC(lam.cand)
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
      
      Sups <- 0
      # get the Lambda_j for the trick 
      lam.j <- colsumC(lam)
      # generate proposed activity centres from bivariate normal distrbution
      S.cand <- cbind(rnorm(N, S[,1], tune.S), rnorm(N, S[,2], tune.S))
      # check if the proposed activity centres are within the state space
      inbox <- S.cand[ ,1] >= xlim[1] & S.cand[ ,1] <= xlim[2] &
               S.cand[ ,2] >= ylim[1] & S.cand[ ,2] <= ylim[2] 
      # compute the proposed euclidean distance
      D.cand <- eudistC(S.cand, X)
      # get the lambda_ij matrix for proposed eculidean distance
      lam.cand <- K * lamC(D.cand, sigma) 
      for(i in which(inbox==1)) {
        # implemenet the efficient trick and the rest is the same as above
        lam.j.cand <- lam.j - lam[i,] + lam.cand[i,] + machine.precision
        llcand <- sum(dpoisC(n, lam0 * lam.j.cand, TRUE) ) 
        if(runif(1)< exp(llcand - llcur)) {
          llcur <- llcand
          lam.j <- lam.j.cand
          lam[i,] <- lam.cand[i,]
          S[i,] <- S.cand[i,]
          D[i,] <- D.cand[i,]
          Sups <- Sups+1
        }
      }
    
      
      # updating N for multiple updates 
      for (it in 1:n.update) {
        # generate epsilon
        eps <- sample(sx, 1)
        # generate a new value for N
        N.cand <- N + eps
        # assume uniform(0, M) for N
        # if a new N > the current N, we add new activity centres
        if(N.cand > 0 & N.cand < M){
          if ((N.cand-N) > 0) {
            #generate new activity centres uniformly
            S.cand <- cbind(runif((N.cand-N), xlim[1], xlim[2]),
                            runif((N.cand-N), ylim[1], ylim[2]))
            # update all matrix, euclidean matrix, lambda_ij for new activity centres
            D.cand <- eudistC(S.cand, X)
            lam.cand <- K * lamC(D.cand, sigma)
            # combine the current lambda_ij and the new one
            lam.cand <- rbind(lam, lam.cand)
            # compute the likelihood for the new values
            llcand <- sum(dpoisC(n, lam0 * colsumC(lam.cand), TRUE))
            
            if(runif(1) < exp(llcand  - llcur  + llN.cand - llN)){
              llcur <- llcand
              lam <- lam.cand
              N <- N.cand
              S <- rbind(S,S.cand)
              D <- rbind(D,D.cand)
              acc.N <- acc.N + 1
            }
          } else {
            # omit the activity centres 
            # where omit.s depends on the stochastic (TRUE or not)
            id <- omit.s(N, N.cand)
            # substract the selected individuals
            lam.cand <- matrix(lam[-id , ], nrow = N.cand)
            llcand <- sum(dpoisC(n, lam0 * colsumC(lam.cand), TRUE)) 
            
            if(runif(1) < exp(llcand  - llcur  + llN.cand - llN)){
              llcur <- llcand
              lam <- lam.cand
              N <- N.cand
              D <- matrix(D[ -id ,], nrow = N.cand)
              S <- matrix(S[ -id ,], nrow = N.cand)
              acc.N <- acc.N + 1
            }
          }
        }
        setTxtProgressBar(pb,iter)    
        out[iter,] <- c(sigma, lam0, N)
        accept.rate[iter,] <- c(acc.sigma/iter, acc.lam0/iter, acc.N/iter)
      }
    }
    close(pb)
    out.mcmc <- mcmc(out[-(1:n.burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}

# RJMCMC algorithm assuming stochastic proposal distribution used for fitting barking deer 
# for irregular state space
rjmcmcIR <- function(data, inits=NULL, tune=NULL, prior, control = list( n.iters=5000,
                     n.chain=3, n.burn=5000, n.update=2, monitor=FALSE, tol=1e-15 ) )
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
  # create a set of step size for RJMCMC
  delta <- seq(2, tune$N, 1)
  sx <- c(delta, -delta)
  
  #extract control parameters
  n.iters <- control$n.iters # number of iterations
  n.burn <- control$n.burn # number of burn-in periods
  n.chain <- control$n.chain # number of chains
  n.update <- control$n.update # number of multiple updates on N
  monitor <- control$monitor # FALSE in default; set TRUE to see trace of updating parameters per 100 iterations
  tol <- control$tol # tolerance set for efficient calculation
  
  # check if prior on all parameters are specified otherwise assumed to be uniform
  I1 <- I2 <- I3 <- 1
  if(is.null(prior$sigma)){
    I1 <- 0
    prior$sigma <- c(1,1)
  }
  if(is.null(prior$N)){
    I2 <- 0
    prior$N <- c(1,0.001)
  }
  if(is.null(prior$lam0)){
    I3 <- 0
    prior$lam0 <- c(0, 1)
  }
  
  result <- mcmc.list()
  niters <-  n.iters + n.burn
  
  for (ch in 1:n.chain) {
    # set initial values for model parameters
    sigma <- inits$sigma
    lam0 <- inits$lam0
    N <- N.cand <- inits$N
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]))
    # check if initial values for activity centres lie within the state space
    inbox <- inside.owin( S[,1], S[,2], wind )
    for (i in 1:N) {
      if(!inbox[i]){
        while (!inbox[i]) {
          S[i,1:2] <- c(runif(1, xlim[1], xlim[2]),
                        runif(1, ylim[1], ylim[2]))
          inbox[i] <- inside.owin( S[i,1], S[i,2], wind )
        }
      }
    }
    # calculate initial euclidean distance matrix
    D <- eudistC(S, X)
    # calculate initial detection matrix lambda_ij
    lam <-  K*lamC(D, sigma) 
    # calculate initial vector Lambda_j by summing over individual i 
    lam.j <- colsumC(lam)
    # calculate initial log likelihood
    llcur <- sum(dpoisC(n, lam0*lam.j, TRUE))
    
    # matrix to hold samples
    out <- matrix( c(sigma, lam0, N), nrow = niters, ncol = 3)
    # out.S <- S
    accept.rate <- matrix( NA, nrow = niters, ncol = 3)
    colnames(out) <- c("sigma", "lam0", "N")
    colnames(accept.rate) <- c("sigma", "lam0", "N")
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    
    acc.sigma <- 0
    acc.lam0 <- 0
    acc.N <- 0
    for(iter in 2:niters) {
      
      if(monitor == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("Acc. S =", Sups/N,"\n")
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
        lam.j.cand <- colsumC(lam.cand)
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
      
      
      # update S
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
        lam.j.cand <- lam.j - lam[i,] + lam.cand[i,] + tol
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
        if(runif(1)< exp(llcand - llcur )) { 
          llcur <- llcand
          lam[i,] <- lam.cand[i,]
          lam.j <- lam.j.cand
          S[i,] <- S.cand[i,]
          D[i,] <- D.cand[i,]
          Sups <- Sups+1
        }
      }
      
      # updating N for multiple updates 
      for (it in 1:n.update) {
        # generate epsilon
        eps <- sample(sx, 1)
        # generate a new value for N
        N.cand <- N + eps
        # compute the prior assuming Negative-Binomial
        llN <- dnbinom(N, size = prior$N[1], prob = prior$N[2], log = TRUE) * I2
        llN.cand <- dnbinom(N.cand, size = prior$N[1], prob = prior$N[2], log = TRUE) * I2
        # if a new N > the current N, we add new activity centres
        if(N.cand > 0 ){
          if ((N.cand-N) > 0) {
            #generate new activity centres uniformly
            S.cand <- cbind(runif((N.cand-N), xlim[1], xlim[2]),
                            runif((N.cand-N), ylim[1], ylim[2]))
            # we need to check if new activity centres are in the state space
            inbox <- inside.owin(S.cand[,1], S.cand[,2], wind)
            for (i in 1:(N.cand-N)) {
              if(!inbox[i]){
                while (!inbox[i]) {
                  S.cand[i,1:2] <- c(runif(1, xlim[1], xlim[2]),
                                     runif(1, ylim[1], ylim[2]))
                  inbox[i] <- inside.owin(S.cand[i,1], S.cand[i,2], wind)
                }
              }
            }
            # update all matrix, euclidean matrix, lambda_ij for new activity centres
            D.cand <- eudistC(S.cand, X)
            lam.cand <-  K*lamC(D.cand, sigma) 
            # adding them to the old Lambda_j
            lam.j.cand <- lam.j + colsumC(lam.cand)
            # compute the likelihood for the new values
            llcand <- sum( dpoisC(n, lam0*lam.j.cand, TRUE))
            
            if(runif(1) < exp(llcand  - llcur + llN.cand - llN)){ 
              llcur <- llcand
              lam.j <- lam.j.cand
              lam <- rbind(lam,lam.cand)
              N <- N.cand
              S <- rbind(S,S.cand)
              D <- rbind(D,D.cand)
              acc.N <- acc.N + 1
            }
          # if a new N < the current N we omit activity centres randomly
          } else {
            # omit the activity centres stochastically
            ne <- sample(1:N, (N-N.cand), replace = FALSE)
            # extract lambda_ij that we omit
            lam.cand <- matrix(lam[ne, ], ncol = J)
            # substract the current Lambda_j with the omitted Lambda_j
            lam.j.cand <- lam.j - colsumC(lam.cand) + tol
            llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
            
            if(runif(1) < exp(llcand  - llcur + llN.cand - llN)){
              llcur <- llcand
              lam.j <- lam.j.cand
              lam <- lam[-ne, ]
              N <- N.cand
              D <- D[-ne, ]
              S <- S[-ne, ]
              acc.N <- acc.N + 1
            }
          }
        }
      }
      
      setTxtProgressBar(pb,iter)
      
      out[iter,] <- c(sigma, lam0, N)
      accept.rate[iter,] <- c(acc.sigma/iter, acc.lam0/iter, acc.N/iter)
      
    }
    close(pb)
    out.mcmc <- mcmc(out[-(1:n.burn),])
    result[[ch]] <- out.mcmc
    # result[[ n.chain + ch ]] <- out.S
  }
  return(result)
}

# RJMCMC algorithm assuming stochastic proposal distribution used for fitting barking deer 
# for irregular state space with trap-level covariates.
# Additonal RJMCMC is performed for model selection
# The function is basically similar to the previous RJMCMC for irregular state space
rjmcmcIRH <- function(data, inits=NULL, tune=NULL, prior, control = list(n.iters=5000,
                      n.chain=3, n.burn=5000, n.update=2, monitor=FALSE, tol=1e-15))
{
  #extract data
  n <- data$n
  K <- data$K
  X <- data$X
  J <- dim(X)[1]
  xlim <- data$xlim
  ylim <- data$ylim
  wind <- data$wind
  # adding covariate
  habitat <- data$habitat
  
  #extract tuning paramaters
  tune.sigma <- tune$sigma
  tune.beta <- tune$beta
  tune.S <- tune$S
  delta <- seq(2, tune$N, 1)
  sx <- c(delta, -delta)
  
  #extract control parameters
  n.iters <- control$n.iters
  n.burn <- control$n.burn
  n.chain <- control$n.chain
  n.update <- control$n.update
  monitor <- control$monitor
  tol <- control$tol
  
  I1 <- I2 <- 1
  if(is.null(prior$sigma)){
    I1 <- 0
    prior$sigma <- c(1,1)
  }
  if(is.null(prior$N)){
    I2 <- 0
    prior$N <- c(1, 0.001)
  }
  
  prior.beta <- prior$beta
  # the number of parameters in log(lambda_0) including the constant beta_0 (P + 1)
  P <- length(inits$beta)
  # indicator variable for model selection
  # 1 = be included; 0 otherwise
  IH <- rep(1, P)
  
  # initial values
  result <- mcmc.list()
  niters <-  n.iters + n.burn
  
  for (ch in 1:n.chain) {
    sigma <- inits$sigma
    N <- N.cand <- inits$N
    beta <- inits$beta
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]))
    inbox <- inside.owin(S[,1], S[,2], wind)
    for (i in 1:N) {
      if(!inbox[i]){
        while (!inbox[i]) {
          S[i,1:2] <- c(runif(1, xlim[1], xlim[2]),
                        runif(1, ylim[1], ylim[2]))
          inbox[i] <- inside.owin(S[i,1], S[i,2], wind)
        }
      }
    }
    # create covariate matrix with a first column corresponding to the constant
    covar <- model.matrix( ~ as.matrix(habitat))
    # obtain the lambda_0 i.e. lambda_0 = exp(beta_0 + beta_1 X_1 + ... + beta_P X_P)
    lam0 <- exp(covar %*% (beta * IH))
    D <- eudistC(S, X)
    lam <-  K*lamC(D, sigma) 
    lam.j <- colsumC(lam)
    llcur <- sum(dpoisC(n, lam0*lam.j, TRUE))
    
    # matrix to hold samples
    out <- matrix(c(sigma, beta, N, IH), nrow = niters, ncol = 2*P + 2)
    accept.rate <- matrix(NA,  nrow = niters, ncol = 2*P + 2)
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    
    acc.beta <- rep(0, P)
    acc.sigma <- 0
    acc.N <- 0
    acc.model <- 0
    iter.beta <- 0
    
    for(iter in 2:niters) {      
      
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
        lam.j.cand <- colsumC(lam.cand)
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
      
      # update S
      # update S
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
        lam.j.cand <- lam.j - lam[i,] + lam.cand[i,] + tol
        llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
        if(runif(1)< exp(llcand - llcur )) { 
          llcur <- llcand
          lam[i,] <- lam.cand[i,]
          lam.j <- lam.j.cand
          S[i,] <- S.cand[i,]
          D[i,] <- D.cand[i,]
          Sups <- Sups+1
        }
      }
     
      
      # updating N for multiple updates 
      for (it in 1:n.update) {
        # generate epsilon
        eps <- sample(sx, 1)
        # generate a new value for N
        N.cand <- N + eps
        # compute the prior assuming Negative-Binomial
        llN <- dnbinom(N, size = prior$N[1], prob = prior$N[2], log = TRUE) * I2
        llN.cand <- dnbinom(N.cand, size = prior$N[1], prob = prior$N[2], log = TRUE) * I2
        # if a new N > the current N, we add new activity centres
        if(N.cand > 0 ){
          if ((N.cand-N) > 0) {
            #generate new activity centres uniformly
            S.cand <- cbind(runif((N.cand-N), xlim[1], xlim[2]),
                            runif((N.cand-N), ylim[1], ylim[2]))
            # we need to check if new activity centres are in the state space
            inbox <- inside.owin(S.cand[,1], S.cand[,2], wind)
            for (i in 1:(N.cand-N)) {
              if(!inbox[i]){
                while (!inbox[i]) {
                  S.cand[i,1:2] <- c(runif(1, xlim[1], xlim[2]),
                                     runif(1, ylim[1], ylim[2]))
                  inbox[i] <- inside.owin(S.cand[i,1], S.cand[i,2], wind)
                }
              }
            }
            # update all matrix, euclidean matrix, lambda_ij for new activity centres
            D.cand <- eudistC(S.cand, X)
            lam.cand <-  K*lamC(D.cand, sigma) 
            # adding them to the old Lambda_j
            lam.j.cand <- lam.j + colsumC(lam.cand)
            # compute the likelihood for the new values
            llcand <- sum( dpoisC(n, lam0*lam.j.cand, TRUE))
            
            if(runif(1) < exp(llcand  - llcur + llN.cand - llN)){ 
              llcur <- llcand
              lam.j <- lam.j.cand
              lam <- rbind(lam,lam.cand)
              N <- N.cand
              S <- rbind(S,S.cand)
              D <- rbind(D,D.cand)
              acc.N <- acc.N + 1
            }
            # if a new N < the current N we omit activity centres randomly
          } else {
            # omit the activity centres stochastically
            ne <- sample(1:N, (N-N.cand), replace = FALSE)
            # extract lambda_ij that we omit
            lam.cand <- matrix(lam[ne, ], ncol = J)
            # substract the current Lambda_j with the omitted Lambda_j
            lam.j.cand <- lam.j - colsumC(lam.cand) + tol
            llcand <- sum(dpoisC(n, lam0*lam.j.cand, TRUE))
            
            if(runif(1) < exp(llcand  - llcur + llN.cand - llN)){
              llcur <- llcand
              lam.j <- lam.j.cand
              lam <- lam[-ne, ]
              N <- N.cand
              D <- D[-ne, ]
              S <- S[-ne, ]
              acc.N <- acc.N + 1
            }
          }
        }
      }
      
      #update beta in the current model
      index <- which(IH == 1)
      for ( p in index ) {
        beta.cand <- beta
        # propose a new value for beta_p
        beta.cand[p] <- beta[p] + tune.beta[p]*rnorm(1)
        # compute the prior distribution assumed to be normal
        ll.beta.cand <- dnorm(beta.cand[p], mean = 0, sd = prior.beta[p], TRUE)
        ll.beta <- dnorm(beta[p], mean = 0, sd = prior.beta[p], TRUE)
        
        # update the lambda_0 vector
        lam0.cand <- exp(covar %*% (beta.cand * IH))
        llcand <- sum(dpoisC(n, lam0.cand*lam.j, TRUE)) 
        
        if(runif(1) < exp( llcand  - llcur + ll.beta.cand - ll.beta ) ){ 
          llcur <- llcand
          lam0 <- lam0.cand
          beta <- beta.cand
          acc.beta[p] <- acc.beta[p] + 1
        }
      }
      
      # update the current model
      for ( p in 2:P ) {
        IH.cand <- IH
        beta.cand <- beta
        # propose a new model either: from 1 to 0 (remove the current parameter) or 
        # from 0 to 1 (adding a new parameter)
        IH.cand[p] <- ifelse(IH.cand[p] == 0, 1, 0)
        # if the proposed model is to remove the current parameter
        # then we simply need to compute the prior and proposal distribution
        if( IH.cand[p] == 0) {
          prop.beta <- dnorm(beta[p], mean = 0, sd = tune.beta[p], TRUE) 
          prio.beta <- dnorm(beta[p], mean = 0, sd = prior.beta[p], TRUE) 
          lam0.cand <- exp(covar %*% (beta * IH.cand))
          llcand <- sum( dpoisC(n, lam0.cand*lam.j, TRUE) )
        } else {
          # we need to propose a new value for a new parameter
          beta.cand[p] <- rnorm(1, mean = 0, sd = tune.beta[p])
          # compute the prior and proposal distribution
          prop.beta <- dnorm(beta.cand[p], mean = 0, sd = tune.beta[p], TRUE) 
          prio.beta <- dnorm(beta.cand[p], mean = 0, sd = prior.beta[p], TRUE)
          lam0.cand <- exp( covar %*% (beta.cand * IH.cand) )
        }
        
        llcand <- sum( dpoisC(n, lam0.cand*lam.j, TRUE))
        
        if(runif(1) < exp( lcand  - llcur + prio.beta - prop.beta)){
          llcur <- llcand
          lam0 <- lam0.cand
          IH <- IH.cand
          beta <- beta.cand
        }
      }
      
      setTxtProgressBar(pb,iter)    
      out[ iter , ] <- c(sigma, beta, N, IH)
      accept.rate[ iter ,] <- c(acc.sigma / iter, acc.beta / iter, acc.N / iter, IH )
    }
    close(pb)
    out.mcmc <- mcmc(out[-(1:n.burn),])
    result[[ch]] <- out.mcmc
  }
  return(result)
}
