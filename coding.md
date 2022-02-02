coding
================
Riki Herliansyah
02/02/2022

## Standard RJMCMC

``` r
srjmcmc <- function(n, X, K=NULL, niters=5000, xlim, ylim, tune=c(0.1, 0.1, 2, 5), prior=NULL,
                   monitorS=FALSE, chain=3, burn=5000, start.values=c(0.1,0.1,60), update=FALSE,
                   thin=2)
{
  K <- dim(n)[2]
  ndot <- rowSums(n)
  result <- mcmc.list()
  niters <- niters + burn
  delta <- seq(1,tune[4],1)
  sx <- c(delta,-delta)
  
  # initial values
  for (j in 1:chain) {
    sigma <- runif(1,0.1, start.values[1])
    lam0 <- runif(1, 0.02, start.values[2])
    N <- N.cand <- rpois(1,start.values[3])
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]))
    D <- e2dist(S, X)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    llcur <- sum(dpois(ndot, K*colSums(lam), log=TRUE))
    
    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=3)
    accept.rate <- matrix(NA, nrow=niters, ncol=3)
    colnames(out) <- c("sigma", "lam0", "N")
    colnames(accept.rate) <- c("sigma", "lam0", "N")
    check <- rep(0, niters)
    
    #cat("\n initial values =", c(sigma, lam0, N, "\n"))
    #cat("\n chain = ", j, "\n")
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    ns <- nl <- nn <-0
    for(iter in 2:niters) {
      if(update == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("  Accept S\n")
          cat("    S =", Sups/N, "\n")
        }
      }
      
      # update sigma
      sigma.cand <- rnorm(1, sigma, tune[1])
      
      if(is.null(prior))
      { 
        dSigma <- dunif(sigma, 0, 100)
        dSigma.cand <- dunif(sigma.cand, 0, 100)
      } else {
        dSigma <- dgamma(sigma, prior[1], prior[2])
        dSigma.cand <- dgamma(sigma.cand, prior[1], prior[2])
      }
      if (sigma.cand > 0){
        lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
        llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
        
        if(runif(1) < exp( llcand - llcur ) ){
          llcur <- llcand
          lam <- lam.cand
          sigma <- sigma.cand
          ns <- ns + 1
        }
      }
      
      # update lam0
      lam0.cand <- rnorm(1, lam0, tune[2])  
      if(lam0.cand>0) {
        lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
        llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
        
        if(runif(1) < exp( llcand  - llcur ) ){
          llcur <- llcand
          lam <- lam.cand
          lam0 <- lam0.cand
          nl <- nl + 1
        }
      }
      
      # update S
      if (iter %% thin ==0) {
        Sups <- 0
        col.lam <- colSums(lam)
        Scand <-cbind(rnorm(N, S[,1], tune[3]), rnorm(N, S[,2], tune[3]))
        
        #check if the new activity centers within the boundary
        inbox <- inside.owin(Scand[,1], Scand[,2], wind)
        dtmp <- e2dist(Scand, X)
        lam.cand <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
        for(i in 1:N) {
          if(!inbox[i])
            next
          c.cand <- abs(col.lam-lam[i,]+lam.cand[i,])
          llcand <- sum(dpois(ndot, K*c.cand, log=TRUE) ) 
          if(runif(1)< exp(llcand - llcur)) {
            llcur <- llcand
            lam[i,] <- lam.cand[i,]
            col.lam <- c.cand
            S[i,] <- Scand[i,]
            D[i,] <- dtmp[i,]
            Sups <- Sups+1
          }
        }
      }
      
      #update N 
      #generate epsilon
      eps <- sample(sx, 1)
      N.cand <- N + eps
      
      if ((N.cand-N) > 0) {
          #generate new activity centres uniformly
          S.cand <- cbind(runif((N.cand-N), xlim[1], xlim[2]),
                          runif((N.cand-N), ylim[1], ylim[2]))
          S.cand <-rbind(S, S.cand)
          D.cand <- e2dist(S.cand, X)
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE))
          
          if(runif(1) < exp( (llcand)  - (llcur)) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            S <- S.cand
            D <- D.cand
            nn <- nn + 1
          }
        } else {
          #omit the activity centres deterministically
          S.cand <- S[1:N.cand,]
          D.cand <- D[1:N.cand,]
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
          
          if(runif(1) < exp( (llcand)  - (llcur)) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            D <- D.cand
            S <- S.cand
            nn <- nn + 1
          }
      }
      setTxtProgressBar(pb,iter)    
      out[iter,] <- c(sigma, lam0, N)
      accept.rate[iter,] <- c(ns/iter, nl/iter, nn/iter)
    }
    close(pb)
    out.mcmc <- mcmc(out[-(1:burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}
```

## Random RJMCMC

``` r
rrjmcmc <- function(n, X, K=NULL, niters=5000, xlim, ylim, tune=c(0.1, 0.1, 2, 5), prior=NULL,
                   monitorS=FALSE, chain=3, burn=5000, start.values=c(0.1,0.1,60), update=FALSE,
                   thin=2)
{
  K <- dim(n)[2]
  ndot <- rowSums(n)
  result <- mcmc.list()
  niters <- niters + burn
  delta <- seq(1,tune[4],1)
  sx <- c(delta,-delta)
  
  # initial values
  for (j in 1:chain) {
    sigma <- runif(1,0.1, start.values[1])
    lam0 <- runif(1, 0.02, start.values[2])
    N <- N.cand <- rpois(1,start.values[3])
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]))
    D <- e2dist(S, X)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    llcur <- sum(dpois(ndot, K*colSums(lam), log=TRUE))
    
    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=3)
    accept.rate <- matrix(NA, nrow=niters, ncol=3)
    colnames(out) <- c("sigma", "lam0", "N")
    colnames(accept.rate) <- c("sigma", "lam0", "N")
    check <- rep(0, niters)
    
    #cat("\n initial values =", c(sigma, lam0, N, "\n"))
    #cat("\n chain = ", j, "\n")
    
    pb = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    ns <- nl <- nn <-0
    for(iter in 2:niters) {
      if(update == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("  Accept S\n")
          cat("    S =", Sups/N, "\n")
        }
      }
      
      # update sigma
      sigma.cand <- rnorm(1, sigma, tune[1])
      
      if(is.null(prior))
      { 
        dSigma <- dunif(sigma, 0, 100)
        dSigma.cand <- dunif(sigma.cand, 0, 100)
      } else {
        dSigma <- dgamma(sigma, prior[1], prior[2])
        dSigma.cand <- dgamma(sigma.cand, prior[1], prior[2])
      }
      if (sigma.cand > 0){
        lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
        llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
        
        if(runif(1) < exp( llcand - llcur ) ){
          llcur <- llcand
          lam <- lam.cand
          sigma <- sigma.cand
          ns <- ns + 1
        }
      }
      
      # update lam0
      lam0.cand <- rnorm(1, lam0, tune[2])  
      if(lam0.cand>0) {
        lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
        llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
        
        if(runif(1) < exp( llcand  - llcur ) ){
          llcur <- llcand
          lam <- lam.cand
          lam0 <- lam0.cand
          nl <- nl + 1
        }
      }
      
      # update S
      if (iter %% thin ==0) {
        Sups <- 0
        col.lam <- colSums(lam)
        Scand <-cbind(rnorm(N, S[,1], tune[3]), rnorm(N, S[,2], tune[3]))
        inbox <- inside.owin(Scand[,1], Scand[,2], wind)
        dtmp <- e2dist(Scand, X)
        lam.cand <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
        for(i in 1:N) {
          if(!inbox[i])
            next
          c.cand <- abs(col.lam-lam[i,]+lam.cand[i,])
          llcand <- sum(dpois(ndot, K*c.cand, log=TRUE) ) 
          if(runif(1)< exp(llcand - llcur)) {
            llcur <- llcand
            lam[i,] <- lam.cand[i,]
            col.lam <- c.cand
            S[i,] <- Scand[i,]
            D[i,] <- dtmp[i,]
            Sups <- Sups+1
          }
        }
      }
      
      #update N 
      eps <- sample(sx, 1)
      N.cand <- N + eps
      
      if ((N.cand-N) > 0) {
          S.cand <- cbind(runif((N.cand-N), xlim[1], xlim[2]),
                          runif((N.cand-N), ylim[1], ylim[2]))
          S.cand <-rbind(S, S.cand)
          D.cand <- e2dist(S.cand, X)
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE))
          
          if(runif(1) < exp( (llcand)  - (llcur)) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            S <- S.cand
            D <- D.cand
            nn <- nn + 1
          }
        } else {
          #omit the activity centres randomly (probabilistic)
          #since the order doesnt matter thus we can add them in such that the probability will cancel out
          #Hence, the acceptance probability stay the same
          ne <- sample(1:N, (N-N.cand), replace = FALSE)
          S.cand <- S[-ne,]
          D.cand <- D[-ne,]
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
          
          if(runif(1) < exp( (llcand)  - (llcur)) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            D <- D.cand
            S <- S.cand
            nn <- nn + 1
          }
      }
      setTxtProgressBar(pb,iter)    
      out[iter,] <- c(sigma, lam0, N)
      accept.rate[iter,] <- c(ns/iter, nl/iter, nn/iter)
    }
    close(pb)
    out.mcmc <- mcmc(out[-(1:burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}
```

## Box algorithm for non-square areas

``` r
brjmcmc <- function(n, X, K=NULL, niters=5000, xlim, ylim, tune=c(0.1, 0.1, 2, 5), prior=NULL,
                   monitorS=FALSE, chain=3, burn=5000, start.values=c(0.1,0.1,60), update=FALSE,
                   thin=1, p1=NULL, p2=NULL, wind=NULL, dS=NULL, box=NULL)
{
  J <- dim(n)[1]
  ndot <- rowSums(n)
  delta <- seq(1,tune[4],1)
  sx <- c(delta,-delta)
  result <- mcmc.list()
  niters <- niters + burn
  size <- box[1]*box[1]
  dA <- sum(dS)
  cx <- seq(xlim[1], xlim[2], length=box[1]+1)
  cy <- seq(ylim[1], ylim[2], length=box[1]+1)
  
  # initial values
  for (j in 1:chain) {
    sigma <- runif(1, 0.02, start.values[1])
    lam0 <- runif(1,0.1, start.values[2])
    N <- N.cand <- rpois(1,start.values[3])
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]), rep(0, N))
    inbox <- inside.owin(S[,1], S[,2], wind)
    #we need to check if all activity centres within the study area
    for (i in 1:N) {
      if(!inbox[i]){
        while (!inbox[i]) {
          S[i,1:2] <- c(runif(1, xlim[1], xlim[2]),
                        runif(1, ylim[1], ylim[2]))
          inbox[i] <- inside.owin(S[i,1], S[i,2], wind)
        }
      }
    }
    D <- e2dist(S[,1:2], X[,1:2])
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    llcur <- sum(dpois(ndot, K*colSums(lam), log=TRUE))
    
    # matrix to hold samples
    out <- matrix(0, nrow=niters, ncol=3)
    accept.rate <- matrix(NA, nrow=niters, ncol=3)
    colnames(out) <- c("sigma", "lam0", "N")
    colnames(accept.rate) <- c("sigma", "lam0", "N")
    check <- rep(0, niters)
    
    #cat("\n initial values =", c(sigma, lam0, N, "\n"))
    #cat("\n chain = ", j, "\n")
    
    close.time = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    ns <- nl <- nn <-0
    for(iter in 2:niters) {
      
      if(update == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("Accept S =", Sups/N, "\n")
        }
      }
      
      # update sigma
      
      sigma.cand <- sigma + rnorm(1,0,tune[1]) 
      
      if(is.null(prior))
      { 
        dSigma <- dunif(sigma, 0, 100)
        dSigma.cand <- dunif(sigma.cand, 0, 100)
      } else {
        dSigma <- dgamma(sigma, prior[1], prior[2])
        dSigma.cand <- dgamma(sigma.cand, prior[1], prior[2])
      }
      if (sigma.cand > 0){
        lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
        llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
        
        if(runif(1) < exp( llcand - llcur + log(dSigma.cand) - log(dSigma)) ){
          llcur <- llcand
          lam <- lam.cand
          sigma <- sigma.cand
          ns <- ns + 1
        }
      }
      
      
      # update lam0
      if (iter %% thin[2] ==0){
        lam0.cand <- lam0 + tune[2]*rnorm(1)
        if(lam0.cand>0) {
          lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
          
          if(runif(1) < exp( llcand  - llcur ) ){
            llcur <- llcand
            lam <- lam.cand
            lam0 <- lam0.cand
            nl <- nl + 1
          }
        }
      } 
      
      # update S
      if (iter %% thin[1] ==0) {
        Sups <- 0
        col.lam <- colSums(lam)
        Scand <-cbind(rnorm(N, S[,1], tune[3]), rnorm(N, S[,2], tune[3]))
        #we need to check if new activity centres within the study area
        inbox <- inside.owin(Scand[,1], Scand[,2], wind)
        dtmp <- e2dist(Scand, X[,1:2])
        lam.cand <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
        for(i in 1:N) {
          if(!inbox[i])
            next
          c.cand <- abs(col.lam-lam[i,]+lam.cand[i,])
          llcand <- sum(dpois(ndot, K*c.cand, log=TRUE) ) 
          if(runif(1)< exp(llcand - llcur)) {
            llcur <- llcand
            lam[i,] <- lam.cand[i,]
            col.lam <- c.cand
            S[i,1:2] <- Scand[i,]
            D[i,] <- dtmp[i,]
            Sups <- Sups+1
          }
        }
      }
      
      #update N
      eps <- sample(sx, 1)
      N.cand <- N + eps
      
      if (N.cand>1) {  
        if (eps > 0) {
          S.cand <- matrix(NA, nrow = eps, ncol=3)
          #we first generate box numbers according assigned probabilities p1
          id <- sample(1:size, eps, prob=p1, replace = TRUE)
          for (i in 1:eps) {
            #we determine the position of associated boxes to generate new activity centres
            ix <- ceiling(id[i]/box[1])
            iy <- id[i] - (ix-1)*box[1]
            S.cand[i,1:2] <- c(runif(1, cx[ix], cx[ix+1]),
                               runif(1, cy[iy], cy[iy+1]))
            inbox <- inside.owin(S.cand[i,1], S.cand[i,2], wind)
            #we need to check if new activity centres within the study area
            if(!inbox){
              while (!inbox) {
                S.cand[i,1:2] <- c(runif(1, cx[ix], cx[ix+1]),
                                   runif(1, cy[iy], cy[iy+1]))
                inbox <- inside.owin(S.cand[i,1], S.cand[i,2], wind)
              }
            }
          }
          #compute the associated probabilities of choosing boxes and generating activity centres
          pb <- p1[id] 
          dSA <- dS[id]
          D.cand <- e2dist(matrix(S.cand[,1:2],nrow=eps), X[,1:2])
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          lam.cand <- rbind(lam, lam.cand)
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE))
          priorS <- eps*log(1/dA)
          propS <- sum(log(pb)) + sum(log(1/dSA))

          
          if(runif(1) < exp( (llcand)  - (llcur) + priorS - propS ) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            S <- rbind(S,S.cand)
            D <- rbind(D,D.cand)
            nn <- nn + 1
          } 
        } else {
          pb <- numeric()
          dSA<-numeric()
  
          #omit acitivity centres randomly
          ne <- sample(1:N, (N-N.cand), replace = FALSE)
          S.cand <- S[-ne,]
          D.cand <- D[-ne,]
    
          ##compute the associated probabilities when reversing the move
          S[,3] <- boxcount(S[,1:2], xlim, ylim, box[1])
          pb <- p1[S[ne,3]]
          dSA <- dS[S[ne,3]]
          
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
          
          priorS <- abs(eps)*log(1/dA)
          propS <- sum(log(pb)) + sum(log(1/dSA))
      
          
          if(runif(1) < exp( (llcand)  - (llcur) + propS - priorS) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            D <- D.cand
            S <- S.cand
            nn <- nn + 1
          }
        }
      }
      
      setTxtProgressBar(close.time,iter)    
      out[iter,] <- c(sigma, lam0, N)
      accept.rate[iter,] <- c(ns/iter, nl/iter, nn/iter)
    }
    close(close.time)
    out.mcmc <- mcmc(out[-(1:burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}
```

## Box algorithm for square areas

``` r
brjmcmc <- function(n, X, K=NULL, niters=5000, xlim, ylim, tune=c(0.1, 0.1, 2, 5), prior=NULL,
                   monitorS=FALSE, chain=3, burn=5000, start.values=c(0.1,0.1,60), update=FALSE,
                   thin=1, p1=NULL, p2=NULL, wind=NULL, dS=NULL, box=NULL)
{
  J <- dim(n)[1]
  ndot <- rowSums(n)
  delta <- seq(1,tune[4],1)
  sx <- c(delta,-delta)
  result <- mcmc.list()
  niters <- niters + burn
  size <- box[1]*box[1]
  dA <- sum(dS)
  cx <- seq(xlim[1], xlim[2], length=box[1]+1)
  cy <- seq(ylim[1], ylim[2], length=box[1]+1)
  
  # initial values
  for (j in 1:chain) {
    sigma <- runif(1, 0.02, start.values[1])
    lam0 <- runif(1,0.1, start.values[2])
    N <- N.cand <- rpois(1,start.values[3])
    S <- cbind(runif(N, xlim[1], xlim[2]),
               runif(N, ylim[1], ylim[2]), rep(0, N))

    D <- e2dist(S[,1:2], X[,1:2])
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    llcur <- sum(dpois(ndot, K*colSums(lam), log=TRUE))
    
    # matrix to hold samples
    out <- matrix(0, nrow=niters, ncol=4)
    accept.rate <- matrix(NA, nrow=niters, ncol=3)
    colnames(out) <- c("sigma", "lam0", "N","logL")
    colnames(accept.rate) <- c("sigma", "lam0", "N")
    check <- rep(0, niters)
    
    #cat("\n initial values =", c(sigma, lam0, N, "\n"))
    #cat("\n chain = ", j, "\n")
    
    close.time = txtProgressBar(min = 1, max = niters, initial = 1, style = 3)
    ns <- nl <- nn <-0
    for(iter in 2:niters) {
      
      if(update == TRUE){
        if(iter %% 100 ==0) {
          cat("\n iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
          cat("current =", out[iter-1,], "\n")
          cat("Accept =", accept.rate[iter-1,], "\n")
          cat("Accept S =", Sups/N, "\n")
        }
      }
      
      # update sigma
      sigma.cand <- sigma + rnorm(1,0,tune[1]) 
      if(is.null(prior))
      { 
        dSigma <- dunif(sigma, 0, 100)
        dSigma.cand <- dunif(sigma.cand, 0, 100)
      } else {
        dSigma <- dgamma(sigma, prior[1], prior[2])
        dSigma.cand <- dgamma(sigma.cand, prior[1], prior[2])
      }
      if (sigma.cand > 0){
        lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
        llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
        
        if(runif(1) < exp( llcand - llcur + log(dSigma.cand) - log(dSigma)) ){
          llcur <- llcand
          lam <- lam.cand
          sigma <- sigma.cand
          ns <- ns + 1
        }
      }
      
      
      # update lam0
      if (iter %% thin[2] ==0){
        lam0.cand <- lam0 + tune[2]*rnorm(1)
        if(lam0.cand>0) {
          lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
          
          if(runif(1) < exp( llcand  - llcur ) ){
            llcur <- llcand
            lam <- lam.cand
            lam0 <- lam0.cand
            nl <- nl + 1
          }
        }
      } 
      
      # update S
      if (iter %% thin[1] ==0) {
        Sups <- 0
        col.lam <- colSums(lam)
        Scand <-cbind(rnorm(N, S[,1], tune[3]), rnorm(N, S[,2], tune[3]))
        inbox <- inside.owin(Scand[,1], Scand[,2], wind)
        dtmp <- e2dist(Scand, X[,1:2])
        lam.cand <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
        for(i in 1:N) {
          if(!inbox[i])
            next
          c.cand <- abs(col.lam-lam[i,]+lam.cand[i,])
          llcand <- sum(dpois(ndot, K*c.cand, log=TRUE) ) 
          if(runif(1)< exp(llcand - llcur)) {
            llcur <- llcand
            lam[i,] <- lam.cand[i,]
            col.lam <- c.cand
            S[i,1:2] <- Scand[i,]
            D[i,] <- dtmp[i,]
            Sups <- Sups+1
          }
        }
      }
      
      #update N
      eps <- sample(sx, 1)
      N.cand <- N + eps
      
      if (N.cand>1) {  
        if (eps > 0) {
          S.cand <- matrix(NA, nrow = eps, ncol=3)
          id <- sample(1:size, eps, prob=p1, replace = TRUE)
          for (i in 1:eps) {
            ix <- ceiling(id[i]/box[1])
            iy <- id[i] - (ix-1)*box[1]
            S.cand[i,1:2] <- c(runif(1, cx[ix], cx[ix+1]),
                               runif(1, cy[iy], cy[iy+1]))
          }
          pb <- p1[id] 
          dSA <- dS[id]
          D.cand <- e2dist(matrix(S.cand[,1:2],nrow=eps), X[,1:2])
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          lam.cand <- rbind(lam, lam.cand)
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE))
          priorS <- eps*log(1/dA)
          propS <- sum(log(pb)) + sum(log(1/dSA))
          
          if(runif(1) < exp( (llcand)  - (llcur) + priorS - propS) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            S <- rbind(S,S.cand)
            D <- rbind(D,D.cand)
            nn <- nn + 1
          } 
        } else {
          pb <- numeric()
          dSA<-numeric()
    
          ne <- sample(1:N, (N-N.cand), replace = FALSE)
          S.cand <- S[-ne,]
          D.cand <- D[-ne,]
        
          S[,3] <- boxcount(S[,1:2], xlim, ylim, box[1])
          pb <- p1[S[ne,3]]
          dSA <- dS[S[ne,3]]
          
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma))
          llcand <- sum(dpois(ndot, K*colSums(lam.cand), log=TRUE)) 
          
          priorS <- abs(eps)*log(1/dA)
          propS <- sum(log(pb)) + sum(log(1/dSA))
          
          if(runif(1) < exp( (llcand)  - (llcur) + propS - priorS) ){
            llcur <- llcand
            lam <- lam.cand
            N <- N.cand
            D <- D.cand
            S <- S.cand
            nn <- nn + 1
          }
        }
      }
      
      setTxtProgressBar(close.time,iter)    
      out[iter,] <- c(sigma, lam0, N, llcur)
      accept.rate[iter,] <- c(ns/iter, nl/iter, nn/iter)
    }
    close(close.time)
    out.mcmc <- mcmc(out[-(1:burn),])
    result[[j]] <- out.mcmc
  }
  return(result)
}
```
