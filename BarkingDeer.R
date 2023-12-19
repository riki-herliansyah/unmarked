library(MASS)
library(coda)
library(readxl)
library(ggplot2)
library(maptools)
library(spatstat)

source('R/C-template.R')
source('R/RJMCMC.R')
source('R/SPA.R')

# import the data set (can requested to the author)
muntjak <- read_excel("muntjak.xlsx")
nj <- as.matrix(muntjak[,4:5])

# import the study area
patch <- read.csv("studyarea.csv")
# set as a window object
win <- as.owin(patch)
# import the trap locations
traps <- as.matrix(muntjak[,2:3])
xlim <- c(-5, 265)
ylim <- c(-5, 245)

# create initial values for parameter
inits = list(sigma = 1,
            lam0 = 1.5,
            N = 3500)
# create tuning parameters
tune.par = list( sigma = 0.05,
                 lam0 = 0.12,
                 S = 65,
                 N = 120)
# set the data set assuming K = 9
dataIn <- list(n = nj[,2],
              K = 9,
              X = traps,
              xlim = xlim,
              ylim = ylim,
              wind = win)
# assume a uniform prior on sigma, a normal prior on log(lambda_0)
# and NB ~ (10, 0.0032) for N
prior <- list(sigma = NULL,
              lam0 = c(0, 10),
              N = c(10, 0.0032))

tstart = proc.time()
model <- rjmcmcIR(data = dataIn, tune = tune.par, prior = prior, control = list(n.iters=1000, 
                  n.chain=1, n.burn=1000, n.update=15, monitor = TRUE, tol=1e-10), 
                  inits = inits)
time <- proc.time() - tstart

summary(model)