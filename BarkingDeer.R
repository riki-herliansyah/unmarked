library(MASS)
library(coda)
library(readxl)
library(ggplot2)
library(maptools)
library(spatstat)

source('R/C-template.R')
source('R/RJMCMC.R')
source('SPA.R')

# import the data set (can requested to the author)
muntjak <- read.delim("muntjak.txt")
nj <- as.matrix(muntjak)

# import the study area
patch <- read.csv("studyarea2.csv")
# set as a window object
win <- as.owin(patch)
# import the trap locations
traps <- as.matrix(read.csv("traps.csv"))
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
# set the data set; the first two column for wet season and the last two are for dry season
# assuming K = 9
data <- list( n = rowSums(nj[, 3:4]),
              K = 9,
              X = traps,
              xlim = xlim,
              ylim = ylim,
              wind = win)
# assume a uniform prior on sigma, a normal prior on log(lambda_0)
# and NB ~ (10, 0.0032) for N
prior <- list( sigma = NULL,
               lam0 = c(0, 10),
               N = c(10, 0.0032))

tstart = proc.time()
model <- rjmcmcIR(data = data, tune = tune.par, prior = prior, control = list(n.iters=500000, 
                  n.chain=1, n.burn=10000, n.update=15, monitor = FALSE, tol=1e-10), 
                  inits = inits)
time <- proc.time() - tstart

summary(model)