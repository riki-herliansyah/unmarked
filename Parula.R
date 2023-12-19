library(coda)
library(Rcpp)
source('R/C-template.R')
source('R/RJMCMC.R')
source('R/SPA.R')
# Parula data modelling
# input the data provided in Candler and Royle (2013)
nopaDat <-
  structure(list(y = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
                                 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0,
                                 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1,
                                 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0,
                                 0, 1, 2, 1, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
                                 1, 0, 2, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
                                 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 1, 1, 1, 0, 0, 0,
                                 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                                 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                                 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0), .Dim = c(105L, 3L), .Dimnames = list(
                                   c("4005", "3055", "3005", "2055", "2005", "1055", "1005",
                                     "4100", "3150", "3100", "2150", "2100", "1150", "1100", "4105",
                                     "3155", "3105", "2155", "2105", "1155", "1105", "4200", "3250",
                                     "3200", "2250", "2200", "1250", "1200", "4205", "3255", "3205",
                                     "2255", "2205", "1255", "1205", "4300", "3350", "3300", "2350",
                                     "2300", "1350", "1300", "4305", "3355", "3305", "2355", "2305",
                                     "1355", "1305", "4400", "3450", "3400", "2450", "2400", "1450",
                                     "1400", "4405", "3455", "3405", "2455", "2405", "1455", "1405",
                                     "4500", "3550", "3500", "2550", "2500", "1550", "1500", "4505",
                                     "3555", "3505", "2555", "2505", "1555", "1505", "4600", "3650",
                                     "3600", "2650", "2600", "1650", "1600", "4605", "3655", "3605",
                                     "2655", "2605", "1655", "1605", "4700", "3750", "3700", "2750",
                                     "2700", "1750", "1700", "4705", "3755", "3705", "2755", "2705",
                                     "1755", "1705"), c("v1", "v2", "v3"))), X = structure(c(5,
                                                                                             5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8,
                                                                                             8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10,
                                                                                             10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 13,
                                                                                             13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15,
                                                                                             15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17,
                                                                                             17, 17, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19,
                                                                                             5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9,
                                                                                             10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7,
                                                                                             8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5,
                                                                                             6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10,
                                                                                             11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8,
                                                                                             9, 10, 11, 5, 6, 7, 8, 9, 10, 11), .Dim = c(105L, 2L), .Dimnames = list(
                                                                                               NULL, c("x", "y"))), xSide = 24, ySide = 16, M = 100, nTraps = 105L,
                 nReps = 3L), .Names = c("y", "X", "xSide", "ySide", "M",
                                         "nTraps", "nReps"))


# we use buffer equal to 3 to create the state space
# the proposal distribution of sigma and lam0 is with 0.3 and 0.15 being tuning parameters respectively
# we use Gaussian random walk for activity centers with respective variance equal to 2.25
#n For updating N, we set epsilon equal to 10 for Parula dataset.
# We generate 300,000 samples and discard the first 10,000
# the number of chain is 3
buffer = 3

xlim <- c(2, nopaDat$xSide + 3)
ylim <- c(2, nopaDat$ySide + 3)

# create initial values for mcmc
inits = list(sigma = runif(1,0.5,1.5),
             lam0 = runif(1,0.1,1),
             N = 50,
             z = rbinom(300, 1, .7),
             psi = runif(1, .2, .8))
# create data set for the model
# assumming sampling occassions K = 3 and the upper limit for the population M = 300
data <- list( n = rowSums(nopaDat$y),
              K = 3,
              X = nopaDat$X,
              xlim = xlim,
              ylim = ylim,
              M = 300)
#assuming uniform prior on sigma and N
prior <- list(sigma = NULL,
              N = NULL)
# tuning parameters for sigma, lambda_0, N and activity centres S
tune = list( sigma = 0.1,
             lam0 = 0.1,
             S = 1,
             N = 10)
# Fitting the model using SPA
set.seed(171)
tstart = proc.time()
spa.parula <- spaR(data = data, tune = tune, inits = inits, prior = prior, 
                    control = list(n.iters = 30000, n.chain = 3, n.burn = 10000, 
                                   monitor = TRUE, tol = 1e-15))
time.spa2 <- proc.time() - tstart

# Fitting the model using RJMCMC for a stochastic removal (stochastic = TRUE)
# Set stochastic = FALSE for fixed RJMCMC
set.seed(171)
tstart = proc.time()
rj.parula <- rjmcmcR(data = data, tune = tune, inits = inits, prior = prior, stochastic = TRUE,
                       control = list(n.iters = 300000, n.chain = 3, n.burn = 10000, monitor = TRUE, 
                                      n.update = 4, tol = 1e-15))
time.rjs2 <- proc.time() - tstart

summary(spa.parula)
plot(spa.parula)