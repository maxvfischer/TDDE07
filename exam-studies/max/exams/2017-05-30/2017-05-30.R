

################################### PROBLEM 1 ###################################
# Rice distribution: Rice(θ, ψ).
# θ: Location parameter, θ >= 0
# ψ: Related to the variance, ψ > 0
# I0(.): Modified Bassel function, first kind, order zero

# Function

rRice <-function(n = 1, theta = 1, psi = 1){
  x <- rnorm(n = n, mean = 0, sd = sqrt(psi))
  y <- rnorm(n = n, mean = theta, sd = sqrt(psi))
  return(sqrt(x^2+y^2))
}

logPostTheta <- function(X, ψ, θ) {
  n <- length(X)

  # Log likelihood
  bessel_func <- besselI(x = X*θ/ψ, nu = 0)
  log.likelihood <- log(X/ψ) - (X^2+θ^2)/(2*ψ) + log(bessel_func)
  
  # Log prior (Scaled inverse Chi squared)
  #log.prior <- dchisq(x = θ, df = (n-1), log = TRUE)
  log.prior <- 0

  return(sum(log.likelihood) + log.prior)
}

# Setup
ψ <- 1
riceData <- c(1.556, 1.861, 3.135, 1.311, 1.877, 0.622, 3.219, 0.768, 2.358, 2.056)
n <- length(riceData)
gridWidth = 0.01
thetas <- seq(0.01, 3, gridWidth)
post_thetas <- numeric()
for (t in 1:length(thetas)) {
  post_thetas[t] <- logPostTheta(X = riceData,
                              ψ = 1, 
                              θ = thetas[t]) 
}


theta_normalized <- exp(post_thetas)/((gridWidth)*sum(exp(post_thetas)))
plot(thetas, theta_normalized, type='l')

# b)
theta_init <- mean(thetas)
optim.res <- optim(par = theta_init,
                   fn = logPostTheta,
                   gr = NULL,
                   X = riceData,
                   ψ = 1,
                   control = list(fnscale = -1),
                   method=c("L-BFGS-B"),
                   lower = 0,
                   hessian = TRUE)

theta_mode <- optim.res$par
theta_sd <- sqrt(-1/optim.res$hessian)
nDraws <- length(theta_normalized)
norm_approx_thetas <- dnorm(x = thetas,
                            mean = theta_mode,
                            sd = theta_sd)

lines(thetas, norm_approx_thetas, col = "red")

# c)
# Compute the predictive distribution for a new observation x_t+1

# Draw from theta posterior
nDraws <- 10000
theta_draws <- rnorm(n = nDraws,
                     mean = theta_mode,
                     sd = theta_sd)
posterior_draws <- numeric()

for (t in 1:nDraws) {
  posterior_draws[t] <- rRice(n = 1, theta = theta_draws[t], psi = 1)
}

hist(posterior_draws, breaks = 20, probability = TRUE)
density_posterior_draws <- density(posterior_draws)
lines(density_posterior_draws, col = "red")

################################### PROBLEM 2 ###################################
# Modeling count data
load(file = 'bids.RData') # Data loaded into "bids"
bidsCounts <- table(bids) # Table with frequency of bids

# a)
# Likelihood: x_1, ..., x_n | θ ~ Poisson(θ)
# Prior: θ ~ Gamma(alpha, beta), alpha = beta = 1

alpha <- 1
beta <- 1
thetaGrid <- seq(3.4, 3.9, length = 1000)
n <- length(bids)
theta.posterior <- dgamma(x = thetaGrid,
                          shape = alpha + sum(bids),
                          rate = beta + n)
plot(thetaGrid, theta.posterior, type = 'l')

# b) Use graphical method to investigate if the Poisson models fits the data well

# i)
dataDistr <- bidsCounts/sum(bidsCounts)
dataGrid <- seq(min(bids), max(bids))
dataMean <- mean(bids)
poisDistr <- dpois(x = bids.grid,
                      lambda = bids.mean)
plot(dataGrid, poisDistr, type = 'l', col = "blue")
lines(dataGrid, dataDistr, type = 'l', col = "red")

# ii)
thetaMean <- mean(theta.posterior)
poisDistr <- dpois(x = dataGrid,
                   lambda = thetaMean)
plot(dataGrid, poisDistr, type = 'l', col = "blue")
lines(dataGrid, dataDistr)

# c) 

#################### SUPPLIED FUNCTION ####################
# Code for Problem 3 - Exam in Bayesian Learning 2017-05-30
GibbsMixPois <- function(x, nComp, alpha, alphaGamma, betaGamma, xGrid, nIter){
  
  # Gibbs sampling for a mixture of Poissons
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   x - vector with data observations (counts)
  #   nComp - Number of mixture components to be fitted
  #   alpha - The prior on the mixture component weights is w ~ Dir(alpha, alpha,..., alpha) 
  #   alphaGamma and betaGamma - 
  #              The prior on the mean (theta) of the Poisson mixture components is 
  #              theta ~ Gamma(alphaGamma, betaGamma) [rate parametrization of the Gamma dist]
  #   xGrid - the grid of data values over which the mixture is evaluated and plotted
  #   nIter - Number of Gibbs iterations
  #
  # OUTPUTS:
  #   results$wSample     - Gibbs sample of mixture component weights. nIter-by-nComp matrix
  #   results$thetaSample - Gibbs sample of mixture component means.   nIter-by-nComp matrix
  #   results$mixDensMean - Posterior mean of the estimated mixture density over xGrid.
  
  
  ####### Defining a function that simulates from a Dirichlet distribution
  rDirichlet <- function(param){
    nCat <- length(param)
    thetaDraws <- matrix(NA,nCat,1)
    for (j in 1:nCat){
      thetaDraws[j] <- rgamma(1,param[j],1)
    }
    thetaDraws = thetaDraws/sum(thetaDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
    return(thetaDraws)
  }
  
  # Simple function that converts between two different representations of the mixture allocation
  S2alloc <- function(S){
    n <- dim(S)[1]
    alloc <- rep(0,n)
    for (i in 1:n){
      alloc[i] <- which(S[i,] == 1)
    }
    return(alloc)
  }
  
  # Initial values for the Gibbs sampling
  nObs <- length(x)
  S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
  theta <- rep(mean(x), nComp) # Each component is initialized at the mean of the data
  
  # Setting up the grid where the mixture density is evaluated.
  mixDensMean <- rep(0,length(xGrid))
  effIterCount <- 0
  
  # Setting up matrices to store the draws
  wSample <- matrix(0, nIter, nComp)
  thetaSample <- matrix(0, nIter, nComp)
  probObsInComp <- rep(NA, nComp)
  
  # Setting up the priors - the same prior for all components
  alpha <- rep(alpha, nComp) 
  alphaGamma <- rep(alphaGamma, nComp) 
  betaGamma <- rep(betaGamma, nComp) 
  
  # HERE STARTS THE ACTUAL GIBBS SAMPLING
  
  for (k in 1:nIter){
    message(paste('Iteration number:',k))
    alloc <- S2alloc(S) # Function that converts between different representations of the group allocations
    nAlloc <- colSums(S)
    
    # Step 1 - Update components probabilities
    w <- rDirichlet(alpha + nAlloc)
    wSample[k,] <- w
    
    # Step 2 - Update theta's in Poisson components
    for (j in 1:nComp){
      theta[j] <- rgamma(1, shape = alphaGamma + sum(x[alloc == j]), rate = betaGamma + nAlloc[j])
    }
    thetaSample[k,] <- theta
    
    # Step 3 - Update allocation
    for (i in 1:nObs){
      for (j in 1:nComp){
        probObsInComp[j] <- w[j]*dpois(x[i], lambda = theta[j])
      }
      S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
    }
    
    # Computing the mixture density at the current parameters, and averaging that over draws.
    effIterCount <- effIterCount + 1
    mixDens <- rep(0,length(xGrid))
    for (j in 1:nComp){
      compDens <- dpois(xGrid, lambda = theta[j])
      mixDens <- mixDens + w[j]*compDens
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
  }
  return(results = list(wSample = wSample, thetaSample = thetaSample, mixDensMean = mixDensMean))
}

# GibbsMixPois <- function(x, nComp, alpha, alphaGamma, betaGamma, xGrid, nIter)
# Theta_k: Mean of Pois(x|theta_k)
# Prior theta ~ Gamma(1, 1) for all k
# Prior weights (w) ~ Uniform(0, 1)

# Estimate mixture of poisson both with K = 2 and K = 3

# Setup
xGrid <- seq(min(bids), max(bids))
nIter <- 500
alphaGamma <- 1
betaGamma <- 1
alpha <- 1
nComp1 <- 2
nComp2 <- 3

gibbs_sampling1 <- GibbsMixPois(x = bids,
                                nComp = nComp1,
                                alpha = alpha,
                                alphaGamma = alphaGamma,
                                betaGamma = betaGamma,
                                xGrid = xGrid,
                                nIter = nIter)

gibbs_sampling2 <- GibbsMixPois(x = bids,
                                nComp = nComp2,
                                alpha = alpha,
                                alphaGamma = alphaGamma,
                                betaGamma = betaGamma,
                                xGrid = xGrid,
                                nIter = nIter)

# d) Use graphical method to investigate if mixture of Poisson with K = 2 fits the data well

dataDistr <- bidsCounts/sum(bidsCounts)


plot(xGrid, gibbs_sampling1$mixDensMean, type="l", col = "red")
lines(xGrid, gibbs_sampling2$mixDensMean, col = "blue")
lines(xGrid, dataDistr)

# e) 
# A true natural born bayesian would do bayesian model comparison over K (which are different models). 
# Pr(K = k | y) = Pr(y | K = k) * Pr(K = k)
# where
# Pr(y | K = k) = integral[ p(y | theta_k, K = k) * p(theta_k | K = k) dtheta]
# is the marginal likelihood of the model with k mixture components and theta_k cointains all parameters in this model.
# p(y | theta_k, K = k) is the usual likelihood for mixture model with k components
# p(theta_k | K = k) is the prior of theta in model with k components
# This integral is hard to compute, but can be solved with MCMC.

################################### PROBLEM 3 ###################################
# Regression

#################### SUPPLIED FUNCTION ####################
# Defining a function that simulates from the scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){
  # Direct sampling from a Gaussian linear regression with conjugate prior:
  #
  # beta | sigma2 ~ N(mu_0, sigma2*inv(Omega_0))
  # sigma2 ~ Inv-Chi2(v_0,sigma2_0)
  # 
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   Omega_0  - prior precision matrix for beta
  #   v_0      - degrees of freedom in the prior for sigma2
  #   sigma2_0 - location ("best guess") in the prior for sigma2
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   results$betaSample     - Posterior sample of beta.     nIter-by-nCovs matrix
  #   results$sigma2Sample   - Posterior sample of sigma2.   nIter-by-1 vector
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  betaHat <- solve(XX,t(X)%*%y)
  Omega_n = XX + Omega_0
  mu_n = solve(Omega_n,XX%*%betaHat+Omega_0%*%mu_0)
  v_n = v_0 + n
  sigma2_n = as.numeric((v_0*sigma2_0 + ( t(y)%*%y + t(mu_0)%*%Omega_0%*%mu_0 - t(mu_n)%*%Omega_n%*%mu_n))/v_n)
  invOmega_n = solve(Omega_n)
  
  # The actual sampling
  sigma2Sample = rep(NA, nIter)
  betaSample = matrix(NA, nIter, nCovs)
  for (i in 1:nIter){
    
    # Simulate from p(sigma2 | y, X)
    sigma2 = rScaledInvChi2(n=1, df = v_n, scale = sigma2_n)
    sigma2Sample[i] = sigma2
    
    # Simulate from p(beta | sigma2, y, X)
    beta_ = rmvnorm(n=1, mean = mu_n, sigma = sigma2*invOmega_n)
    betaSample[i,] = beta_
    
  }
  return(results = list(sigma2Sample = sigma2Sample, betaSample=betaSample))
}

# beta | sigma2 ~ N(mu0, sigma2 * Zigma-1)
# sigma2 ~ Inv-X(v0, sigma0_2)

# mgp: Miles per gallon
# Weight: 
# sixcyl/eightcyl:
# If four: (0, 0)
# If six: (1, 0)
# If eight: (0, 1)

# BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter)
load("cars.RData")

library(mvtnorm)

# Setup
mu_0 <- c(0, 0, 0, 0)
Omega_0 <- 0.01*diag(4)
v_0 <- 1
sigma2_0 <- 36
nIter <- 1000

linearReg.res <- BayesLinReg(y = cars$mpg, 
                             X = as.matrix(cars[, 2:5]), 
                             mu_0 = mu_0, 
                             Omega_0 = Omega_0,
                             v_0 = v_0,
                             sigma2_0 = sigma2_0,
                             nIter = nIter)

beta.draws <- linearReg.res$betaSample
sigma2.draws <- linearReg.res$sigma2Sample

# i)
par(mfrow = c(2, 3))
hist(beta.draws[,1], breaks = 20, xlab = "Beta 1")
hist(beta.draws[,2], breaks = 20, xlab = "Beta 2")
hist(beta.draws[,3], breaks = 20, xlab = "Beta 3")
hist(beta.draws[,4], breaks = 20, xlab = "Beta 4")
hist(sigma2.draws, breaks = 20, xlab = "Sigma2")

# ii) 
# As it's a linear loss function, the most appropriate point estimate is the median
beta.draws.median <- apply(beta.draws, MARGIN = 2, median)

# iii)
beta.CI_intervals <- t(apply(beta.draws, MARGIN = 2, function(x) quantile(x, probs = c(0.025, 0.975))))
sigma2.CI_interval <- t(apply(X = as.matrix(sigma2.draws), 
                            MARGIN = 2, 
                            FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
CI_intervals <- rbind(beta.CI_intervals, sigma2.CI_interval)
rownames(CI_intervals) <- c("β1", "β2", "β3", "β4", "σ2")

par(mfrow = c(2, 3))
hist(beta.draws[,1], breaks = 20, xlab = "Beta 1")
abline(v = CI_intervals[1,], col = "red")
hist(beta.draws[,2], breaks = 20, xlab = "Beta 2")
abline(v = CI_intervals[2,], col = "red")
hist(beta.draws[,3], breaks = 20, xlab = "Beta 3")
abline(v = CI_intervals[3,], col = "red")
hist(beta.draws[,4], breaks = 20, xlab = "Beta 4")
abline(v = CI_intervals[4,], col = "red")
hist(sigma2.draws, breaks = 20, xlab = "Sigma2")
abline(v = CI_intervals[5,], col = "red")

# b)
#  Investigate if the effect on mpg is different in cars with six cylinders compared to 
# cars with 8 cylinder

# If we look at the medians of the coefficients of sixcyl and eightcyl, it seems like
# eightcyl with the median coefficient of -6.213359 affects mpg in a more decreasing way than
# sixcyl with the median coefficient of -4.323130.

# c) 
# mpg ~ N[ β0 + β1*weight + β2*sixcyl + β3*eightcyl, sigma2]
# rScaledInvChi2 <- function(n, df, scale){
# mu_0 <- c(0, 0, 0, 0)
# Omega_0 <- 0.01*diag(4)
# v_0 <- 1
# sigma2_0 <- 36
# nIter <- 1000

target <- c(1, 3.5, 0, 0)

predictive.mean <- beta.draws %*% target

mpg.simulation <- beta.draws %*% target + rnorm(n = nIter, mean = 0, sd = sqrt(sigma2.draws))

hist(mpg.simulation, breaks = 50)

################################### PROBLEM 4 ###################################


# c)
# W, L, L, W, W, L, L, L, W, 
n <- 10
