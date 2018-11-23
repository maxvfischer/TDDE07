######################## QUESTION 1 ########################

load(file = 'yProportions.RData')
n <- length(yProp)

# a) 
thetaGrid <- seq(0.01, 15, length = 1000)
theta.posterior <- numeric()

for (t in 1:length(thetaGrid)) {
  
  # log.likelihood
  log.likelihood <- sum(dbeta(x = yProp, shape1 = thetaGrid[t], shape2 = thetaGrid[t], log = TRUE))
  
  # log.prior
  log.prior <- dexp(x = thetaGrid[t], rate = 1, log = TRUE)
  
  theta.posterior[t] <- log.likelihood + log.prior
}

plot(thetaGrid, exp(theta.posterior), type = 'l')

# b) 
# Function
logPost <- function(thetas, Y) {
  
  # log.likelihood
  log.likelihood <- sum(dbeta(x = Y, shape1 = thetas[1], shape2 = thetas[2], log = TRUE))
  
  # log.priors
  log.prior1 <- dexp(x = thetaGrid[1], rate = 1, log = TRUE)
  log.prior2 <- dexp(x = thetaGrid[2], rate = 1, log = TRUE)
  
  return(log.likelihood + log.prior1 + log.prior2)
}

# Setup
thetas.init <- c(1, 1)

optim.res <- optim(par = thetas.init,
                   fn = logPost,
                   gr = NULL,
                   control = list(fnscale = -1),
                   method=c("L-BFGS-B"),
                   lower = c(0.000001, 0.000001),
                   hessian = TRUE,
                   Y = yProp)

thetas.mean <- optim.res$par
thetas.covariate <- -solve(optim.res$hessian)

# c) 
# A natural born bayesian would calculate the marginal likelihood of the data for both models, then
# calculate the Bayes factor.

######################## QUESTION 2 ########################
# Regression

########## Given code ##########

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

library(MASS)
library(mvtnorm)
BostonHousing = Boston
y = BostonHousing$medv
X = cbind(1,BostonHousing[,1:13]) # Adding a column of ones for the intercept
names(X)[1] <- "intercept"
covNames <- names(X)
y <- as.numeric(y)
X <- as.matrix(X)
nCovs <- dim(X)[2]
# a)
# (y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){

# Prior betas
mu_0 <- rep(0, nCovs)
Omega_0 <- (1/(10^2))*diag(nCovs)

# Prior sigma2
sigma2_0 <- 6^2
v_0 <- 1

nIter <- 5000

# 5000 simulations of betas and sigma2
posterior.simulation <- BayesLinReg(y = y,
                                    X = X,
                                    mu_0 = mu_0,
                                    Omega_0 = Omega_0,
                                    v_0 = v_0,
                                    sigma2_0 = sigma2_0,
                                    nIter = nIter)
colnames(posterior.simulation$betaSample) <- covNames

lstat.beta.sim <- posterior.simulation$betaSample[,"lstat"]
lstat.density <- density(lstat.beta.sim)
lstat <- data.frame(x = lstat.density$x, 
                    y = lstat.density$y)

lstat.sorted <- lstat[order(lstat$y),]


