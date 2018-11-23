################################### PROBLEM 1 ###################################
# 

load(file = 'CauchyData.RData')

dCauchy <- function(x, theta = 0, gamma = 1){
  return(dens = (1/(pi*gamma))*(1/(1+((x-theta)/gamma)^2)))
}

# Setup
gamma <- 1
mu <- 0
sigma2 <- 10^2

gridWidth <- 0.01
thetaGrid <- seq(2, 6.5, gridWidth)
log.posterior.theta <- numeric()
  
for (t in 1:length(thetaGrid)) {
  
  # Log likelihood
  log.likelihood <- log(dCauchy(yVect,
                                theta = thetaGrid[t],
                                gamma = gamma))
  
  # Log prior
  log.prior <- dnorm(x = thetaGrid[t],
                     mean = mu,
                     sd = sqrt(sigma2),
                     log = TRUE)
  
  log.posterior.theta[t] <- sum(log.likelihood + log.prior)
  
}
posterior.theta <- exp(log.posterior.theta)

plot(thetaGrid, posterior.theta/(sum(posterior.theta) * gridWidth), type= 'l')

# b)
# Functions
dlognormal <- function(x, mu, sigma2){
  return(dens = (1/(sqrt(2*pi*sigma2)*x))*exp((-1/(2*sigma2))*(log(x)-mu)^2))
}

logPostCauchy <- function(θγ, Y) {
  θ <- θγ[1]
  γ <- θγ[2]
  
  # Log likelihood 
  log.likelihood <- log(dCauchy(Y,
                            theta = θ,
                            gamma = γ))
  
  # Theta log prior
  theta.log.prior <- dnorm(x = θ,
                           mean = 0,
                           sd = sqrt(10^2),
                           log = TRUE)
  
  # Gamma log prior
  gamma.log.prior <- log(dlognormal(x = γ,
                                mu = 0,
                                sigma2 = 1))
  
  return (sum(log.likelihood + theta.log.prior + gamma.log.prior))
}

# Also gamma is unknown
# Likelihood: p(y | γ, θ) ~ Cauchy()
# Prior: p(θ) ~ Normal(0, 10^2)
# Prior: p(γ) ~ lognormal(0, 1)

# Setup
θγ.init <- c(1, 1)

optim.res <- optim(par = θγ.init,
                   fn = logPostCauchy,
                   Y = yVect,
                   gr = NULL,
                   lower = c(-Inf, 0.000001),
                   upper = c(Inf, Inf),
                   method = c("L-BFGS-B"),
                   control = list(fnscale = -1),
                   hessian = TRUE)

hessian <- optim.res$hessian
mean <- optim.res$par

# c)
nDraws <- 10000

posterior.draws <- rmvnorm(n = nDraws,
                           mean = mean,
                           sigma = -solve(hessian))
margPostCauchy <- posterior.draws[,1] + posterior.draws[,2]*tan(pi * (0.99 - 0.5))
margPostCauchy.density <- density(margPostCauchy)
hist(margPostCauchy, breaks = 20, freq = FALSE)
lines(margPostCauchy.density, col = "red")

################################### PROBLEM 2 ###################################
# Response variable: medv (Median value of the house in 1000$)

library(mvtnorm)

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

# BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){
library(MASS)

# Setup
BostonHousing = Boston
y = BostonHousing$medv
X = cbind(1,BostonHousing[,1:13]) # Adding a column of ones for the intercept
names(X)[1] <- "intercept"
covNames <- names(X)
y <- as.numeric(y)
X <- as.matrix(X)
nIter <- 5000

# Prior beta
nCovs <- dim(X)[2]
mu_0 <- rep(0, nCovs)
Omega_0 <- (1/(10^2))*diag(nCovs)

# Prior sigma2
v_0 <- 1
sigma2_0 <- 6^2




linRegSim <- BayesLinReg(y = y,
                         X = X,
                         mu_0 = mu_0,
                         Omega_0 = Omega_0,
                         v_0 = v_0,
                         sigma2_0 = sigma2_0,
                         nIter = nIter)

# Point estimate for quadratic loss function
posterior.mean.betas <- t(apply(linRegSim$betaSample, MARGIN = 2, mean))
posterior.mean.sigma2 <- mean(linRegSim$sigma2Sample)
posterior.mean <- cbind(posterior.mean.betas, posterior.mean.sigma2)

# 95% equal tail credible interval
posterior.CI_interval.betas <- apply(linRegSim$betaSample, MARGIN = 2, function(x) quantile(x, c(0.025, 0.975)))
posterior.CI_interval.sigma2 <- apply(as.matrix(linRegSim$sigma2Sample), MARGIN = 2, function(x) quantile(x, c(0.025, 0.975)))

# Combine Credible intervals of beta and sigma2
posterior.CI_intervals <- cbind(posterior.CI_interval.betas, posterior.CI_interval.sigma2)
colnames(posterior.CI_intervals) <- c("β1", "β2", "β3", "β4", "β5", "β6", "β7", "β8", "β9", "β10", "β11", "β12", "β13", "β14", "σ2")

# Analyse CI for rm
hist(linRegSim$betaSample[, 7], breaks = 50)
abline(v = posterior.CI_intervals[, "β7"])

# b)

# Setup
x381 <- X[381,]
y381 <- y[381]
x381["crim"] <- 10
n <- length(linRegSim$sigma2Sample)

y_pred <- linRegSim$betaSample %*% x381 + rnorm(n = n, mean = 0, sd = sqrt(linRegSim$sigma2Sample))

y_pred_CI_99 <- quantile(y_pred, c(0.005, 0.995))
y_pred_mean <- mean(y_pred)
hist(y_pred, breaks = 100)
abline(v = y_pred_mean, col = "blue")
abline(v = y_pred_CI_99[1], col = "red")
abline(v = y_pred_CI_99[2], col = "red")

prob_above_20 <- sum(y_pred >= 20)/n
prob_above_30 <- sum(y_pred >= 30)/n

# c)
# 

################################### PROBLEM 3 ###################################
# On paper

################################### PROBLEM 4 ###################################
# Firm produces products
# X_t denotes quality demanded of the product in quarter t
# X_t | θ (iid) ~ Poisson(θ)

# Demand previous year:
X <- c(220, 323, 174, 229)

# a)
# Likelihood: X_t | θ (iid) ~ Poisson(θ)
# Conjugate prior: θ ~ Gamma(α, β)
# Mean = 250, Sd = 50
# α = 25, β = 1/10
 
# Setup
nDraws = 1000
n <- length(X)

# Prior
prior.alpha <- 25
prior.beta <- 1/10

# Posterior
posterior.simulations <- rgamma(n = nDraws,
                                shape = sum(X) + prior.alpha,
                                rate = n + prior.beta)

hist(posterior.simulations)

# b)
# Setup
pred.draws <- numeric()
for (i in 1:length(posterior.simulations)) {
  pred.draws[i] <- rpois(n = 1,
                         lambda = posterior.simulations[i])
}
hist(pred.draws, freq = FALSE, breaks = 30) 
sum(pred.draws <= 200)/length(pred.draws)


## Real distribution
gridWidth <- 1
xGrid <- seq(100, 400, gridWidth)

nDraws <- 1000
pred.posterior <- numeric()
for (t in 1:length(xGrid)) {
  pred.posterior[t] <- (sum(X) + prior.alpha)*log(n + prior.beta) + lgamma(xGrid[t] + sum(X) + prior.alpha)
  pred.posterior[t] <- pred.posterior[t] - lgamma(sum(X) + prior.alpha) - lfactorial(xGrid[t]) - (xGrid[t] + sum(X) + prior.alpha)*log(1 + n + prior.beta)
}
  
lines(xGrid, exp(pred.posterior), col = "red")
sum(exp(pred.posterior[which(xGrid <= 300)]))
mean = xGrid[which.max(exp(pred.posterior))]
sd <- xGrid[which(cumsum(exp(pred.posterior) = 3/18))]
##