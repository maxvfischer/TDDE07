########################## Lab1 ##########################
### Meta-info
## Beta distribution
## Simulations
## Gini coefficient
## Credible interval
## Highest Posterior Density HPD

# Task 1: Bernoulli ... again
# a) Draw random numbers from Beta distribution and graphical verification
# of posterior
# b) Simulate to compute posterior prob of Pr(theta < 0.4)
# c) Compute log-oods posterior distribution

# Task 2: Log-normal distribution and the Gini coefficient
# a) Simulate 1000 draws from posterior or theta2. Compare with real value
# b) Compute posterior distribution of Gini coefficient G
# c) Compute 95% equal tail credible interval of Gini coefficient G.
# Doing a kernal density estimate
# Compute 95% Highest Posterior Density interval (HPD) of G

# Task 3: Bayesian inference for the concentration parameter in the von Mises distributio
# a) Plot posterior distribution of kappa for wind direction data
# b) Find approximate posterior mode of kappa

#####################################################################################################

############# Task 1 #############
# Bernoulli ... again
#a)
# Instrucitons
# Likelihood: y_1, ..., y_n | θ ~ Bern(θ)
# Prior: θ ~ Beta(alpha_0, beta_0), alpha_0 = beta_0 = 2
# Posterior: θ|y_1, ..., y_n ~ Beta(alpha_0 + s, beta_0 + f)
# s = 14
# n = 20
# f = 6

# Setup
n = 20
s = 14
f = n - s
nDraws = 50000
drawsInterval = 10
intervalVec <- seq(10, nDraws, drawsInterval)

# Prior
prior.alpha = 2
prior.beta = 2

# Posterior
posterior.alpha <- prior.alpha + s
posterior.beta <- prior.beta + f
posterior.draws_means = numeric()
posterior.draws_sd = numeric()

# For-loop - (10, 20, 30, ..., 49980, 49990)
for (i in intervalVec) {
  posterior.draws <- rbeta(n = i, 
                           shape1 = posterior.alpha, 
                           shape2 = posterior.beta) # Draw from beta
  
  posterior.draws_means <- c(posterior.draws_means, mean(posterior.draws)) # Add mean to vector of means
  posterior.draws_sd <- c(posterior.draws_sd, sd(posterior.draws)) # Add sd to vector of sd
}

# True values
posterior.true_mean <- posterior.alpha/(posterior.alpha + posterior.beta)
posterior.true_sd <- sqrt((posterior.alpha*posterior.beta)/
                            ((posterior.alpha + posterior.beta)^2 * (posterior.alpha + posterior.beta + 1)))

# Plot
par(mfrow = c(1, 2))
plot(x = intervalVec, 
     y = posterior.draws_means, 
     type = 'l',
     xlab = 'No. of draws',
     ylab = 'Mean of data') # Plot means
abline(h = posterior.true_mean, col = 'red') # Add line of real mean to plot

plot(x = intervalVec, 
     y = posterior.draws_sd, 
     type = 'l',
     xlab = 'No. of draws',
     ylab = 'Sd of data') # Plot sd's
abline(h = posterior.true_sd, col = 'red')

# b) 
# Setup
n = 20
s = 14
f = n - s
nDraws = 10000

# Prior
prior.alpha = 2
prior.beta = 2

# Posterior
posterior.alpha <- prior.alpha + s
posterior.beta <- prior.beta + f
posterior.draws_means = numeric()
posterior.draws_sd = numeric()

# Draws
posterior.draws_10000 <- rbeta(n = nDraws,
                               shape1 = posterior.alpha,
                               shape2 = posterior.beta)

# Calculate probability
posterior.prob_0.4 <- length(which(posterior.draws_10000 < 0.4))/length(posterior.draws_10000)
posterior.prob_real<- pbeta(q = 0.4,
                            shape1 = posterior.alpha,
                            shape2 = posterior.beta)

# c)
# Functions
logOdds <- function(theta) {
  return (log(theta/(1-theta)))
}

# Setup
n = 20
s = 14
f = n - s
nDraws = 10000

# Prior
prior.alpha = 2
prior.beta = 2

# Posterior
posterior.alpha <- prior.alpha + s
posterior.beta <- prior.beta + f

# Draws
posterior.draws_10000 <- rbeta(n = nDraws,
                               shape1 = posterior.alpha,
                               shape2 = posterior.beta)

# Log-odds the draws
logOdds.draws <-logOdds(posterior.draws_10000)

# Hist and plot the density function of the log-odds draws
hist(logOdds.draws, probability = TRUE)
lines(density(logOdds.draws), col = 'red')

############# Task 2 #############
# Log-normal distribution and the Gini coefficient
# Likelihood: y_1, ..., y_n | µ, σ2 ~ log[ N(µ, σ2) ], µ known, σ2 unknown
# Prior: p(σ2) ∝ 1/σ2

# Posterior of σ2: Inv-X(n, tao^2)
# Tao^2 - The sample variance. Calculated as following:
# sum[ (log(y_i) - µ)^2 ]/n

# If X ~ N(0, 1)
# Y = exp(X) ~ log[ N(0,1) ] (lognormal)

# Setup
y <- c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
nDraws = 10000
mu = 3.5

# a)
# Functions
scaled_rchisq <- function(Y, mu, nDraws) {
  n <- length(Y) # Length of data
  Z <- rchisq(n = nDraws, df = n) # Draw from Chi-squared distribution
  tao2 <- sum((log(Y) - mu)^2)/n # Calculate tau2 (sample standard diviation)
  
  return (n*tao2/Z)
}

σ2.draws <- scaled_rchisq(y, mu, nDraws) # Draw from Scaled inverse chi-squared distribution
hist(σ2.draws, breaks=1000)

# b)
# Gini-coefficient measuers inequality (0: equal, 1: unequal)
# Uses the Lorenz curve where:
# x-axis: Cumulative share of people from low to high income
# y-axis: Cumulative share of income earned
#
# If a straight line, it's 100% equal
#
# The Gini-coefficient is the ratio between the area between the straight line and the Lorenz curve
# divided by the total area
#
# If the data follows an lognormal distribution (e.g wealth), the Gini-coef is calculated as follows:
# G = 2 * Φ(σ/√2) - 1

σ <- sqrt(σ2.draws) # Square σ2 to get σ
gini_coef <- 2 * pnorm(σ/sqrt(2), mean = 0, sd = 1) - 1 # Calculate the Gini-coefficients for each σ
hist(gini_coef, breaks=200, probability = TRUE) # Hist Gini-coefficients
lines(density(gini_coef), col='red') # Plot density curve

# c)
# Functions
eqtail_CI <- function(..X, interval) {
  lower <- (1-interval)/2
  upper <- 1 - lower
  n <- length(..X)
  X <- sort(..X) # Sort from smallest to largest value
  
  return (list(lower=X[n*lower], upper=X[n*upper]))
}

HPD <- function(density, interval) {
  gini_df <- data.frame(x = density$x, y = density$y)
  gini_df <- gini_df[with(gini_df, order(y)),]
  gini_df <- data.frame(x = gini_df$x, y = gini_df$y)
  n <- dim(gini_df)[1]
  lower <- 1 - interval
  print(lower)
  HPD_cumsum <- cumsum(gini_df$y)/sum(gini_df$y)
  HPD_lower <- which(HPD_cumsum >= lower)[1]
  gini_df_95 <- gini_df[(HPD_lower + 1):n, ]
  HPD_interval <- c(min(gini_df_95$x), max(gini_df_95$x))
  return (list(lower = HPD_interval[1], upper = HPD_interval[2]))
}

# 95% equal tail credible interval
gini_coef_CI <- eqtail_CI(gini_coef, 0.95)

hist(gini_coef, breaks=200, probability = TRUE) # Hist Gini-coefficients
lines(density(gini_coef), col='red') # Plot density curve
abline(v = gini_coef_CI$lower, col='blue') # Plot lower CI line
abline(v = gini_coef_CI$upper, col='blue') # Plot upper CI line

# Highest Posterior Density Interval
gini_density <- density(gini_coef)
HPD_interval <- HPD(gini_density, 0.95)

# Plot histogram of Gini coefficients, 95% credible interval (blue) and 95% HPD interval (green)
hist(gini_coef, breaks=200, probability = TRUE) # Hist Gini-coefficients
lines(density(gini_coef), col='red') # Plot density curve
abline(v = gini_coef_CI$lower, col='blue') # Plot lower CI line
abline(v = gini_coef_CI$upper, col='blue') # Plot upper CI line
abline(v = HPD_interval$lower, col='green') # Plot lower HPD interval
abline(v = HPD_interval$upper, col='green') # Plot upper HPD interval

############# Task 3 #############
# von Mises distribution looks like a normal distribution with a spiky top and
# is a continues probability distribution on the circle, where theta is an angle.
# Kappa (κ): Concentration parameter. Large κ gives small variance around µ.
# Likelihood: p(y_1, ..., y_n | µ, κ) = exp[ κ * cos(y - µ) ]/(2πIo(κ))
# Prior: κ ~ Exp(λ = 1), mean = 1/λ

# Setup
# Wind-angles in degrees on 10 different days
# North is zero
# 
y.degrees <- c(40, 303, 326, 285, 296, 314, 20, 308, 299, 296) 
y.radians <- c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mu <- 2.39 # Mean directon

# a)
# Prior
kappa <- data.frame(seq = seq(0, 10, 0.01),
                    posterior = 0
)

for (i in 1:dim(kappa)[1]) { # Loop over every kappa
  k <- kappa$seq[i] # Extract current kappa
  prior <- exp(-k) # Calculate prior with current kappa
  bessel <- besselI(x = k,
                    nu = 0) # Bessel-function
  likelihood <- prod(exp(k * cos(y.radians - mu))/(2*pi*bessel)) # Calculate von Mises probability
  kappa$posterior[i] <- likelihood * prior # Calculate posterior with current kappa
}

# Plot posterior for different kappas
plot(kappa$seq, kappa$posterior, type='l')

# b)
index <- which.max(kappa$posterior) # Finds index with maximum posterior
kappa.mode <- kappa$seq[which.max(kappa$posterior)] # Extract kappa with maximum posterior

plot(kappa$seq, kappa$posterior, type='l')
abline(v = kappa.mode, col='red')

#####################################################################################################

########################## Lab 2 ##########################
## Meta-info
# Linear regression
# Polynomial regression
# Logistic regressions
# Credible interval
# Maximum likelihood
# Optim
# Hessian
# Mode of beta
# Predictive distributions

# Task 1: Linear and polynomial regression
# a) Set the prior hyperparameters µ0, Ω0, ν0 and σ2 to sensible values
# b) Check if your prior from a) is sensible
# Simulate draws from joint prior and compute regression curve of each draw
# c) Simulates from the joint posterior distribution of β0, β1,β2 and σ2
# Plot: 
# • Posterior mean of simulations
# • Curve of lower and upper credible interval of f(time)
# d) Simulate highest expected temperatures
# Simulate from posterior distribution of time with highest expected temperatures
# What to do to mitigate risk of overfitting high order polynomial regression?

# Task 2: Posterior approximation for classification with logistic regression
# a) Fit logistic regression using maximum likelihood estimations
# b) Approximate the posterior distribution of the 8-dim parameter vector β with a multivariate normal distribution
# c) Simulates from the predictive distribution of the response variable in a logistic regression

#####################################################################################################

# Functions
scal_inv_schsq <- function(v, σ_2, nDraws) {
  X <- rchisq(n = nDraws, df = v)
  return (v*σ_2/X)
}

############# Task 1 #############
# Respons variable: temp = β0 + β1*time + β2 * time^2 + ε, ε ~ N(0, σ2)
# Covariate: time = No. of days since beginning of year/366

# a)
# Conjugate priors
# β | µ ~ N(µ0, σ2*Ω0_(-1))
# σ2 ~ Inv-X(v0, σ0_2)
#
# Set prior hyperparameters:
# µ0: The expected value of the betas [vector]
# Ω0: How sure we are about the betas (scales the variance) [matrix]
# ν0: How sure we are about our prior knowledge of the sigmas [scalar]
# σ0_2: The variance of the betas [vector]

dataset <- read.table("TempLinkoping.txt", header = TRUE)

prior.mu0 <- c(-5, 100, -100)
prior.v0 <- 10
prior.σ0_2 <- (7/1.96)^2 # 10 = 1.95 * σ -> 10 degrees are in the confidence interval 95% of the times
prior.Ω0 <- matrix(c(0.5, 0, 0, 0, 0.1, 0, 0, 0, 0.1),
                   nrow = 3,
                   ncol = 3)
prior.inv_Ω0 <- solve(prior.Ω0)

# b)
nDraws = 1000
beta_draws <- matrix(nrow = nDraws, ncol = 3)

plot.new()
plot.window(xlim=c(0,1), ylim=c(-20, 30))
axis(side=1)
axis(side=2)

# Sample betas & plot regression curves
for (i in 1:nDraws) {
  sigma2 <- scal_inv_schsq(prior.v0, prior.σ0_2, 1) # Sample sigma2 from scaled inverse chi squared
  beta_draws[i,] <- rmvnorm(n = 1, mean = prior.mu0, sigma = sigma2*prior.inv_Ω0) # Sample betas from Multinormal
  
  # Generates regression point for each dataset$time. All of them are then plotted as a line
  lines(dataset$time, 
        beta_draws[i,1] + beta_draws[i,2]*dataset$time + beta_draws[i,3]*dataset$time^2,
        col=rgb(0, 0, 0, 0.1))
}

lines(dataset$time,
      mean(beta_draws[,1]) + mean(beta_draws[,2])*dataset$time + mean(beta_draws[,3])*dataset$time^2,
      col=rgb(1, 0, 0, 1))

# Priors not sensible. New priors
prior.mu0 <- c(-5, 100, -100)
prior.v0 <- 10
prior.σ0_2 <- (7/1.96)^2 # 7 = 1.95 * σ -> 7 degrees are in the confidence interval 95% of the times
prior.Ω0 <- matrix(c(0.5, 0, 0, 0, 0.1, 0, 0, 0, 0.1),
                   nrow = 3,
                   ncol = 3)
prior.inv_Ω0 <- solve(prior.Ω0)
n <- dim(dataset)[1]

# c) 

# Setup
X <- cbind(1, dataset$time, dataset$time^2) # Create X with constant row for beta_0
Y <- dataset$temp
n <- dim(dataset)[1]
CI <- 0.90
CI_lower <- (1-CI)/2 # 0.05
CI_upper <- 1-(1-CI)/2 # 0.95

# Posterior
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%Y # Beta_hat by classic by OLS (Ordinary Least Square)
posterior.mun <- solve(t(X) %*% X + prior.Ω0) %*% (t(X) %*% X %*% beta_hat + prior.Ω0 %*% prior.mu0) # Mu_n
posterior.Ωn <- t(X)%*%X + prior.Ω0 # Omega_n
posterior.inv_Ωn <- solve(posterior.Ωn) # Inverse Omega_n
posterior.vn <- prior.v0 + n # v_n
posterior.sigman_2 <- (t(Y) %*% Y + t(prior.mu0) %*% prior.Ω0 %*% prior.mu0 - t(posterior.mun) %*% posterior.Ωn %*% posterior.mun)/posterior.vn

nDraws <- 10000
beta_post_draws <- matrix(nrow = nDraws,
                          ncol = 3)
posterior_y <- matrix(nrow = nDraws,
                      ncol = n)

# Draws of sigma2 and beta posteriors
for (i in 1:nDraws) {
  sigma2 <- as.vector(scal_inv_schsq(posterior.vn, posterior.sigman_2, 1))
  beta_post_draws[i,] <- rmvnorm(n = 1,
                                 mean = posterior.mun,
                                 sigma = sigma2*posterior.inv_Ωn)
}

# Generates a large matrix (10000 x 366), For each time/column (366) -> 10000 predicted temps/rows, one for each beta-draw
temp_pred_each_time <-  t(X %*% t(beta_post_draws))
temp_pred_each_time_sorted = apply(X = temp_pred_each_time, 
                                   MARGIN = 2, 
                                   sort) # Sort each time/row

# Extract lower and upper
temp_pred_CI90 <- rbind(temp_pred_each_time_sorted[round(nDraws*CI_lower),], 
                        temp_pred_each_time_sorted[round(nDraws*CI_upper),]) 

# Plot data and mean regression line
plot(dataset)
lines(dataset$time, mean(beta_post_draws[,1]) + mean(beta_post_draws[,2]) * dataset$time + mean(beta_post_draws[,3]) * dataset$time^2,
      col=rgb(1, 0, 0, 1))
lines(dataset$time, temp_pred_CI90[1,], col=rgb(0, 1, 0, 1))
lines(dataset$time, temp_pred_CI90[2,], col=rgb(0, 1, 0, 1))

# d)
# Time with highest expected temperature: time = -B1/2B2
# Calculated from the derivation of f(time) set to 0. 

posterior_max_time <- -beta_post_draws[,2]/(2*beta_post_draws[,3])
abline(v = mean(posterior_max_time), col=rgb(0, 0, 1, 1))

# e)
# To mitigate the risk of overfitting due to higher order polynomials we want to have a small
# mu_n for larger betas. This will lead to a smaller coefficient for the higher polynomials,
# thus making them affect the end result less.
# 1) Set my_0 corresponding to higher polynomials low
# 2) Set the diagonal indices of Ometa_0 corresponding higher polynomial high

############# Task 2 #############
# Posterior approximation for classification with logistic regression
# Dataset:
# Response variable: Work
# Covariates: Constant  HustbandInc EducYear  ExpYear   ExpYear2  Age   NSmallChild   NBigChild

dataset <- read.table("WomenWork.dat", header=TRUE)

# a)

logRegFit <- glm(formula = Work ~.-Constant, data = dataset, family = "binomial")

# b)
# Approx posterior distribution of 8-dim parameter vector β
# Posterior: β | y, X ~ N(Beta_mode, Inv_Hessian_At_Beta_bode)
# Likelihood: y | β, X = exp(x_i*β)^(y_i)/[ 1+exp(x_i*β) ]
# Prior: β ~ N(0, τ_2*I), τ = 10

# Functions
# !!!!! Important to remember !!!!!
## 1) Always use log posterior as it's more stable and avoids problems with to small or large numbers
## 2) Don't forget that in log -> posterior = log.likelihood + log.prior
## 3) Don't forget to handle Infinity
postLogReg <- function(β, mu_0, X, Y, tau) {
  no_of_betas <- length(β) # Number of covariates
  sigma <- tau^2 * diag(no_of_betas) # Calculate sigma tau^2 * I
  
  # Likelihood
  # Logarithm of prod[ exp(x*β)^Y / (1 + exp(x*β))]
  log.likelihood <- sum(t(Y) %*% X %*% β) - log(prod(1 + exp(X %*% β)))
  
  if (abs(log.likelihood) == Inf) log.likelihood = -20000;
  
  # Prior
  log.prior <- dmvnorm(β, mean = mu_0, sigma = sigma, log = TRUE)
  
  return (log.likelihood + log.prior)
}

# Setup
n_parameters <- dim(dataset[,-1])[2] # No. of covariates
X <- as.matrix(dataset[, -1])
Y <- as.matrix(dataset[, 1])
nDraws = 10000
CI_interval = c(0.025, 0.975)

# Prior
β0 <- as.matrix(rep(0, n_parameters)) # Initial beta-values
beta.prior.mu_0 <- rep(0, n_parameters) # Mu_0
beta.prior.tau <- 10 # Tau: Given in the task

# Optim
# par: Initial values for parameter to me optimized
# fn: Function to be minimized/maximized
# Variables: Pass all variables except the one being maximized
optim.res <- optim(par = β0, 
                   fn = postLogReg, 
                   gr = NULL,
                   mu_0 = beta.prior.mu_0, 
                   X = X, 
                   Y = Y, 
                   tau = beta.prior.tau,
                   method = "BFGS", 
                   control = list(fnscale=-1), 
                   hessian = TRUE
)

β.mode <- optim.res$par # Mode of beta
β.hessian.neg.inv <- -solve(optim.res$hessian) # Negative inverse hessian of beta

# Beta draws
β.draws <- matrix(nrow = nDraws,
                  ncol = n_parameters)

for (i in 1:nDraws) {
  β.draws[i, ] <- mvrnorm(n = 1, mu = β.mode, Sigma = β.hessian.inv)
}

# Calculate Credible Interval of NSmallChild
NSmallChild.draws <- sort(β.draws[, 7]) # Sort all draws in ascending order
NSmallChild.CI <- c(NSmallChild.draws[round(nDraws*CI_interval[1])],
                    NSmallChild.draws[round(nDraws*CI_interval[2])]) # Extract Credible intervals

# Hist of draws of NSmallChild
breaks <- 200
h <- hist(β.draws[,7], breaks = breaks, plot = FALSE)
cut <- cut(h$breaks, c(NSmallChild.CI[1], NSmallChild.CI[2]))
plot(h, 
     col = cut, 
     main = "Draws of NSmallChild",
     xlab = "Value")
abline(v = NSmallChild.CI[1], col = rgb(1, 0, 0, 1))
abline(v = NSmallChild.CI[2], col = rgb(1, 0, 0, 1))

# c)

# Functions
sigmoid <- function(x) {
  return (exp(x) / (1 + exp(x)))
}

drawPredDist <- function(βmode, negInvHess, y, nDraws) {
  beta_draws <- rmvnorm(n = nDraws,
                        mean = βmode,
                        sigma = negInvHess) # Draw from posterior beta distribution
  
  return (sigmoid(beta_draws %*% y))
}

# Setup
target <- c(1, 10, 8, 10, (10/10)^2, 40, 1, 1) # Target woman covariates
pred <- numeric() # Vector to collect results
nDraws <- 10000 # No. of draws

# Distribution of logistic regression of target
pred_work <- drawPredDist(βmode = β.mode, 
                           negInvHess = β.hessian.neg.inv,
                           y = target,
                           nDraws = nDraws)

# Histogram of logistic regression
hist(pred_work, breaks=100)

# ~99.5% of the values are below 0.5.
# The data is indicating that the target woman doesn't work.
percent_below_05 <- sum(ifelse(pred_work < 0.5, 1, 0))/length(pred_work)

#####################################################################################################

########################## Lab 3 ##########################
## Meta-info
# Gibbs sampling
# MCMC
# Mixure normal models
# STAN
# Credible interval
# Effective samples

# Task 1: Normal model, mixture of normal model with semi-conjugate prior
# a) Implement Gibbs sampler that simulate from a joint posterior. Data is normal
# Evaluate Gibbs model graphically
# b) Implement Gibbs sampler that simulates form a joint posterior when mixture normal models
# Evaluate convergance of model

# Task 2: Time series models in Stan
# a) Write a function in R that simulates data from given AR(1)-process
# b) Treat µ, φ and σ2 as unknowns and estimate them using MCMC. 
# Report credible interval and effective samples
# c) Same as a) and b) but with new data
# d) Change σ2 to formative prior

#####################################################################################################

############# Task 1 #############

### Setup
x = read.table("rainfall.dat", header=TRUE)[,1]
n = length(x)

## prior
mu0 = mean(x)
tau0 = 10
v0 = 5 # << n so should not matter too much (w close to 1)
sigma0 = 40

### Functions

### Implementation
# a)
# Assume: y_1, ..., y_n | µ, σ2 ~ N(µ, σ2), (iid), µ and σ2 are unknown
#
# Prior:
# µ ~ N(µ0, τ0_2)
# σ2 ~ Inv-X(v0, σ0_2)
#
# Posterior:
# µ | σ2, y ~ N(µn, taun_2)
# σ2 | µ, y ~ Inv-X[ vn, (v0 * σ0_2 + sum[ (y_i - µ)^2 ]) / (n + v0) ]

# to keep track of samples
mu.samples = numeric()
sigma2.samples = numeric()

# initialize with prior values (any starting point is sufficient)
mu.sample = mu0
sigma2.sample = sigma0

for (i in 1:10000) {
  ## mu sampling
  # mun
  w = (n / sigma2.sample) / ( n / sigma2.sample + 1 / tau0^2 )
  mun = w * mean(x) + (1 - w) / mu0
  # taun2
  taun2 = 1 / (n / sigma2.sample + 1 / tau0^2)
  # sample mu
  mu.sample = rnorm(1, mean = mun, sd = sqrt(taun2))
  #mu.sample = rexp(1, rate = 1 / mun) if you want to try exponential
  mu.samples[i] = mu.sample
  
  ## sigma sampling
  # vn
  vn = n + v0
  # sigma2n
  sigma2n = (v0 * sigma0^2 + sum((x - mu.sample)^2)) / (n + v0)
  chi2 = rchisq(vn, n=1)
  sigma2.sample = vn * sigma2n / chi2
  sigma2.samples[i] = sigma2.sample
}
# plot of the trajectories
plot(mu.samples, type="l")
plot(sigma2.samples, type="l")

# histogram of the distribution
hist(mu.samples)
hist(sigma2.samples)

#b)
# Daily precipitation follows an iid two-component mixture of normal distributions
# p(y_i | µ, σ2, π) = π * N(y_i | µ1, σ1_2) + (1 - π) * N(y_i | µ2, σ2_2)
# µ = (µ1, µ2)
# σ2 = (σ1_2, σ2_2)
#
# π ~ Beta(α, β)
#
# Gibbs sampling
# π | I, x ~ Beta (α1 + n1, α2 + n2)
#
# Model 1:
#     σ1_2 | I, x ~ Inv-X(Vn1, σn1_2)
#     µ1 | I, σ1_2, x ~ N(µn1, σ1_2/κn1)
#
# Model 2:
#     σ2_2 | I, x ~ Inv-X(Vn_2, σn2_2)
#     µ2 | I, σ2_2, x ~ N(µn2, σ2_2/κn1)
#     I | π, µ1, σ1_2, µ2, σ2_2, x ~ Bern(θ_i), i = 1, ..., n

# Data options
x = as.matrix(x)
# x = as.matrix(rexp(6000, rate=0.025)), interesting comparison

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10 * rep(1,nComp) # Dirichlet(alpha)
muPrior <- mu0 * rep(1,nComp) # Prior mean of mu
tau2Prior <- tau0^2 * rep(1,nComp) # Prior std of mu
sigma2_0 <- sigma0^2 * rep(1,nComp) # Best guess
nu0 <- v0 * rep(1,nComp) # degrees of freedom for prior on sigma2

# MCMC options
# Reaches reasonable distribution at around 35 iterations
nIter <- 100 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 2 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
# param = c(alpha1 + n1, alpha2 + n2)
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
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

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))

mus = numeric()
for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print("nAlloc")
  print(nAlloc)
  # Update components probabilities
  # alpha = (alpha, beta)
  # nAlloc = (s, f)
  # As nAlloc is changing, pi will tend to move towards the model that includes the most
  # observations
  pi <- rDirichlet(alpha + nAlloc)
  print("Pi")
  print(pi)
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
    mus[k] = mu[j]
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    #Sys.sleep(sleepTime)
  }
  
}

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

# c)
par(mfrow = c(1, 1))
hist(x, freq=FALSE, breaks=20)
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
mu = mean(mu.samples)
sigma2 = mean(sigma2.samples)
lines(xGrid, dnorm(xGrid, mean = mu, sd=sqrt(sigma2)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram", "Mixture density", "Regular Gibbs density"), col = c("black", "red", "blue"), lwd = 2)



############# Task 2 #############

# a) Write a function that simulates from an AR(1)-process:
# x_t = µ + Φ(x_t-1 - µ) + ε_i, ε_i ~ N(0, σ2)
#
# AR(1)-process is stable when -1 <= Φ <= 1
# When Φ is either -1 or 1, the process seems to follow a pattern,
# when Φ is around 0 it seems ot be random. This is reasonable:
# When Φ becomes smaller and smaller, the previous time-series doesn't
# affect the future time-series as much. Instead, the error ε
# will affect the future time-series. As ε is a normally distributed
# error, it will appear to be random.

# Function
AR1_process <- function(mu, sigma2, phi, T) {
  X <- numeric()
  X[1] <- mu
  
  for (t in 2:T) {
    ε <- rnorm(n = 1,
               mean = 0,
               sd = sqrt(sigma2))
    X[t] <- mu + phi*(X[t-1] - mu) + ε
  }
  return(X)
}

# Setup
µ <- 10
σ2 <- 2
T <- 200
Φ <- seq(-1, 1, 0.1)

X <- matrix(nrow = length(Φ),
            ncol = T)

for (t in 1:length(Φ)) {
  X[t, ] <- AR1_process(µ, σ2, Φ[t], T)
}

par(mfrow = c(4, 4))
for (t in 1:length(Φ)) {
  main_legend <- paste(c("Φ:",Φ[t]), collapse= " ")
  plot(X[t,], 
       ylab = "x_t",
       xlab = "Index: t",
       main = main_legend)
}

# b)
library(rstan)
# Setup
µ <- 10
σ2 <- 2
T <- 200
Φ1 <- 0.3
Φ2 <- 0.95

# AR(1)-processes for each phi.
# Create ground truth dataset
X_Φ1 <- AR1_process(µ, σ2, Φ1, T)
X_Φ2 <- AR1_process(µ, σ2, Φ2, T)

# Stan-code
# You need to define the variables 
stan_program = '
data {
int<lower=0> T;
vector[T] X;
}
parameters {
real mu;
real<lower=0> sigma;
real phi;
}
model {
for (t in 2:T) 
X[t] ~ normal(mu + phi * (X[t-1] - mu), sigma);
}
'
nIter <- 4000
burnIn <- 1000
# To find information within the stan.res-variable:
# summary(stan.res.Φ1): All information
# summary(stan.res.Φ1)$summary: Just the summary
# extract(stan.res.Φ1): The estimastes of variables

# The percentages in the summary are the credible intervals of the
# covariates. When phi is large, the CI-interval of mu will be very wide.
# This is because it will be affected of previous time-series alot.

# Effective posterior samples:
# When the observations in a sample distribution is correlated or weighted,
# the effective sample size is calculated.
# Effective sample size is thought of as an "Exchange rate" between
# dependent and independent samples. For example, 1000 dependent samples
# are worth as much as 80 independent samples.
# If the samples are independent, the effective sample size equals the
# actual sample size.
stan.res.Φ1 <- stan(model_code = stan_program,
                    data = list(X = X_Φ1, N = T),
                    warmup = burnIn,
                    iter = nIter)
stan.res.Φ2 <- stan(model_code = stan_program,
                    data = list(X = X_Φ2, N = T),
                    warmup = burnIn,
                    iter = nIter)

# Plot
par(mfrow = c(1, 2))
plot(extract(stan.res.Φ1)$mu, extract(stan.res.Φ1)$phi)
plot(extract(stan.res.Φ2)$mu, extract(stan.res.Φ2)$phi)
plot(extract(stan.res.Φ1)$mu, type='l')
plot(extract(stan.res.Φ2)$mu, type='l')

# c)
# campy.dat: No. of cases of campylobacter infections in the north of
# province Quebec (Canada). Data recorded in 4 weeks intervals from
# January 1990 to end of October 2000.
# 13 observations/year
# 140 observations in total

# c_t: No. of infections at each timepoint.
# c_t follows an indenepdent Poisson distribution when conditioned on a
# latent AR(1)-process x_t:
#
#     c_t | x_t ~ Poisson(exp(x_t))
#
# x_t: Same AR(1)-process as in a)

X.campy <- read.table("campy.dat", header = TRUE)
T <- dim(X.campy)[1]

stan_program = '
data {
int<lower=0> T;
int c[T];
}
parameters {
real phi;
real sigma;
real mu;
vector[T] X;
}
model {
phi ~ uniform(-1, 1);
for (t in 2:T)
X[t] ~ normal(mu + phi*(X[t-1] - mu), sigma);

for (t in 1:T)
c[t] ~ poisson(exp(X[t]));
}
'

# Setup STAN
nIter <- 4000
burnIn <- 1000

# Stan
stan.res.campy <- stan(model_code = stan_program,
                       data = list(c = X.campy[, 1], T = T),
                       warmup = burnIn,
                       iter = nIter)

# Mean of c_i. Calculated by taking the mean of exp(X), where X is generated by STAN.
mean_stan <- exp(apply(extract(stan.res.campy)$X, MARGIN = 2, mean))

# Credible interval
stan.summary <- summary(stan.res.campy)$summary
CI_lower <- exp(stan.summary[-(1:3), "2.5%"])
CI_upper <- exp(stan.summary[-(1:3), "97.5%"])

# Plot
plot(X.campy[,1]) # Plot data
lines(mean_stan, col = rgb(1, 0, 0, 1)) # Plot mean of stan
lines(CI_lower, col = rgb(0, 1, 0, 1)) # Plot lower Credible interval
lines(CI_upper, col = rgb(0, 1, 0, 1)) # Plot upper Credible interval

# Hist of theta_i
hist(log(mean_stan), breaks = 200)

# d)
# A Scaled inverse chi square prior has been added to sigma2.

stan_program = '
data {
int<lower=0> T;
int c[T];
}
parameters {
real phi;
real sigma2;
real mu;
vector[T] X;
}
model {
phi ~ uniform(-1, 1);
sigma2 ~ scaled_inv_chi_square(T, 0.03);

for (t in 2:T)
X[t] ~ normal(mu + phi*(X[t-1] - mu), sqrt(sigma2));

for (t in 1:T)
c[t] ~ poisson(exp(X[t]));
}
'

# Setup STAN
nIter <- 4000
burnIn <- 1000

# Stan
stan.res.campy <- stan(model_code = stan_program,
                       data = list(c = X.campy[, 1], T = T),
                       warmup = burnIn,
                       iter = nIter)

# Mean of c_i. Calculated by taking the mean of exp(X), where X is generated by STAN.
mean_stan <- exp(apply(extract(stan.res.campy)$X, MARGIN = 2, mean))

# Credible interval
stan.summary <- summary(stan.res.campy)$summary
CI_lower <- exp(stan.summary[-(1:3), "2.5%"])
CI_upper <- exp(stan.summary[-(1:3), "97.5%"])

# Plot
plot(X.campy[,1]) # Plot data
lines(mean_stan, col = rgb(1, 0, 0, 1)) # Plot mean of stan
lines(CI_lower, col = rgb(0, 1, 0, 1)) # Plot lower Credible interval
lines(CI_upper, col = rgb(0, 1, 0, 1)) # Plot upper Credible interval

# Histogram
hist(log(mean_stan), breaks = 200)

#####################################################################################################

########################## Lab 4 ##########################
## MCMC
## Metropolis algorithm
## Poisson regression
## Hessian
## Mode
## Maximum likelihood estimator
## Program general functions

# Task 1: Poisson regression - the MCMC way
# a) Obtain the maximum likelihood estimator of β in the Poisson regression model
# Find significant covariates
# b) Bayesian analysis of the Poisson regression
# Find mode and hessian of Beta
# c) Simulate from the actual posterior of β using the Metropolis algorithm and 
# compare with the approximate results in b)
# Program general function
# d) Use MCMC draws from c) to simulate from the predictive distribution of
# the number of bidders in a new auction

#####################################################################################################

library(mvtnorm)

############# Task 1 #############
# Consider the following Poisson regression model:
#   y_i | β ~ Poisson[ exp(t(x_i) * β ) ], i = 1, ..., n

dataset <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
X <- dataset[, -1]
Y <- dataset[, 1]

# a)
# Obtain the maximum likelihood estimatior of β

# Fit Poisson model
nBids.fitted <- glm(formula = nBids ~.-Const,
                    data = dataset,
                    family = poisson)

# Find significant coefficients
# 99.9%: (Intercept), VerifyID, Sealed, LogBook, MinBidShare
# 95%: MajBlen
nBids.summary <- summary(nBids.fitted)

# b)
# Bayesian Analysis of the Poisson regression
# Zellner's prior:
#     β ~ N(0, 100 * (t(X)*X)^(-1)), X: n x p covariate matrix
#
# Approximate posterior density:
#     β | y ~ N[ B_hat, Jy(β_hat)^(-1) ]

# Functions
logPostPois <- function(β, mu, sigma, y, X) {
  n <- length(Y)
  
  # log likelihood
  log.likelihood <- sum(y * X %*% β - exp(X %*% β))
  
  # If log.likelihood -Inf or Inf
  if (abs(log.likelihood) == Inf) log.likelihood = -20000;
  
  # log prior
  # OBS: Use dmvnorm!! You want the likelihood of the betas!
  log.prior <- dmvnorm(x = β, mean = mu, sigma = sigma, log = TRUE)
  
  return (log.likelihood + log.prior)
}

# Setup
X.matrix <- as.matrix(X)
beta_init <- as.matrix(10 * rep(1, dim(X.matrix)[2]))
sigma <- 100 * solve(t(X.matrix)%*%(X.matrix))
mu <- rep(0, dim(sigma)[2])
optim.res <- optim(par = beta_init,
                   fn = logPostPois,
                   gr = NULL,
                   mu = mu,
                   sigma = sigma,
                   y = Y,
                   X = X.matrix,
                   method = "BFGS",
                   control = list(fnscale = -1),
                   hessian = TRUE
)

β.mode <- optim.res$par # Mode of betas
β.neg.inv.hessian <- -solve(optim.res$hessian) # Negative inverse hessian, when betas = mode

# c) Metropolis

# Functions
Metropolis <- function(theta, logPostFunc, c, Σ, nDraws, warmUp, ...) {
  theta_c <- theta # theta_c: Current thetas
  thetas <- matrix(nrow = nDraws,
                   ncol = length(theta))
  sigma_c <- c * Σ
  for (i in -warmUp:nDraws) {
    
    # Generate new thetas in the surounding area of theta_c
    theta_p <- as.vector(rmvnorm(n = 1, mean = theta_c, sigma = sigma_c)) 
    
    log.post.p <- logPostFunc(theta_p, ...)
    log.post.c <- logPostFunc(theta_c, ...)
    
    # Calculate the acceptance probability.
    # If there's a high probability of theta_p, accept_prob will be equal to 1 and the bern trial will 
    # draw in its favour.
    accept_prob <- min(1, exp(log.post.p - log.post.c))
    
    # Bern trial to decide if theta_c should continue as before or if a "jump" to newly generated theta_p should be made
    bern_trial <- rbinom(1, size = 1, prob = accept_prob)
    
    # Set theta_c to new value
    if (bern_trial == 1) {
      theta_c <- theta_p
    }
    
    # When burn-in is passed, start to save the theta-values
    if(i > 0) thetas[i, ] <- theta_c;
  }
  
  return(thetas)
}

# Setup
X.matrix <- as.matrix(X)
beta_init <- rep(0, dim(X.matrix)[2])
sigma <- 100 * solve(t(X.matrix)%*%(X.matrix))
mu <- rep(0, dim(sigma)[2])
c = 0.5

nDraws = 4000 # No. of draws
warmUp = round(nDraws/5) # No. of Burn-in

metropolis.res <- Metropolis(theta = beta_init,
                             logPostFunc = logPostPois,
                             c = c,
                             Σ = β.neg.inv.hessian,
                             nDraws = nDraws,
                             warmUp = warmUp,
                             mu = mu,
                             sigma = sigma,
                             y = Y,
                             X = X.matrix)

postPhi <- exp(metropolis.res) # Phi = exp(β) is often more interprentable

# Plot β
par(mfrow = c(3, 3))
for (i in 1:dim(metropolis.res)[2]) {
  main_legend <- paste(c("β", i), collapse = "")
  plot(density(metropolis.res[, i]),
       main = main_legend,
       type = 'l',
       ylab = main_legend)
}

# Plot Phis
par(mfrow = c(3, 3))
for (i in 1:dim(postPhi)[2]) {
  main_legend <- paste(c("Phi:", i), collapse = " ")
  plot(density(postPhi[, i]),
       main = main_legend,
       type = 'l',
       ylab = main_legend)
}

# d) Use draws in c) to simulate the predictive distribution of a new bid

new_bid <- c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)

# Simulation from predictive distribution by using Metropolis betas
# Simulate from Poisson distribution with lambda = exp(x_i * β)
lambda <- metropolis.res %*% new_bid
pred <- rpois(n = nDraws,
              exp(lambda))

# Hist draws
hist(pred, breaks = 30)

# Calculate probability of no bidders
prob_no_bidder <- sum(pred == 0)/length(pred)