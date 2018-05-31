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