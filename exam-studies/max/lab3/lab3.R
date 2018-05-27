########################## Lab3 ##########################

############# Task 1 #############
# rainfall.dat:
# Daily recordings of percipitation from 1948 to 1983 (rain or snow in units of 1/100 inch,
# zero recordings are excluded).

X <- read.table("rainfall.dat")

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

# Function
scal_inv_schsq <- function(vn, s2) {
  X <- rchisq(n = 1, df = vn)
  draw <- vn*s2/X
  return (draw)
}

# Prior
µ0 <- 30
τ0_2 <- 1500
v0 <- 40
vn <- v0 + n
σ0_2 <- 100

# Setup
nDraws = 10000
n <- dim(X)[1]
sigma2.samples <- numeric()
sigma2.sample <- σ0_2
mu.samples <- numeric()
mu.sample <- µ0
X.mean <- mean(X[,1])
w_vec <- numeric()

# Loop doing Gibbs sampling
for (i in 1:nDraws) {
  w <- (n/sigma2.sample)/( (n/sigma2.sample) + (1/τ0_2) ) # Calcuate weights
  w_vec[i] <- w # Insert weights into vector
  µn <- w * X.mean + (1 - w)*µ0 # Calculate my_n
  
  taun_2 <- 1/( (n/sigma2.sample) + (1/τ0_2)) 

  mu.sample <- rnorm(n = 1, 
                        mean = µn,
                        sd = taun_2) # Calculate sample of mu
  mu.samples[i] <- mu.sample

  # Calculate sample of sigma2 by inserting sample of mu
  sigma2.sample <- scal_inv_schsq(vn, ((v0*σ0_2 + sum((X - mu.sample)^2)) / (n + v0)))
  sigma2.samples[i] <- sigma2.sample
}

# Hist and plots of mu and sigma2 samples.
hist(mu.samples)
hist(sigma2.samples)
plot(sigma2.samples, type='l')
plot(mu.samples, type='l')

# b) ############################## FEL NÅGONSTANS ##############################
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

# Setup
nDraws <- 200
n <- dim(X)[1]

## Priors

# Mu
µ0 <- c(30, 30)
w <- c(0.5, 0.5)
kappa0 <- c(10, 10)
tao0_2 <- c(40, 40)
µ.sample <- c(30, 30)

# Sigma2
v0 <- c(40, 40)
σ0_2 <- c(1500, 1500)
σ2.sample <- c(1500, 1500)

µ.samples <- matrix(data = 0,
                    nrow = nDraws,
                    ncol = 2)
σ2.samples <- matrix(data = 0,
                     nrow = nDraws,
                     ncol = 2)

# Pi
π <- 0.5
alphas <- c(10, 10)

thetas <- numeric()

nModels <- c(n*π, n*π)

for (i in 1:nDraws) {
  print("....")
  print(i)
  print(nModels)
  # Calculate theta for I ~ Bern(theta)
  for (j in 1:n) {
    dnorm_model1 <- dnorm(x = X[j, 1], mean = µ.sample[1], sd = σ2.sample[1])
    dnorm_model2 <- dnorm(x = X[j, 1], mean = µ.sample[2], sd = σ2.sample[2])
    thetas[j] <- (1 - π)*dnorm_model2/(π*dnorm_model1 + (1-π)*dnorm_model2)
  }
  
  # Sample I for each data point
  I <- rbinom(n = n,
              size = 1,
              prob = thetas)
  
  # Create matrix with col1 = 1 if data was assigned to model1, 
  # col2 = 1 of data was assigned to model2
  X.I <- matrix(data = 0,
                nrow = n,
                ncol = 2)
  X.I[which(I == 0), 1] <- 1
  X.I[which(I == 1), 2] <- 1
  
  # Sum X.I[, 1] and X.I[, 2] to find out how many data points that are assigned to each model
  n1 <- sum(X.I[, 1])
  n2 <- sum(X.I[, 2])
  nModels <- c(n1, n2) # No. of datapoints in each odel

  # Sample π
  π <- rbeta(n = 1, shape1 = nModels[1] * alphas[1], shape2 = nModels[2] * alphas[2])
  
  # Calculate mean of each model. This is done my only calculating the mean of the data points
  # included in each model
  X.mean <- c(mean(X[which(I == 0), 1]), mean(X[which(I == 1), 1]))
  
  # Calculate my_n
  w <- (nModels/σ2.sample)/((nModels/σ2.sample) + (1/tao0_2))
  myn <- w*µ0 + (1-w)*X.mean
  
  # Calculate tao_n^2
  taon_2 <- ((nModels/σ2.sample) + (1/tao0_2))
  
  # Sample µ
  µ.sample[1] <- rnorm(n = 1, 
                       mean = myn[1],
                       sd = sqrt(taon_2))
  µ.sample[2] <- rnorm(n = 1,
                       mean = myn[2],
                       sd = sqrt(taon_2))
  
  µ.samples[i,] <- µ.sample
  
  # Calculate vn
  vn <- v0 + nModels
  
  # Calculate sigma2_n
  sq_sum_of_diff1 <- sum( (X[which(I == 0), ] - µ.sample[1])^2 )
  sq_sum_of_diff2 <- sum( (X[which(I == 1), ] - µ.sample[2])^2 )  
  sq_sum_of_diff <- c(sq_sum_of_diff1, sq_sum_of_diff2)
  σn_2 <- (v0*σ0_2 + sq_sum_of_diff)/(nModels + v0)
  
  # Sample σ2
  σ2.sample <- scal_inv_schsq(vn, σn_2)
  σ2.samples[i, ] <- σ2.sample
}

# c) 

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

     