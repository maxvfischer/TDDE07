#### Lab 3 Task 2
### Libraries
library(rstan)

### Setup
mu = 10
sigma2 = 2
T = 200

### Functions
AR.simulation = function(mu, sigma2, T, phi) {
  x.sample = numeric()
  x.sample[1] = mu
  for(i in 2:T) {
    x.sample[i] = mu + phi * (x.sample[i-1] - mu) + rnorm(1, mean = 0, sd = sqrt(sigma2))
  }
  return(x.sample)
}

### Implementation
# a)

phis = seq(-1,1,0.05)
n = length(phis)
x.samples = matrix(nrow=n, ncol=T, 0)
for(i in 1:n) {
  x.samples[i,] = AR.simulation(mu, sigma2, T, phis[i])
}
# plot results
par(mfrow=c(2,2))
plot(x.samples[1,], xlab="t", ylab="x_t", main="phi=-1")
plot(x.samples[15,], xlab="t", ylab="x_t", main="phi=--0.3")
plot(x.samples[27,], xlab="t", ylab="x_t", main="phi=0.3")
plot(x.samples[40,], xlab="t", ylab="x_t", main="phi=0.95")

# b)
# sample AR processes with phi=0.3 and phi=0.95
x.03 = AR.simulation(mu, sigma2, T, phi=0.3)
x.95 = AR.simulation(mu, sigma2, T, phi=0.95)

ARStanModel = 'data {
  int<lower=0> N;
  vector[N] x;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}
model {
  for (n in 2:N)
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma);
}'

# perform MCMC
x.fit.03 = stan(model_code=ARStanModel, data=list(x=x.03, N=T))
x.fit.95 = stan(model_code=ARStanModel, data=list(x=x.95, N=T))

# extract summary
summary.03 = summary(x.fit.03)
summary.95 = summary(x.fit.95)
# extract n_eff
n_eff.03 = summary.03$summary[,"n_eff"]
n_eff.95 = summary.95$summary[,"n_eff"]

## number of effective samples (out of 4000):
n_eff.mu.03 = n_eff.03["mu"] #1345
n_eff.phi.03 = n_eff.03["phi"]# 1351
n_eff.sigma.03 = n_eff.03["sigma"] #1836
n_eff.mu.95 = n_eff.95["mu"] #1369
n_eff.phi.95 = n_eff.95["phi"] #1365
n_eff.sigma.95 = n_eff.95["sigma"] #1937

# extract params
x.params.03 = extract(x.fit.03)
x.params.95 = extract(x.fit.95)

## credible intervals
probs = c(0.025,0.975)
CI.mu.03 <- apply(as.matrix(x.params.03$mu), 2, quantile, probs=probs) # (9.709384, 10.208119)
CI.mu.95 <- apply(as.matrix(x.params.95$mu), 2, quantile, probs=probs) # (-5.017947, 21.271056)
CI.phi.03 <- apply(as.matrix(x.params.03$phi), 2, quantile, probs=probs) # (0.1847434, 0.4413590)
CI.phi.95 <- apply(as.matrix(x.params.95$phi), 2, quantile, probs=probs) # (0.7994717, 0.9406423)
CI.sigma.03 <- apply(as.matrix(x.params.03$sigma), 2, quantile, probs=probs) # (1.328897, 1.617287)
CI.sigma.95 <- apply(as.matrix(x.params.95$sigma), 2, quantile, probs=probs) # (0.5284631, 1.9843986)

## estimation of true values:
# posterior means
post.mean.x.03 = get_posterior_mean(x.fit.03)
post.mean.x.95 = get_posterior_mean(x.fit.95)

# post mu, phi, sigma: 0.3
mu.post.03 = post.mean.x.03[1,5] # 6.934868
phi.post.03 = post.mean.x.03[2,5] # 0.3133354
sigma.post.03 = post.mean.x.03[3,5] # 1.466067

# post mu, phi, sigma: 0.95
mu.post.95 = post.mean.x.95[1,5] # 1.249812
phi.post.95 = post.mean.x.95[2,5] # 0.8716227
sigma.post.95 = post.mean.x.95[3,5] # 1.405657

#ii)
# plot density
par(mfrow=c(1,1))
plot(x=x.params.03$mu, y=x.params.03$phi,
     xlab="mu", ylab="phi", main="Posterior density, AR with phi=0.3")
plot(x=x.params.95$mu, y=x.params.95$phi,
     xlab="mu", ylab="phi", main="Posterior density, AR with phi=0.95")

#c)

data.campy = read.table("campy.dat", header=TRUE)[,1]
N = length(data.campy)

PoissonARStanModel = '
data {
  int<lower=0> N;
  int c[N];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  vector[N] x;
}
model {
  // mu ~ normal();
  // sigma2 ~ scaled_inv_chi_square();
  // phi ~ normal(0, 0.6);
  phi ~ uniform(-1,1);
  for (n in 2:N)
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma);
  for (n in 1:N)
    c[n] ~ poisson(exp(x[n]));
}'
c.fit = stan(model_code=PoissonARStanModel, data=list(N = N, c = data.campy))

par(mfrow=c(1,1))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
c.params = extract(c.fit)
c.summary = summary(c.fit)

X.summary = summary(c.fit)$summary[-c(1,2,3),]

# extract CI and mean of x
x.upper = X.summary[,"97.5%"]
x.lower = X.summary[,"2.5%"]
x.mean = X.summary[,"mean"]

# plot data, mean and CI of the latent intensity
plot(data.campy, xlab="time", ylab="value")
lines(exp(x.upper), col="green")
lines(exp(x.lower), col="green")
lines(exp(x.mean), col="red")
# very volatile, a smoothness factor seems reasonable

# d)
SmoothPoissonARStanModel = '
data {
  int<lower=0> N;
  int c[N];
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
  vector[N] x;
}
model {
  // mu ~ normal();
  // phi ~ normal(0, 0.6);
  sigma2 ~ scaled_inv_chi_square(N, 0.03);
  phi ~ uniform(-1,1);
  for (n in 2:N)
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sqrt(sigma2));
  for (n in 1:N)
    c[n] ~ poisson(exp(x[n]));
}'
c.fit.smooth = stan(model_code=SmoothPoissonARStanModel, data=list(N = N, c = data.campy))
X.summary.smooth = summary(c.fit.smooth)$summary[-c(1,2,3),]

# extract CI and mean of x
x.upper.smooth = X.summary.smooth[,"97.5%"]
x.lower.smooth = X.summary.smooth[,"2.5%"]
x.mean.smooth = X.summary.smooth[,"mean"]

# plot data, mean and CI of the latent intensity
plot(data.campy, xlab="time", ylab="value")
lines(exp(x.upper.smooth), col="blue")
lines(exp(x.lower.smooth), col="blue")
lines(exp(x.mean.smooth), col="red")

