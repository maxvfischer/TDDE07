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
plot(x.samples[1,])
plot(x.samples[15,])
plot(x.samples[25,])
plot(x.samples[41,])

# b)
# sample AR processes with phi=0.3 and phi=0.95
x.03 = AR.simulation(mu, sigma2, T, phi=0.3)
x.95 = AR.simulation(mu, sigma2, T, phi=0.95)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
    x[n] ~ normal(mu + phi * x[n-1], sigma);
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
CI.mu.03 <- apply(as.matrix(x.params.03$mu), 2, quantile, probs=probs) # (5.59146, 8.25748)
CI.mu.95 <- apply(as.matrix(x.params.95$mu), 2, quantile, probs=probs) # (0.5284631, 1.9843986)
CI.phi.03 <- apply(as.matrix(x.params.03$phi), 2, quantile, probs=probs) # (0.1847434, 0.4413590)
CI.phi.95 <- apply(as.matrix(x.params.95$phi), 2, quantile, probs=probs) # (0.7994717, 0.9406423)
CI.sigma.03 <- apply(as.matrix(x.params.03$sigma), 2, quantile, probs=probs) # (1.328897, 1.617287)
CI.sigma.95 <- apply(as.matrix(x.params.95$mu), 2, quantile, probs=probs) # (0.5284631, 1.9843986)

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
par(mfrow=c(1,2))
plot(x=x.params.03$mu, y=x.params.03$phi,
     xlab="mu", ylab="phi", main="Posterior density, AR with phi=0.3")
plot(x=x.params.95$mu, y=x.params.95$phi,
     xlab="mu", ylab="phi", main="Posterior density, AR with phi=0.3")

