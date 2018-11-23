#### Lab 4 Task 1 - Poisson regression, the MCMC way.
### Libraries
library(mvtnorm)

### Setup
data.ebay = read.table("eBayNumberOfBidderData.dat", header=TRUE)

### Functions
logPostPoisson = function(β, mu0, Sigma0, X, y) {
  p = length(β)
  
  # log of the likelihood
  log.likelihood = sum(y * X %*% β - exp(X %*% β))
  
  # if likelihood is very large or very small, stear optim away
  if (abs(log.likelihood) == Inf) log.likelihood = -20000;
  
  # log of the prior
  log.prior = dmvnorm(β, mean = mu0, sigma = Sigma0, log = TRUE)
  
  return(log.likelihood + log.prior)
}

### Implementation
# a)
# Fit the model
glm.fit = glm(nBids ~ .-Const, data=data.ebay ,family="poisson")

# Check significance
summary(glm.fit)
# The significant covariates are: (Intercept), VerifyID, Sealed, MajBlem (0.05), LogBook, MinBidShare

# b)
y.ebay = data.ebay[,1]
X.ebay = as.matrix(data.ebay[,-1])
# Prior
β0 = rep(0, 9) # mean prior
Sigma0 = 100 * solve(t(X.ebay) %*% X.ebay) # Sigma prior

# Optimize posterior of β given priors and data
optim.res = optim(β0, logPostPoisson, gr=NULL, 
                   β0, Sigma0, X.ebay, y.ebay,
                   method="BFGS", control=list(fnscale=-1), hessian=TRUE)

# Results
β.mode = optim.res$par
β.invhessian = -solve(optim.res$hessian)

#c)
Metropolis = function(nBurnIn, nSamples, theta, c, logPostFunc, ...) {
  # Setup
  theta.c = theta
  Sigma.c = c * β.invhessian
  nAccepted = 0
  p = length(theta)
  theta.samples = matrix(NA, nSamples, p)
  
  # Iterations
  for(i in -nBurnIn : nSamples) {
    # Sample from proposal distribution
    theta.p = as.vector(rmvnorm(1, mean = theta.c, sigma = Sigma.c))
    
    # Calculate log posteriors
    log.post.p = logPostFunc(theta.p, ...)
    log.post.c = logPostFunc(theta.c, ...)
    
    # Calculate alpha
    alpha = min(1, exp(log.post.p - log.post.c))
    
    # Select sample with probability alpha
    u = runif(1)
    if (u <= alpha){
      theta.c = theta.p
      nAccepted = nAccepted + 1
    }
    
    # Save sample if not burnin
    if (i>0) theta.samples[i,] = theta.c
  }
  cat("Sampled", nSamples, "samples with an acceptance rate of", nAccepted/nSamples)
  return(theta.samples)
}

# Setup sampling
c = 0.5
nSamples = 4000
nBurnIn = 1000

# Samples from posterior using metropolis
β.samples = Metropolis(nBurnIn, nSamples, β.mode, c, logPostPoisson, β0, Sigma0, X.ebay, y.ebay)

# Estimate parameters using posterior mean
β.post.mean = apply(β.samples, 2, mean)

# Plot parameter differences for comparision
plot(rep(0, 9), ylim = c(-0.015, 0.015), col="blue",
     xlab="β", ylab="Parameter difference", main="Parameter value differences of the different methods")
points(glm.fit$coefficients-β.mode, col="green")
points(glm.fit$coefficients-β.post.mean, col="red")
legend("bottomright", legend=c("glm", "posterior mode", "posterior mean MCMC"), lwd=2,
       col=c("blue", "green", "red"))

# Plot posterior distribution of all β
plot(density(β.samples[,1]), main="Const", xlab=expression(beta[1]), ylab="Density")
plot(density(β.samples[,2]), main="PowerSeller", xlab=expression(beta[2]), ylab="Density")
plot(density(β.samples[,3]), main="VerifyID", xlab=expression(beta[3]), ylab="Density")
plot(density(β.samples[,4]), main="Sealed", xlab=expression(beta[4]), ylab="Density")
plot(density(β.samples[,5]), main="MinBlem", xlab=expression(beta[5]), ylab="Density")
plot(density(β.samples[,6]), main="MajBlem", xlab=expression(beta[6]), ylab="Density")
plot(density(β.samples[,7]), main="LargNeg", xlab=expression(beta[7]), ylab="Density")
plot(density(β.samples[,8]), main="LogBook", xlab=expression(beta[8]), ylab="Density")
plot(density(β.samples[,9]), main="MinBidShare", xlab=expression(beta[9]), ylab="Density")

# Plot posterior distribution of all phi=exp(β)
phis = exp(β.samples)
plot(density(phis[,1]), main="Const", xlab=expression(phi[1]), ylab="Density")
plot(density(phis[,2]), main="PowerSeller", xlab=expression(phi[2]), ylab="Density")
plot(density(phis[,3]), main="VerifyID", xlab=expression(phi[3]), ylab="Density")
plot(density(phis[,4]), main="Sealed", xlab=expression(phi[4]), ylab="Density")
plot(density(phis[,5]), main="MinBlem", xlab=expression(phi[5]), ylab="Density")
plot(density(phis[,6]), main="MajBlem", xlab=expression(phi[6]), ylab="Density")
plot(density(phis[,7]), main="LargNeg", xlab=expression(phi[7]), ylab="Density")
plot(density(phis[,8]), main="LogBook", xlab=expression(phi[8]), ylab="Density")
plot(density(phis[,9]), main="MinBidShare", xlab=expression(phi[9]), ylab="Density")

# d)
# Auction vector
x = c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)

# Setup for sampling
nBids.samples = numeric()

# Sample nBids
for(i in 1:nSamples) {
  nBids.samples[i] = rpois(1, exp(x %*% β.samples[i,]))
}

# Plot distribution
barplot(table(nBids.samples))

# Calculate probabily of no bidders
prob = sum(nBids.samples == 0) / nSamples # 0.36
