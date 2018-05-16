#### Lab 4 Task 1 - Poisson regression, the MCMC way.
### Libraries
library(mvtnorm)

### Setup
data.ebay = read.table("eBayNumberOfBidderData.dat", header=TRUE)

### Functions
logPostPoisson = function(β, y, X, Sigma) {
  p = length(β)
  
  # log of the likelihood
  log.likelihood = sum(y * (X %*% β) - exp(X %*% β))
  
  # if likelihood is very large or very small, stear optim away
  if (abs(log.likelihood) == Inf) log.likelihood = -20000;
  
  # log of the prior
  log.prior = dmvnorm(β, mean = matrix(0, p, 1), sigma = Sigma, log = TRUE)
  
  return(log.likelihood + log.prior)
}

### Implementation
# a)
# Fit the model
glm.fit = glm(nBids ~ .-Const, data=data.ebay ,family="poisson")

# Check significance
summary(glm.fit)
# The significant are: (Intercept), VerifyID, Sealed, MajBlem (0.05), LogBook, MinBidShare

# b)
y.ebay = data.ebay[,1]
X.ebay = as.matrix(data.ebay[,-1])
# Prior
β0 = rep(0, 9) # mean prior
Sigma0 = 100 * solve(t(X.ebay) %*% X.ebay)

# Optimize parameter values
optim.res = optim(β0, logPostPoisson, gr=NULL, 
                   y.ebay, X.ebay, Sigma0,
                   method="BFGS", control=list(fnscale=-1), hessian=TRUE)

# Results
β.mode = optim.res$par
β.invhessian = -solve(optim.res$hessian)

#c)


