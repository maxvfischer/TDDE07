#Bernoulli ... again.
#Let y1, ..., yn|θ ∼ Bern(θ), and assume that you have obtained a sample with s = 14
#successes in n = 20 trials. Assume a Beta(α0, β0) prior for θ and let α0 = β0 = 2.

# a)
# Draw random numbers from the posterior θ|y ∼ Beta(α0 + s, β0 + f), y = (y1, . . . , yn), 
# and verify graphically that the posterior mean and standard deviation
# converges to the true values as the number of random draws grows large.

# Data
s <- 14
n <- 20
f <- n - s

# Prior
prior_alpha <- 2
prior_beta <- 2

# Posterior
posterior_alpha <- prior_alpha + s 
posterior_beta <- prior_beta + f

theta.means <- c()
theta.sd <- c()
theta.n <- c()
for (n in seq(10, 10000, 100)) {
  thetas <- rbeta(n, posterior_alpha, posterior_beta) # Generate Thetas
  theta.means <- append(theta.means, mean(thetas)) # Calculate mean and add to vector
  theta.sd <- append(theta.sd, sd(thetas)) # Calculate standard deviation and add to vector
  theta.n <- append(theta.n, n) # Add number of draws to vector
}

# Plot means of thetas
plot(theta.n, theta.means, type='l', xlab='Number of draws', ylab='Mean of Theta')

# Plot standard deviation of thetas
plot(theta.n, theta.sd, type='l', xlab='Number of draws', ylab='Standard dev of Theta')

# ------------------------------------------------------------------------------------------------- #
# b)
# Use simulation (nDraws = 10000) to compute the posterior probability 
# Pr(θ < 0.4|y) and compare with the exact value [Hint: pbeta()].
