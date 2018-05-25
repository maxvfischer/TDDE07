########################## Lab1 ##########################

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

# Draws of sigma2 and beta posteriors
for (i in 1:nDraws) {
  sigma2 <- as.vector(scal_inv_schsq(posterior.vn, posterior.sigman_2, 1))
  beta_post_draws <- rmvnorm(n = 1,
                             mean = posterior.mun,
                             sigma = sigma2*posterior.inv_Ωn)
}

# Plot data and mean regression line
plot(dataset)
lines(dataset$time, mean(beta_post_draws[,1]) + mean(beta_post_draws[,2]) * dataset$time + mean(beta_post_draws[,3]) * dataset$time^2,
          col=rgb(1, 0, 0, 1))
