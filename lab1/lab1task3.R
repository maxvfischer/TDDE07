
# a)
# Input: y-vector, one kapa value, my
# Output: vonMises probability of input. Returns a scalar

likelihood2VonMises <- function(y, k, my) {
  return(
    prod(
      exp(k*cos(y-my))/(2*pi*besselI(k, nu=0))
    )
  )
}

# -----
likelihoodVonMises <- function(y, k, my) {
    return(prod(exp(k*cos(y-my))/(2*pi*besselI(k, nu=0))))
}

y.values <- c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
kapa.values <- seq(0, 10, 0.01)
my <- 2.39

# Vector of prior values, given kapas
kapa.prior <- dexp(kapa.values)

# Calculate likelihood for each kapa, given y-vector
likelihoods = numeric()
i = 1
for(k in kapa.values) {
  likelihoods[i] = likelihoodVonMises(y.values, k, my)
  i = i + 1
}

posterior = likelihoods * kapa.prior # Vector of posterior values
plot(kapa.values, posterior, type='l') # Plot

# b)
max.posterior = which.max(posterior) # Find maximum posterior value
kapa.max.posterior = kapa.values[max.posterior] # Find maximum kapa value given max value in posterior
# Max kapa: 2.12

