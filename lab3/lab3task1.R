#### Lab 3 Task 1

### Setup
x = read.table("rainfall.dat", header=TRUE)[,1]
n = length(x)

## prior
mu0 = mean(x)
tau0 = 10
v0 = 100 # << n so should not matter too much (w close to 1)
sigma0 = 10

### Functions

### Implementation
# a) i)
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

