# Lab 2 Task 1
# Linear and polynomial regression
## Setup
data = read.table("TempLinkoping.txt", header=TRUE)
library(MASS)

## Implementation
# a)
# my0
B0 = -5 # prior knowledge: Mean temp of -5 degrees 1st of January
B1 = 20 # curve for 20 degrees in 30th of June
B2 = 60 # curve for 20 degrees in 30th of June
my0 = c(B0, B1, B2)
curve(-5+B1*x+B2*x^2)
points(data) # for comparision

# sigma0
sigma2_0 = 12.5771 # diff within 7 degress 95% of the cases, results in sigma of 7/1.96 --> sigma0 = 12.5771

# v0
v0 = 60 # we have lived long

# omega0
omega0 = diag(c(0.5, 0.1, 0.1)) # more certain about B0 than the others

# b)
# new hyper parameters
B0 = -5
B1 = 100
B2 = -100
my0 = c(B0, B1, B2)
sigma2_0 = 2
v0 = 100
omega0 = diag(c(0.5, 0.5, 0.5))
nDraws = 1000

# simulation
chi2 = rchisq(v0, n=nDraws)
sigma2_draws = v0*sigma2_0/chi2
omega0_inv = solve(omega0)
B_draws = matrix(0, nDraws, 3)

# start new plot
plot.new()
plot.window(xlim=c(0,1), ylim=c(-20,30))
axis(side=1)
axis(side=2)

# simulate draws
for(i in 1:nDraws){
  B_draws[i,] = mvrnorm(n = 1, mu = my0, Sigma=sigma2_draws[i]*omega0_inv)
  lines(data$time, B_draws[i,1]+B_draws[i,2]*data$time+B_draws[i,3]*data$time^2, col=rgb(0,0,0,0.1))
}
lines(data$time, mean(B_draws[,1])+mean(B_draws[,2])*data$time+mean(B_draws[,3])*data$time^2, col=rgb(1,0,0,1))

# c) 
B0 = -5
B1 = 100
B2 = -100
my0 = c(B0, B1, B2)
sigma2_0 = 2
v0 = 100
omega0 = diag(c(0.5, 0.5, 0.5))
X = cbind(1, data$time, I(data$time)^2)
Y = data$temp
B_hat = solve(t(X)%*%X)%*%t(X)%*%Y

# n
myn = solve((t(X)%*%X+omega0)) %*% (t(X)%*%X%*%B_hat + omega0 %*% my0)
omegan = t(X)%*%X + omega0
omegan_inv = solve(omegan)
vn = v0 + dim(X)[1]
sigma2_n = ((v0*sigma2_0 + t(Y)%*%Y + t(my0)%*%omega0%*%my0 - t(myn)%*%omegan%*%myn)/vn)[1,1]

# Posterior draws
chi2_n = rchisq(vn, n=nDraws)
sigma2_draws = vn*sigma2_n/chi2_n


beta_draws = matrix(0, nDraws, 3)

plot(data)
posterior_y = matrix(0, nDraws, 366)
for (i in 1:nDraws) {
  beta_draws[i,] = mvrnorm(n = 1, mu = myn, Sigma = sigma2_draws[i]*omegan_inv)
  lines(data$time, beta_draws[i,1] + beta_draws[i,2]*data$time + beta_draws[i,3]*data$time^2, col=rgb(0,0,0,0.1))
  posterior_y[i,] = beta_draws[i,1] + beta_draws[i,2]*data$time + beta_draws[i,3]*data$time^2
}

posterior_CI = apply(X = posterior_y, MARGIN=2, FUN=function(x) quantile(x,c(0.05, 0.95), na.rm=T))
lines(data$time, posterior_CI[2,], col="green")
lines(data$time, posterior_CI[1,], col="green")

# d)
# f'(time) = B_1 + 2*B_2*time = 0
# time_* = (-B_1)/(2*B_2)
plot(data)
time_star = numeric()
posterior_y = numeric()
for (i in 1:nDraws) {
  beta_draw = mvrnorm(n = 1, mu = myn, Sigma = sigma2_draws[i]*omegan_inv)
  time_star[i] = (-beta_draw[2])/(2*beta_draw[3])
  posterior_y[i] = beta_draw[1] + beta_draw[2]*time_star + beta_draw[3]*time_star^2
}
plot(time_star, posterior_y)
hist(posterior_y, breaks=20) # around 15 degrees

# e)
# replace prior with ridge regression (or laplace prior) prior, where mu0 = 0 and
# omega0 is diagonal matrix with lambda values, larger lamba results in more shrinkage = less overfitting
     