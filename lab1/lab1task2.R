library(geoR)
# Log-normal distribution and the Gini coefficient.
# Assume that you have asked 10 randomly selected persons about their monthly income
# (in thousands Swedish Krona) and obtained the following ten observations: 14,
# 25, 45, 25, 30, 33, 19, 50, 34 and 67. A common model for non-negative continuous
# variables is the log-normal distribution

# a)
# Simulate 10, 000 draws from the posterior of σ2
# (assuming µ = 3.5) and compare with the theoretical 
# Inv − χ2(n, τ 2) posterior distribution.

tau2 <- function(..Y, my) {
  n <- length(..Y)
  return(sum((log(..Y)-my)^2)/n)
}

data <- c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
my <- 3.5
n <- length(data)
nDraws = 10000

draw <- rchisq(n=nDraws, df=n)
sigma2 <- n*tau2(data, my)/draw
hist(sigma2, breaks=100)

# b)
sigma = sqrt(sigma2)
G <- 2 * pnorm(sigma/sqrt(2), mean = 0, sd = 1) - 1
hist(G, breaks=100)

# c)
G.sorted <- sort(G)[(nDraws*0.025+1):(nDraws*0.975)] 
G.credible_interval <- c(min(G.sorted), max(G.sorted)) # Credible interval 
print(G.credible_interval)

G.density <- density(G)
G.density.df <- data.frame(x = G.density$x, y = G.density$y)

G.density.ordered <- G.density.df[with(G.density.df, order(y)),]
G.density.ordered <- data.frame(x = G.density.ordered$x, y = G.density.ordered$y)
print(G.density.ordered)
G.density.ordered.95 = G.density.ordered[(512*0.05+1):512,]
print(G.density.ordered.95)
