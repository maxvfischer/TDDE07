# libraries
library("mvtnorm")

# data
data = read.table("WomenWork.dat", header=TRUE)

# a)
glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)
head(glmModel)
plot(glmModel$fitted.values)
pred = glmModel$fitted.values > 0.5
plot(data$NSmallChild, data$ExpYears2, col=rgb(0,1-pred,pred))

# b)
# Approximate the posterior distribution of the 8-dim parameter vector Beta
# setup data
y = data[,1]
X = data[,-1]
names = names(X)
p = dim(X)[2]
X = as.matrix(X)

# setup prior
my = rep(0, p)
tau = 10
Sigma = tau^2*diag(p)

# log posterior of a logistic regression
LogPostLogistic = function(beta, y, X, my, Sigma) {
  p = length(beta)
  ypred = X %*% beta # make prediction
  
  # Logarithm of the likelihood seen on page 5 in lecture 6.
  log.likelihood = sum(ypred * y - log(1 + exp(ypred)))
  
  # If likelihood is very large or very small, the abs(log.likelihood) will be close to infinity. In those cases we set
  # log.likelihood to -20000
  if (abs(log.likelihood) == Inf) log.likelihood = -20000;
  
  # Logarithm of the prior
  log.prior <- dmvnorm(beta, matrix(0, p, 1), Sigma, log=TRUE);
  
  return(log.likelihood + log.prior)
}

# use optim to find max
beta.init = rep(0, p)
optim.res = optim(beta.init, LogPostLogistic, gr=NULL, y, X, my, Sigma, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
beta.mode = optim.res$par #0.62672884 -0.01979113  0.18021897  0.16756670 -0.14459669 -0.08206561 -1.35913317 -0.02468351
beta.invhessian = -solve(optim.res$hessian)
#  2.266022568  3.338861e-03 -6.545121e-02 -1.179140e-02  0.0457807243 -3.029345e-02 -0.1887483542 -0.0980239285
#  0.003338861  2.528045e-04 -5.610225e-04 -3.125413e-05  0.0001414915 -3.588562e-05  0.0005066847 -0.0001444223
# -0.065451206 -5.610225e-04  6.218199e-03 -3.558209e-04  0.0018962893 -3.240448e-06 -0.0061345645  0.0017527317
# -0.011791404 -3.125413e-05 -3.558209e-04  4.351716e-03 -0.0142490853 -1.340888e-04 -0.0014689508  0.0005437105
#  0.045780724  1.414915e-04  1.896289e-03 -1.424909e-02  0.0555786706 -3.299398e-04  0.0032082535  0.0005120144
# -0.030293450 -3.588562e-05 -3.240448e-06 -1.340888e-04 -0.0003299398  7.184611e-04  0.0051841611  0.0010952903
# -0.188748354  5.066847e-04 -6.134564e-03 -1.468951e-03  0.0032082535  5.184161e-03  0.1512621814  0.0067688739
# -0.098023929 -1.444223e-04  1.752732e-03  5.437105e-04  0.0005120144  1.095290e-03  0.0067688739  0.0199722657

# simulate draws
nDraws = 10000
beta.draws = rmvnorm(n = nDraws, mean = beta.mode, sigma = beta.invhessian) # draws from posterior

# CI
beta.NSmallChild.CI = apply(X = beta.draws, MARGIN=2, FUN=function(x) quantile(x,c(0.025, 0.975), na.rm=T))[,7] #-2.1315425 -0.6067602

# c)

