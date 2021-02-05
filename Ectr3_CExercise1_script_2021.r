rm(list=ls())


### Problem 1 

### (a) simulate AR(1) process with AR coefficient 0.9 

beta.1 <- 0.9
sigma2 <- 2

T.here <- 100
burn.in <- 20

set.seed(123456)   ## setting a seed gives the same set of random numbers every time you run the script
u.vector <- rnorm(T.here+burn.in, mean=0, sd=sqrt(sigma2))  ## T + burn-in observations of the disturbances

y.vector1 <- numeric(T.here+burn.in)  ## empty vector to be filled with observations
y.vector1[1] = 0  ## initialize by the unconditional mean of the process

for(t in 2:(T.here+burn.in)){
  y.vector1[t] = beta.1 * y.vector1[t-1] + u.vector[t]
}

y.vector = y.vector1[-c(1:burn.in)]
plot(ts(y.vector))




## (b) estimation of AR coefficient beta_1
y.t <- matrix(y.vector[-1], ncol=1) # vector of dependent variable observations y_t, t=2,..., T
y.lag <- matrix(y.vector[-T.here], ncol=1) # vector of regressor obervations y_{t-1}

beta.hat <- (t(y.lag)%*%y.lag)^{-1} %*% t(y.lag) %*% y.t  # ^{-1} works when inverting a scalar, not a matrix



## (c) repeat estimation but include second lag of y_t
y.t <- matrix(y.vector[-c(1:2)], ncol=1)  # vector of dependent variable observations y_t, t=3,...,T
y.lag <- cbind(
  y.vector[-c(1,T.here)],
  y.vector[-c((T.here-1):T.here)]) # matrix of regressor observations (y_{t-1} and y_{t-2})


beta.hat2 <- solve(t(y.lag)%*%y.lag) %*% t(y.lag) %*% y.t  # least squares (X'X)^{-1}X'y. %*% computes matrix product


T.here1 <- length(y.t)  # adjust sample size because of pre-sample observations 
resids <- y.t - y.lag %*% beta.hat2 # vector of residuals
sigma2.hat2 <- 1/(T.here1-2) * t(resids) %*% resids  # estimate sigma^2

cov.beta.hat <- as.numeric(sigma2.hat2) * solve(t(y.lag) %*% y.lag) # estimated covariance matrix of beta_hat
se.beta2.hat <- sqrt(diag(cov.beta.hat)) # extract standard errors of coefficient estimates

## compute t-statistic
t.beta2 <- beta.hat2[2]/se.beta2.hat[2]

## compute p-value
2*pnorm(-abs(t.beta2))


## compare output with lm
summary(lm(y.t ~ y.lag -1))  # -1 supresses the constant in least squares estimation 









############ Problem 2 ############ 

rm(list=ls())

## a) 
## first: magic numbers

T.here <- 100  ## number of observations
burn.in <- 20  ## number of burn-in obs.

C <- matrix(c(0.1, -0.2), ncol=1) # store intercept vector as (2x1)-dimensional matrix
A1 <- matrix(c(4/5, 0, 8/3, 3/10), ncol=2) # matrix is filled column-wise
Sigma.mat <- matrix(c(1, 0.9, 0.9, 2), ncol=2)

K <- nrow(A1)  ## dimension of the VAR

mu.y <- solve(diag(2)-A1) %*% C  # unconditional mean of the process. solve() computes matrix inverse

set.seed(23456)
e.matrix <- matrix(rnorm((T.here+burn.in)*K), nrow=K)  # draw standard normally distributed iid random numbers
y.matrix <- matrix(0, nrow=K, ncol=(T.here+burn.in)) # empty (2x(T+burn.in))- matrix of zeros

y.matrix[,1] = mu.y ## initialize with unconditional mean

chol.Sigma <- chol(Sigma.mat)  # note: chol() produces an upper diagonal matrix 
t(chol.Sigma)%*%chol.Sigma # check whether LL'=Sigma_u



## recursion to fill matrix y.matrix
for(t in 2:(T.here+burn.in)){
  y.matrix[,t] = C + A1 %*% y.matrix[,(t-1)] + t(chol.Sigma) %*% e.matrix[,t]
}

ts.y.matrix <- ts(t(y.matrix), start =c(1990, 1), frequency=4) # store y.matrix as time series object

plot(ts.y.matrix) ## two separate graphs
plot(ts.y.matrix, plot.type="single", col=c(1,3))



#### Alternative for simulating y_t: draw standard normal random vectors along the way
set.seed(23456)
y.matrix <- matrix(0, nrow=K, ncol=(T.here+burn.in))
for(t in 2:(T.here+burn.in)){
  y.matrix[,t] = C + A1 %*% y.matrix[,(t-1)] + t(chol.Sigma) %*% rnorm(2)
}

plot(ts(t(y.matrix)))



## b) 
eigen(A1)  ## eigendecomposition; read up on the function using ?eigen
eigen(A1)$values  ## extract eigenvalues


## c)
# empirical mean
apply(ts.y.matrix, 2, mean)
mu.y # compare with theoretical mean. do they get closer as you increase T?






