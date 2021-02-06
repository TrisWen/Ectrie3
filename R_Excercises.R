

####################################################################
######################### Week 1 ###################################
####################################################################

####################### Problem 1 ##################################
# Question A
beta <- 0.9
sigma2 <- 2
Tt <- 100
burn_in <- 20

set.seed(123456)

u <- rnorm(Tt + burn_in, 0, sqrt(sigma2))
y <- numeric(Tt+burn_in)
y[1] = 0

for (i in 2:(burn_in+Tt)){
  y[i] = beta*y[i-1] + u[i] 
}
    
y = y[-c(1:burn_in)]
plot(y,
     type = "l")


# Question B
y1 <- matrix(y[-1], ncol = 1)
x <- matrix(y1[-Tt], ncol =1)

beta_est = solve((t(x)%*%x)) %*% (t(x)%*%y1)



#Question C
y2 <- matrix(y[-c(1,2)],ncol = 1)
x <- matrix(c(y[-c(1,Tt)], y[-c(Tt-1, Tt)]), ncol = 2)

beta_est = solve((t(x)%*%x))%*% (t(x)%*%y2)

u_hat <- y2 - x %*% beta_est
Tt_new <- length(y2)
K <- 2

u_var <- (Tt-K)^{-1} * t(u_hat)%*%u_hat

beta_cov <- as.numeric(u_var) * solve(t(x)%*%x)
se_beta <- sqrt(diag(beta_cov))

test_statistic<- beta_est[2]/se_beta[2]
2*pnorm(-abs(test_statistic))
summary(lm(y2~x -1))


#Question D
simulation <- function(t,beta0, beta1, variance, burn){
  set.seed(123456)
  error = rnorm(t + burn, 0,sqrt(variance))
  y = numeric(t+burn)
  
  for (i in 2:(t+burn)){
    y[i] = beta0 + beta1 * y[i-1] + error[i]
  }
  y = y[-c(1:burn)]
  plot(y,
       type = 'l')
  return (y)
}

y <- simulation(100,0,0.9,2,20)



####################### Problem 2 ##################################
# Question A
cov_matrix <- matrix(c(1,0.9,0.9,2), nrow=2, ncol=2)
c <- matrix(c(0.1,-0.2), ncol = 1)
A1 <- matrix(c(4/5,0,8/3,3/10), nrow=2, ncol=2)
t <- 100
burn <- 20

K <- nrow(A1)

mean_y <- solve(diag(2)-A1) %*% c

set.seed(23456)
e_matrix <- matrix(rnorm((t+burn)*K), ncol = K)
y_matrix <- matrix(0, nrow = ,(t+burn), ncol = K)

y_matrix[,1] = mean_y


chol_sigma <- chol(cov_matrix)
t(chol_sigma) %*% chol_sigma

for (i in 2:(t+burn)){
  y_matrix[i,] = c + A1 %*% y_matrix[(i-1),] + t(chol_sigma) %*% e_matrix[i,]
}

y_ts <- ts(y_matrix[-c(1:burn),], start =c(1990, 1), frequency=4)

plot(y_ts)
plot(y_ts, plot.type="single", col=c(1,3))


# Question B
eigenvalues <- eigen(A1)
eigenvalues$values



# Question C
apply(y_matrix, 2, mean)
mean_y


####################################################################
######################### Week 2 ###################################
####################################################################
