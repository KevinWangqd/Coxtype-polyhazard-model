## Frailty model
source("Cox_competing_correct.R")
set.seed(123)
n = 1500
gamma = 0.5
u <- rgamma(n, shape = 1/gamma, scale = gamma)  # mean=1, var=gamma

X_init <- matrix(rnorm(n * 6), ncol = 6)

X1 <- X_init[,c(1,2,3)]
X2 <- X_init[,c(3,4)]
X3 <- X_init[,c(3,5,6)]

X <- cbind(X1,X2,X3)

# True beta values per group

true_beta <- list(c(2, 2, -1.5), c(-1.5,2), c(-2, -2, 2.5))

# Compute LPs under each group's model
lp1 <- u*exp(X1 %*% true_beta[[1]])  # n x 1
lp2 <- u*exp(X2 %*% true_beta[[2]])  # n x 1
lp3 <- u*exp(X3 %*% true_beta[[3]])  # n x 1

lambda1 <- 1
lambda2 <- 1
lambda3 <- 1


lambda11 <- 2
lambda22 <- 0.5
lambda33 <- 1.5

# Generate hazard and survival times
time1 <- rexp(n, rate = lambda1*lp1)
time2 <- rexp(n, rate = lambda2*lp2)
time3 <- rexp(n, rate = lambda3*lp3)

time11 <- rexp(n, rate = lambda11*lp1)
time22 <- rexp(n, rate = lambda22*lp2)
time33 <- rexp(n, rate = lambda33*lp3)

time = pmin(time1, time2, time3)
summary(time)

cut=0.3

time <- ifelse(time<cut,time, cut + pmin(time11,time22,time33))


status <- rep(1,n)
# Censoring indicator (1 = event, 0 = censored)
## right censorship 
pi <- 0.3 # censorship rate

f <- function(lambda){
  mean(punif(time, min = 0, max =lambda)-pi)
}

lambda = uniroot(f, interval = c(0.0001, 10^6))$root
C <- runif(n, 0, lambda)
status = rep(1,n)

for (i in 1:n) {
  if(time[i]>C[i]){
    time[i] = C[i]
    status[i] = 0
  }
}

sum(status)/n

label_true = list(1,1,1,2,2,3,3,3)
# Print transformed data to check the result
transformed_data = data_transform(time, status, X, label_true)

result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-4, frailty = T)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta


cov <- cov_theta_profile(fit_result = result, transformed_data,
                              sigma2 = result$sigma2, frailty = T, 
                              eps = 1/sqrt(n), tol = 5e-3,
                              fast = T, max_iter = 100)


cov$se_by_group
cov$se_sigma2



## CCR
true_cause <- apply(cbind(time1,time2,time3), 1, which.min)
true_cause <- true_cause[order(time)]
estimated_cause <- apply(result$eta_matrix, 1, which.max)


L <- ncol(result$eta_matrix)

CCR_by_cause <- sapply(1:L, function(l) {
  idx <- transformed_data$status == 1 & true_cause == l
  mean(estimated_cause[idx] == l)
})

CCR_by_cause

