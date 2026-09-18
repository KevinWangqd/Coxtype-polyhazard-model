source("Cox_competing_correct_cluster.R")


## Frailty model
set.seed(123)
n = 600

gamma = 1
u <- rgamma(n, shape = 1/gamma, scale = gamma)  # mean=1, var=gamma

X_init <- matrix(rnorm(n * 3), ncol = 3)

X1 <- X_init[,c(1)]
X2 <- X_init[,c(2)]
X3 <- X_init[,c(3)]

X <- cbind(X1,X2,X3)

# True beta values per group
true_beta <- list(c(2), c(1.5), c(-2))

# Compute LPs under each group's model
lp1 <- exp(X1 * true_beta[[1]])  # n x 1
lp2 <- exp(X2 * true_beta[[2]])  # n x 1
lp3 <- exp(X3 * true_beta[[3]])  # n x 1

lambda1 <- 2
lambda2 <- 2
lambda3 <- 2
# Cause-specific latent times
# --------------------
time1 <- rexp(n, rate = u * lambda1 * lp1)
time2 <- rexp(n, rate = u * lambda2 * lp2)
time3 <- rexp(n, rate = u * lambda3 * lp3)

# --------------------
# Observed competing risks outcome
# --------------------
time   <- pmin(time1, time2, time3)
status <- rep(1, n)

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

label_true = list(1,2,3)
# Print transformed data to check the result
transformed_data = data_transform(time, status, X, label_true)

result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 5e-3, frailty = T)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta
#result$b_hat


cov_mat <- cov_theta_profile(fit_result=result, transformed_data, sigma2 = NULL, frailty = T, eps = 1/sqrt(n), tol = 1e-2, max_iter = 100)
cov_mat$se_by_group

## CCR
true_cause <- apply(cbind(time1,time2, time3), 1, which.min)
true_cause <- true_cause[order(time)]

eta_matrix <- result$eta_matrix

estimated_cause <- apply(eta_matrix, 1, which.max)

estimated_cause <- estimated_cause[which(transformed_data$status == 1)]
true_cause <- true_cause[which(transformed_data$status == 1)]



L <- ncol(result$eta_matrix)

CCR_by_cause <- sapply(1:L, function(l) {
  idx <-  estimated_cause == l
  mean(true_cause[idx] == l)
})

CCR_by_cause
