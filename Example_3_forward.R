source("Cox_competing_correct.R")
source("Cox_competing_algorithm.R")

set.seed(123)
n = 1500
X_init <- matrix(rnorm(n * 6), ncol = 6)

X1 <- X_init[,c(1,2,3)]
X2 <- X_init[,c(3,4)]
X3 <- X_init[,c(3,5,6)]

X <- cbind(X1,X2,X3)

# True beta values per group

true_beta <- list(c(2, 2, -1.5), c(-1.5,2), c(-2, -2, 2.5))

# Compute LPs under each group's model
lp1 <- exp(X1 %*% true_beta[[1]])  # n x 1
lp2 <- exp(X2 %*% true_beta[[2]])  # n x 1
lp3 <- exp(X3 %*% true_beta[[3]])  # n x 1

lambda1 <- 1
lambda2 <- 1
lambda3 <- 1


# Generate hazard and survival times
time1 <- rexp(n, rate = lambda1*lp1)
time2 <- rexp(n, rate = lambda2*lp2)
time3 <- rexp(n, rate = lambda3*lp3)

time = pmin(time1, time2, time3)

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

result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta


#cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)

#cov$se_by_group

# Other variables
X_omit <- matrix(rnorm(n *4), ncol = 4)

X_select <- cbind(X_init, X_omit)


res_forward <- forward_select(time, status, X_candidate=X_select, max_add = 10, lambda1 = log(n), lambda2 = log(n))
res_forward$label
res_forward$AIC
res_forward$X

transformed_data = data_transform(time, status, res_forward$X, res_forward$label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)
cov$se_by_group

compute_iAUC(transformed_data,result = result)
compute_iBS(transformed_data,result = result)


res_forward3 <- forward_select(time, status, X_candidate=X_select, max_add = 10, lambda1 = 2, lambda2 = 2)
res_forward3$label
res_forward3$AIC
res_forward3$X

#X3 X6 X4 X5 X1 X2 X9 X7

#X3 X1 X2
#X3 X6 X5 X7
#X3 X4 X9

transformed_data = data_transform(time, status, res_forward3$X, res_forward3$label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)
cov$se_by_group

compute_iAUC(transformed_data,result = result)
compute_iBS(transformed_data,result = result)
