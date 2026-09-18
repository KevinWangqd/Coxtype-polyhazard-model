source("Cox_competing_correct_cluster_LASSO.R")
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

result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F, penalty = 0)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta


#cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)

#cov$se_by_group



##high dim
X_omit <- matrix(rnorm(n * 194), ncol = 194)
X_select <- cbind(X_init, X_omit)



label_true = c(list(1,1,c(1,2,3),2,3,3),sample(rep(1:3, length.out=194)))
# Print transformed data to check the result
transformed_data = data_transform(time, status, X_select, label_true)

cvfit <- glmnet::cv.glmnet(x = X_select, y = Surv(time, status),family = "cox", nfolds=5, alpha = 1)
cvfit$lambda.min

result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F, penalty = 0.03)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)
cov$se_by_group

## Try1
res_high <- high_select(time = time , status = status, X = X_select, n_group = 3, size_group = 3, 
                        init_iauc = 0.6, in_iauc = 0.02, tar_P = 30, n_draw = 100)
#--- Iteration 6 ---
#Current candidate pool size: 47
#Current candidate pool: 1, 2, 3, 4, 5, 6, 7, 13, 14, 19, 21, 22, 30, 32, 41, 55, 61, 66, 67, 70, 88, 97, 102, 108, 113, 118, 121, 122, 126, 128, 130, 136, 139, 142, 144, 151, 153, 165, 169, 173, 180, 181, 184, 186, 187, 193, 197
#Current iAUC threshold: 0.7

X_select <- X_select[,c(1, 2, 3, 4, 5, 6, 7, 13, 14, 19, 21, 22, 30, 32, 41, 55, 61, 66, 67, 70, 88, 97, 102, 108, 113, 118, 121, 122, 126, 128, 130, 136, 139, 142, 144, 151, 153, 165, 169, 173, 180, 181, 184, 186, 187, 193, 197)]
res_high <- high_select(time = time , status = status, X = X_select, n_group = 3, size_group = 3, 
                        init_iauc = 0.72, in_iauc = 0.02, tar_P = 30, n_draw = 100)

#Final candidate pool size: 28
#Final candidate pool: 1  2  3  4  5  6  7  8 12 15 16 17 20 22 24 26 27 28 29 30 34 37 39 40 42 43 44 47
#Final iAUC threshold: 0.78
#Indices in original X_select: 1   2   3   4   5   6   7  13  22  41  55  61  70  97 108 118 121 122 126 128 142 153 169 173 181 184 186 197

length(res_high$selected_indices)
res_high$X_selected

res_forward <- forward_select(time, status, X_candidate=res_high$X_selected, max_add = 10, lambda1 = log(n), lambda2 = log(n))
res_forward$label
res_forward$AIC
res_forward$X

#X3 X1 X2
#X3 X6 X5 X16
#X3 X4 X7

transformed_data = data_transform(time, status, res_forward$X, res_forward$label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)
cov$se_by_group

compute_iAUC(transformed_data,result = result)
compute_iBS(transformed_data,result = result)


### Try2
X_select <- cbind(X_init, X_omit)
res_high <- high_select(time = time , status = status, X = X_select, n_group = 4, size_group = 4, 
                        init_iauc = 0.6, in_iauc = 0.02, tar_P = 30, n_draw = 100)
sort(res_high$selected_indices)
res_high$X_selected

#Final candidate pool size: 25
#Final candidate pool: 1   2   3   5   6  15  18  28  36  51  59  65  77  96 106 110 126 130 155 158 159 174 178 179 181
#Final iAUC threshold: 0.82

res_forward <- forward_select(time, status, X_candidate=res_high$X_selected, max_add = 10, lambda1 = log(n), lambda2 = log(n))
res_forward$label
res_forward$AIC
res_forward$X[1:2,]

#X3 
#X3 X6 X5 X158 
#X3 X2 X1

transformed_data = data_transform(time, status, res_forward$X, res_forward$label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)
cov$se_by_group

compute_iAUC(transformed_data,result = result)
compute_iBS(transformed_data,result = result)



### Try3
# By iAUC threshold:0.82
X_select <- X_select[,c(1, 2, 3, 4, 5, 6, 7, 12, 13, 19, 20, 22, 25, 26, 31, 33, 53, 55, 57, 72, 75, 88, 91, 93, 101, 105, 107, 108, 112, 123, 135, 137, 142, 146, 149, 150, 151, 153, 157, 160, 192, 193, 196)]
res_high <- high_select(time = time , status = status, X = X_select, n_group = 4, size_group = 4, 
                        init_iauc = 0.84, in_iauc = 0.02, tar_P = 30, n_draw = 300)
sort(res_high$selected_indices)
res_high$X_selected

#Final candidate pool size: 23
#Final candidate pool:  1   2   3   4   5   6   7  25  26  31  88  91  93 107 108 112 123 142 146 149 150 160 196
#Final iAUC threshold: 0.86
#1  2  3  4  5  6  7 13 14 15 22 23 24 27 28 29 30 33 34 35 36 40 43

res_forward <- forward_select(time, status, X_candidate=res_high$X_selected, max_add = 10, lambda1 = log(n), lambda2 = log(n))
res_forward$label
res_forward$X[1:2,]
#X3 X1 X2
#X3 X6 X5 
#X3 X4 X7 X8


transformed_data = data_transform(time, status, res_forward$X, res_forward$label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(n), tol = 5e-3, fast = T, max_iter = 100)
cov$se_by_group

compute_iAUC(transformed_data,result = result)
compute_iBS(transformed_data,result = result)
