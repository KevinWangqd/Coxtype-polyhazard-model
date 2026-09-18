source("Cox_competing_correct_cluster.R")

set.seed(123)
n = 600
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

# Generate hazard and survival times
time1 <- rexp(n, rate = lambda1*lp1)
time2 <- rexp(n, rate = lambda2*lp2)
time3 <- rexp(n, rate = lambda3*lp3)

time = pmin(time1, time2, time3)
status <- rep(1,n)
# Censoring indicator (1 = event, 0 = censored)
## right censorship 
pi <- 0.1 # censorship rate

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

result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta
cov_mat <- cov_theta_profile(fit_result=result, transformed_data, sigma2 = NULL, frailty = FALSE, eps = 1/sqrt(n), tol = 1e-2, max_iter = 100)
cov_mat$se_by_group

jpeg("Bhazard_Ex1_cen1.jpeg", quality = 100, units = "in", width = 10, height = 6, res = 300)
plot(
  result$param_new[[1]]$bh_times, result$param_new[[1]]$cum_lambda,
  type = "s", lwd = 2, col = "blue",
  xlab = "Time",
  ylab = "Estimated baseline hazard",
  xlim = c(quantile(result$param_new[[1]]$bh_times,0.05),quantile(result$param_new[[1]]$bh_times,0.90)),
  ylim = c(0,quantile(result$param_new[[1]]$cum_lambda,0.90)),
  main = "Case 1 (censorship rate = 10%)"
)

lines(
  result$param_new[[1]]$bh_times, result$param_new[[2]]$cum_lambda,
  col = "orange", lwd = 2
)

lines(
  result$param_new[[1]]$bh_times, result$param_new[[3]]$cum_lambda,
  col = "red", lwd = 2
)


lines(
  result$param_new[[1]]$bh_times, 2 * result$param_new[[1]]$bh_times,
  col = "black", lwd = 3, lty = 2
)

legend("bottomright", legend = c("CF1", "CF2", "CF3"), col = c( "blue","orange","red"), lwd = 4, cex=1.5)
dev.off()


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



## Cox-Snell Residual
eval_cumhaz <- function(t, bh_times, cum_lambda) {
  idx <- findInterval(t, bh_times)
  out <- numeric(length(t))
  keep <- idx > 0
  out[keep] <- cum_lambda[idx[keep]]
  out
}

n <- length(time)
L <- ncol(result$eta_matrix)

coxsnell <- numeric(n)
for(l in 1:L){
  beta_hat <- result$param_new[[l]]$beta
  lp <- as.vector(transformed_data$X[[l]] %*% beta_hat)
  H0 <- eval_cumhaz(transformed_data$time, result$param_new[[l]]$bh_times,result$param_new[[l]]$cum_lambda)
  
  coxsnell <- coxsnell + H0 * exp(lp)
}


library(survival)
cs_fit <- survfit(Surv(coxsnell, transformed_data$status) ~ 1)


jpeg("CS_Ex1_cen3.jpeg", quality = 100, units = "in", width = 10, height = 6, res = 300)

plot(cs_fit, xlab="Cox-Snell Residual", ylab="Survival",
     main="Case 1 (censorship rate = 30%)",xlim = c(0,5.5),col="blue",lwd=1.5)

# Theoretical Exp(1): S(r) = exp(-r)
curve(exp(-x), col="black", lwd=2, add=TRUE)

legend("topright", legend=c("KM of residuals", "Exp(1)"),
       col=c("blue","black"), lwd=2)

dev.off()


