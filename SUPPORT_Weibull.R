library(survival)
library(casebase)
data(support)
attach(support)
library(riskRegression)
library(MASS)

generate_hazard_matrix <- function(model, data, time_points, model_type="cox") {
  
  n_samples <- nrow(data)
  
  hazard_matrix <- matrix(
    NA,
    nrow=n_samples,
    ncol=length(time_points)
  )
  
  
  if(model_type=="cox"){
    
    base_hazard <- basehaz(
      model,
      centered=FALSE
    )
    
    lp <- predict(
      model,
      newdata=data,
      type="lp"
    )
    
    
    for(j in seq_along(time_points)){
      
      t <- time_points[j]
      
      idx <- max(
        c(0, which(base_hazard$time <= t))
      )
      
      H0 <- if(idx==0){
        0
      } else {
        base_hazard$hazard[idx]
      }
      
      hazard_matrix[,j] <- H0 * exp(lp)
    }
    
    
  } else if(model_type=="weibull"){
    
    scale <- model$scale
    
    lp <- predict(
      model,
      newdata=data,
      type="linear"
    )
    
    
    for(j in seq_along(time_points)){
      
      t <- time_points[j]
      
      hazard_matrix[,j] <-
        (t / exp(lp))^(1/scale)
      
    }
    
  }
  
  
  colnames(hazard_matrix) <- paste0(
    "t_",
    round(time_points,1)
  )
  
  hazard_matrix
}
## coma 
data <- support[support$dzgroup == "Coma", ]
data$sex <- ifelse(data$sex=="male",1,0)
data[,c("age","scoma", "hrt","temp","adlp",
  "adlsc","sps","aps", "meanbp","wblc","resp",
  "pafi","alb","bili","sod", "crea","ph","glucose", "bun","urine")] <- scale(data[,c("age","scoma", "hrt","temp","adlp",
                                                                               "adlsc","sps","aps", "meanbp","wblc","resp",
                                                                               "pafi","alb","bili","sod", "crea","ph","glucose", "bun","urine")])

cox_full <- coxph(Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                 resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                 bun + urine + adlp + adlsc+ sps + aps, data=data)

summary(cox_full)

Weibull_full <- survreg(Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                    resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                    bun + urine + adlp + adlsc+ sps + aps, data=data)

summary(Weibull_full)


time_points <- quantile(data$d.time, probs = seq(0.01, 0.99, by = 0.07))
time_points[3] <- 4.5

hazard_cox <- generate_hazard_matrix(cox_full, data, time_points
                                     , "cox")

time_AUC <- Score(
  object   = list(Cox = hazard_cox),
  formula  = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$d.time,
    status = data$death
  ),
  times    = time_points,
  metrics  = "AUC",
  cens.model = "km"
)

mean(time_AUC$AUC$score$AUC)

# Calculate Brier Score
bs_res <- Score(
  object      = list(Cox = 1 - exp(-hazard_cox)),
  formula     = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$d.time,
    status = data$death
  ),
  times    = time_points,
  metrics     = "Brier",
  cens.model  = "km"
)

mean(bs_res$Brier$score$Brier[which(bs_res$Brier$score$model=="Cox")])



hazard_W <- generate_hazard_matrix(Weibull_full, data, time_points
                                     , "weibull")

time_AUC <- Score(
  object   = list(Weibull = hazard_W),
  formula  = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$d.time,
    status = data$death
  ),
  times    = time_points,
  metrics  = "AUC",
  cens.model = "km"
)

mean(time_AUC$AUC$score$AUC)

# Calculate Brier Score
bs_res <- Score(
  object      = list(Weibull = 1 - exp(-hazard_W)),
  formula     = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$d.time,
    status = data$death
  ),
  times    = time_points,
  metrics     = "Brier",
  cens.model  = "km"
)

mean(bs_res$Brier$score$Brier[which(bs_res$Brier$score$model=="Weibull")])


















source("Functions_new.txt")

Y <- log(data$d.time)
delta <- data$death

load("X_can.RData")
load("label_can.RData")
transformed_data = data_transform(Y, delta, X_can, label_can)


init_param = param_init(X_can, label_can)

result <- EM_fit(transformed_data, init_param, lambda1 = 0, lambda2 = 0, maxit = 10^3, tolerance = 1e-6)

#result$param_new

summary_result <- summary_EM_fit(result, transformed_data)

# Print the summary for each group
summary_result
#scoma age sps adlsc sod alb
#aps hrt adlsc bili urine bun meanbp
#age sps temp pafi bili urine ph


#expected_time = sur_time_all(result, transformed_data,method = "MR", M = 200, n_samples = 1000)
expected_time = sur_time_all(result, transformed_data,method = "IS_exp", M = 100, n_samples = 10^4)

surv_obj <- Surv(time = transformed_data$Y, event = transformed_data$delta)

# Calculate the concordance index
c_index <- concordance(Surv(transformed_data$Y, transformed_data$delta) ~ expected_time)

# Print the result
print(c_index$concordance)

library(timeROC)
library(survival)




calculate_total_hazard_weibull <- function(time_points, result, transformed_data) {
  
  L <- length(result$param_new)
  n_samples <- nrow(transformed_data$X[[1]])
  
  hazard_matrix <- sapply(time_points, function(t) {
    
    total_hazard <- numeric(n_samples)
    
    for (l in 1:L) {
      
      sigma <- result$param_new[[l]][[1]]
      coef  <- result$param_new[[l]][[2]]
      
      alpha <- coef[1]
      beta  <- coef[-1]
      
      X_beta <- as.numeric(transformed_data$X[[l]] %*% beta)
      
      linpred <- alpha + X_beta
      
      # Weibull hazard contribution
      #h_l <- (1 / (sigma * t)) *
       # (t / exp(linpred))^(1 / sigma)
      
      #total_hazard <- total_hazard + h_l
      H_l <- (t / exp(linpred))^(1/sigma)
      
      total_hazard <- total_hazard + H_l
    }

    
    total_hazard
  })
  
  colnames(hazard_matrix) <- paste0("t_", round(time_points, 2))
  
  S_matrix <- exp(-hazard_matrix)
  
  return(list(
    hazard_matrix = hazard_matrix,
    S_matrix = S_matrix
  ))
}






compute_iAUC_weibull <- function(transformed_data, result, timepoint) {
  
  time_points <- timepoint
  total <- calculate_total_hazard_weibull(time_points, result, transformed_data)
  
  time_AUC <- riskRegression::Score(
    object = list(competingWeibull = total$hazard_matrix),
    formula = Surv(time, status) ~ 1,
    data = data.frame(
      time = exp(transformed_data$Y),
      status = transformed_data$delta
    ),
    times = time_points,
    metrics = "AUC",
    cens.model = "km"
  )
  
  mean(time_AUC$AUC$score$AUC,na.rm = T)
}

compute_iBS_weibull <- function(transformed_data, result, timepoint) {
  
  time_points <- timepoint
  total <- calculate_total_hazard_weibull(time_points, result, transformed_data)
  
  bs_res <- riskRegression::Score(
    object = list(competingWeibull = 1 - total$S_matrix),
    formula = Surv(time, status) ~ 1,
    data = data.frame(
      time = exp(transformed_data$Y),
      status = transformed_data$delta
    ),
    times = time_points,
    metrics = "Brier",
    cens.model = "km"
  )
  
  mean(bs_res$Brier$score$Brier[
    bs_res$Brier$score$model == "competingWeibull"
  ],na.rm = T)
}



compute_iAUC_weibull( transformed_data, result,time_points )
compute_iBS_weibull( transformed_data, result,time_points )



## Diagnostic plots
L <- length(result$param_new)

sigma <- numeric(L)
alpha <- numeric(L)
beta  <- vector("list", L)

for (l in 1:L) {
  sigma[l] <- result$param_new[[l]][[1]]
  alpha[l] <- result$param_new[[l]][[2]][1]
  beta[[l]] <- result$param_new[[l]][[2]][-1]
}
compute_survival <- function(t, sigma, alpha, beta, X_list) {
  n <- nrow(X_list[[1]])
  L <- length(sigma)
  
  H <- rep(0, n)  # cumulative hazard
  
  for (l in 1:L) {
    H <- H + (t / exp(alpha[l] + X_list[[l]] %*% beta[[l]]))^(1 / sigma[l])
  }
  S <- exp(-H)
  return(S)
}

time_grid <- seq(0.01, exp(max(transformed_data$Y)), length.out = 100)

S_mat <- sapply(time_grid, function(t) {
  compute_survival(t, sigma, alpha, beta, transformed_data$X)
})


library(survival)




compute_cs_residual <- function(time, sigma, alpha, beta, X_list) {
  n <- length(time)
  L <- length(sigma)
  
  H <- rep(0, n)
  
  for (l in 1:L) {
    eta_l <- alpha[l] + X_list[[l]] %*% beta[[l]]
    lambda_l <- exp(eta_l)
    
    H <- H + (time / lambda_l)^(1 / sigma[l])
  }
  
  return(H)
}



r_cs <- compute_cs_residual(
  time = exp(transformed_data$Y),
  sigma = sigma,
  alpha = alpha,
  beta = beta,
  X_list = transformed_data$X
)

cs_fit <- survfit(Surv(r_cs, transformed_data$delta) ~ 1)


jpeg("diag_sepsis2.jpeg", quality = 100, units = "in", width = 8, height = 6, res = 300)

plot(cs_fit, col="blue", xlab="Cox-Snell Residual", ylab="Survival Probability")

# Theoretical Exp(1): S(r) = exp(-r)
curve(exp(-x), col="black", lwd=2, add=TRUE)

legend("topright", legend=c("KM of Cox-Snell residuals", "Exp(1)"),
       col=c("blue","black"), lwd=2)
dev.off()


