library(survival)
library(casebase)
data(support)
attach(support)
source("Cox_competing_correct_cluster.R")
source("Cox_competing_algorithm.R")
library(riskRegression)
library(MASS)
generate_hazard_matrix <- function(model, data, time_points, model_type = "cox") {
  
  n_samples <- nrow(data)
  
  hazard_matrix <- matrix(
    0,
    nrow = n_samples,
    ncol = length(time_points)
  )
  
  
  if (model_type == "cox") {
    
    ############################################################
    # Cox PH model
    ############################################################
    
    lp <- predict(
      model,
      newdata = data,
      type = "lp"
    )
    
    
    # Standard Cox: direct baseline hazard
    if (is.null(model$frail)) {
      
      base_hazard <- basehaz(
        model,
        centered = FALSE
      )
      
    } else {
      
      ##########################################################
      # Frailty Cox: reconstruct Breslow baseline hazard
      ##########################################################
      
      train_lp <- predict(
        model,
        type = "lp"
      )
      
      fit_offset <- coxph(
        Surv(model$y[,1], model$y[,2]) ~ offset(train_lp),
        data = model$model
      )
      
      base_hazard <- basehaz(
        fit_offset,
        centered = FALSE
      )
      
    }
    
    
    ############################################################
    # Add frailty contribution if available
    ############################################################
    
    if (!is.null(model$frail)) {
      
      b_hat <- model$frail
      
      # assuming data corresponds to subjects with estimated frailty
      lp <- lp + b_hat
      
    }
    
    
    ############################################################
    # Cumulative hazard
    ############################################################
    
    for (j in seq_along(time_points)) {
      
      t <- time_points[j]
      
      idx <- max(
        c(0, which(base_hazard$time <= t))
      )
      
      H0_t <- if(idx == 0){
        0
      } else {
        base_hazard$hazard[idx]
      }
      
      hazard_matrix[,j] <- H0_t * exp(lp)
    }
    
  }
  
  
  colnames(hazard_matrix) <- paste0(
    "t_",
    round(time_points,1)
  )
  
  return(hazard_matrix)
}

generate_hazard_matrix <- function(model, data, time_points, model_type = "cox") {
  
  n_samples <- nrow(data)
  hazard_matrix <- matrix(NA, nrow = n_samples, ncol = length(time_points))
  
  if (model_type == "cox") {
    # Cox PH model
    base_hazard <- basehaz(model, centered = FALSE)
    lp <- predict(model, type = "lp")
    
    for (j in 1:length(time_points)) {
      t <- time_points[j]
      idx <- max(which(base_hazard$time <= t))
      baseline_cum_hazard <- ifelse(length(idx) > 0, base_hazard$hazard[idx], 0)
      hazard_matrix[, j] <- baseline_cum_hazard * exp(lp)
    }
    
  } else if (model_type == "weibull") {
    # Weibull AFT model
    scale <- model$scale
    lp <- predict(model, newdata = data, type = "linear")
    
    for (j in 1:length(time_points)) {
      t <- time_points[j]
      z <- (log(t) - lp) / scale
      hazard_matrix[, j] <- (1/(scale * t)) * exp(z)^(1/scale) / gamma(1/scale)
    }
  }
  
  colnames(hazard_matrix) <- paste0("t_", round(time_points, 1))
  return(hazard_matrix)
}


## coma 
data <- support[support$dzgroup == "Coma", ]
data$sex <- ifelse(data$sex=="male",1,0)

data[,c("age","scoma", "hrt","temp","adlp",
  "adlsc","sps","aps", "meanbp","wblc","resp",
  "pafi","alb","bili","sod", "crea","ph","glucose", "bun","urine")] <- scale(data[,c("age","scoma", "hrt","temp","adlp",
                                                                               "adlsc","sps","aps", "meanbp","wblc","resp",
                                                                               "pafi","alb","bili","sod", "crea","ph","glucose", "bun","urine")])
summary(data$slos)

#data$cluster <- as.factor(ifelse(data$race=="white",1,ifelse(data$race=="black",2,3)))
#data$cluster <- as.factor(ifelse(data$num.co>1,0,1))

#data$cluster <- as.factor(ifelse(data$slos<6,1,ifelse(data$slos<17,2,3)))



cox_full <- coxph(Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                 resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                 bun + urine + adlp + adlsc+ sps + aps , data=data)

summary(cox_full)


cox_aic_result <- cox_full
## COMA: iAUC 0.809; iBS 0.141

cox_aic_result <- stepAIC(cox_full, upper = Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                            resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                            bun + urine + adlp + adlsc+ sps + aps,
                          direction = "both",
                          trace = FALSE, k=log(596))
## iAUC 0.798; iBS 0.143
cox_aic_result$formula

cox_aic_result <- stepAIC(cox_full, upper = Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                            resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                            bun + urine + adlp + adlsc+ sps + aps,
                          direction = "both",
                          trace = FALSE, k=2)
summary(cox_aic_result)
## iAUC 0.807; iBS 0.142


timepoints <- quantile(data$d.time, probs = seq(0.01, 0.99, by = 0.07))
timepoints[3] <- 4.5

hazard_cox <- generate_hazard_matrix(cox_aic_result, data, timepoints, "cox")

time_AUC <- Score(
  object   = list(Cox = hazard_cox),
  formula  = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$d.time,
    status = data$death
  ),
  times    = timepoints,
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
  times    = timepoints,
  metrics     = "Brier",
  cens.model  = "km"
)

mean(bs_res$Brier$score$Brier[which(bs_res$Brier$score$model=="Cox")])



## Forward selection
#X_candidate = cbind(data$sex,data[,c("age","scoma", "hrt","temp","adlp",
#                                     "adlsc","sps","aps", "meanbp","wblc","resp",
#                                     "pafi","alb","bili","sod", "crea","ph","glucose", "bun","urine")])

#res_forward <- forward_select(data$d.time, data$death, X_candidate=X_candidate, max_add = 25, lambda1 = 6, lambda2 = log(596))
#res_forward$label
#res_forward$X[1:3,]

#X_can <- res_forward$X[,1:15]
#scoma aps age sps hrt
#temp pafi adlsc bili urine
# bun ph sod meanbp alb

#label_can <- res_forward$label[1:15]
#label_can[[8]] <- c(1,2)
#label_can[[7]] <- c(3)
#label_can[[10]] <- c(2,3)
#label_can[[11]] <- c(2)
#label_can[[9]] <- c(2,3)


#1 3 4 8 13 15
#2 5 8 9 10 11 14
#3 4 6 7 9 10 12

#scoma age sps adlsc sod alb
#aps hrt adlsc bili urine bun meanbp
#age sps temp pafi bili urine ph


load("X_can.RData")
load("label_can.RData")


transformed_data = data_transform(data$d.time, data$death, X_can, label_can,cluster = data$cluster)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

result$b_hat
#result$var_b
result$sigma2


cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = T, eps = 1/sqrt(length(data$d.time)), tol = 1e-3, fast = T, max_iter = 100)
cov$se_by_group

cov$se_sigma2

compute_iAUC(transformed_data,result = result, timepoints)
compute_iBS(transformed_data,result = result, timepoints)

#iAUC 0.830; iBS 0.136













result$eta_matrix[transformed_data$status==1,]

eta_matrix <- result$eta_matrix
# Subset eta matrix for event patients
eta_event <- eta_matrix[transformed_data$status == 1, ]

# Number of competing groups
L <- ncol(eta_event)

# Function to compute normalized entropy for one patient
compute_H <- function(eta_i) {
  H_i <- -sum(eta_i * log(eta_i)) / log(L)
  return(H_i)
}

# Apply to all selected patients
H_values <- apply(eta_event, 1, compute_H)

# Optional: attach to data frame
result <- data.frame(
  H = H_values
)

# Quick summary
summary(H_values)

H_values

class <- ifelse(H_values < 0.3333, "single",
                ifelse(H_values < 0.666, "double", "triple"))

max_col <- max.col(eta_matrix, ties.method = "first")
max_col

top2_cols <- t(apply(eta_matrix, 1, function(x) order(x, decreasing = TRUE)[1:2]))


patient_strata <- as.data.frame(cbind(time=transformed_data$time[transformed_data$status==1],
                                      result$eta_matrix[transformed_data$status==1,],
                                      class,group = ifelse(class == "single", 
                                                            max_col, 
                                                            ifelse(class == "double", 
                                                                   apply(eta_matrix, 1, function(x) paste(order(x, decreasing = TRUE)[1:2], collapse = ",")), 
                                                                   "1,2,3")) ))

## Cox-Snell Residual
eval_cumhaz <- function(t, bh_times, cum_lambda) {
  idx <- findInterval(t, bh_times)
  out <- numeric(length(t))
  keep <- idx > 0
  out[keep] <- cum_lambda[idx[keep]]
  out
}

n <- length(data$d.time)
L <- ncol(result$eta_matrix)

coxsnell <- numeric(n)
for(l in 1:L){
  beta_hat <- result$param_new[[l]]$beta
  lp <- as.vector(transformed_data$X[[l]] %*% beta_hat)
  H0 <- eval_cumhaz(transformed_data$time, result$param_new[[l]]$bh_times,result$param_new[[l]]$cum_lambda)
  
  coxsnell <- coxsnell + H0 * exp(lp)
}

cs_fit_cox <- survfit(Surv(coxsnell, transformed_data$status) ~ 1)


library(survival)
load("COMA_Weibull.Rdata")
cs_fit_weibull <- cs_fit






jpeg("CS_comp_COMA.jpeg", quality = 100, units = "in", width = 10, height = 8, res = 300)

plot(cs_fit_cox, xlab="Cox-Snell Residual", ylab="Survival",
     main="Competing Cox model",xlim = c(0,4))

# Theoretical Exp(1): S(r) = exp(-r)
curve(exp(-x), col="red", lwd=2, add=TRUE)

legend("topright", legend=c("KM of residuals", "Exp(1)"),
       col=c("black","red"), lwd=2)

dev.off()



jpeg("CS_comp_COMA.jpeg", quality = 100, units = "in", width = 10, height = 8, res = 300)
plot(cs_fit_cox,
     xlab = "Cox-Snell Residual",
     ylab = "Survival",
     main = "Cox-Snell Residual Diagnostics",
     xlim = c(0.1, 5),
     col = "blue",
     lwd = 2)

lines(cs_fit_weibull,
      col = "red",
      lwd = 2)

# Theoretical exp(-r) line
curve(exp(-x),
      from = 0, to = 5,
      add = TRUE,
      col="black",
      lwd = 2)

legend("topright",
       legend = c("KM of Cox-type Polyhazard residuals",
                  "KM of Poly-Weibull residuals",
                  "Exp(1)"),
       col = c("blue", "red", "black"),
       lwd = c(2, 2, 2))
dev.off()




## One choice
cox_full <- coxph(Surv(d.time,death)~ age+ sex+ scoma + hrt+
                    temp  + adlp + adlsc+ sps + aps, data=data)

X = cbind(data$age,data$sex,data$scoma, data$hrt,data$temp,data$adlp, data$adlsc,data$sps,data$aps)

label = list(c(1,2,3), c(2),c(3), c(2), c(1),c(1), c(3), c(1,3),c(2))

transformed_data = data_transform(data$d.time, data$death, X, label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$eta_matrix
#result$param_new[[1]]$beta
#result$param_new[[2]]$beta
#result$param_new[[3]]$beta

#cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(length(data$time_day)), tol = 1e-3, fast = F, max_iter = 100)
#cov$se_by_group

compute_iAUC(transformed_data,result = result, timepoints)
compute_iBS(transformed_data,result = result, timepoints)
## iAUC=0.818; iBS =0.139







