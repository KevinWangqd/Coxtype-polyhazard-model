source("Cox_competing_correct.R")
source("Cox_competing_algorithm.R")
load("GSE72094.RData")
library(riskRegression)
library(MASS)
#load("final_gene_compare.RData")
#Function for Cox model
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


## Standardization
#colnames(data)
gene_names <- c("merck-NM_001089_at", "merck-NM_003034_at", "merck2-AI858819_at", "merck-AA129758_a_at",
  "merck-AK025097_at", "merck-BC023330_at", "merck-BG772810_s_at", "merck-ENST00000380176_a_at", "merck-NM_001033567_at",
  "merck-NM_153838_a_at", "merck-AF499000_a_at", "merck-AK057093_at", "merck-BM670173_s_at", "merck-BX107093_at",
  "merck-CR610205_at", "merck-ENST00000327625_a_at", "merck-hCT1653457.1_at", "merck-HSG00446072_x_at",
  "merck-NM_000165_at", "merck-NM_000782_at", "merck-NM_001007229_at",
  "merck-NM_003901_a_at", "merck-NM_017888_at", "merck-NM_020925_s_at", "merck-NM_033125_at", "merck2-AB014565_at",
  "merck2-AJ417818_at", "merck2-BG739747_s_at", "merck2-CB306793_at", "merck-AB209100_x_at"
)

X_candidate = scale(data[,gene_names])
X_candidate <- as.matrix(X_candidate)
X_candidate=cbind(scale(data$age),X_candidate)
colnames(X_candidate) <- c("age", "ABCA3",  "ST8SIA1","merck2-AI858819_at",  "EIF4E3",  "AK025097",
                           "BC023330",  "EIF6",  "PDLIM5", "RHOT1",  "GPR115",
                           "CD302 LY75-CD302",  "merck-AK057093_at",  "merck-BM670173_s_at",  
                           "AC096670",  "PPIE",
                           "SLCO4C1","merck-hCT1653457.1_at", "SNORD12", "GJA1" ,  "CYP24A1",
                           "SPOP",  "SGPL1", "ACSM5",  "CACHD1", "SLC22A16",
                           "RAB11FIP3",   "MUC15",  "merck2-BG739747_s_at",  "INPP5B",   "RPL21P19")
                        
X_candidate <- X_candidate[,c(1,2,3,5,
                              8:11,16:17,19:20, 22:28,30,31)]  
data <- data[,1:6]

GSE72094 <- read.csv("final_data.csv")
GSE72094$CCNA2 <- GSE72094$merck.NM_001237_a_at
GSE72094$AURKA <- GSE72094$merck2.BE856617_at
GSE72094$AURKB <- GSE72094$merck.NM_004217_at
GSE72094$FEN1 <- GSE72094$merck2.XM_937756_a_at
GSE72094$CCND3 <- GSE72094$merck2.BQ669293_at
GSE72094$NCALD <- GSE72094$merck2.AF251061_at
GSE72094$MACF1 <- GSE72094$merck2.AB029290_at
GSE72094$LRC4 <- GSE72094$merck2.NM_007360_at
GSE72094$NLRC4 <- GSE72094$merck.NM_021209_s_at
GSE72094$PLEKHN1 <- GSE72094$merck.NM_032129_at
GSE72094$RASIP1 <- GSE72094$merck.NM_017805_at
GSE72094$SPP1 <- GSE72094$merck2.CA447290_at
GSE72094$GPT2 <- GSE72094$merck2.BX099266_at
#GSE72094$SGPL1 <- GSE72094$merck.ENST00000299297_at
GSE72094$PCOLCE2 <- GSE72094$merck2.AK223633_at

GSE72094[,64:77] <- scale(GSE72094[,64:77])
X_candidate <- cbind(X_candidate,GSE72094[,64:77])

## Cox competing selection
n <- nrow(data)
res_forward <- forward_select(data$time_day, data$status, X_candidate=X_candidate, max_add = 25, lambda1 = log(n), lambda2 = log(n))
res_forward$label
res_forward$X[1:3,]


#BIC_new
#1 5 6 11 14    
#2 4 8 9 10 13    
#3 4 7 12 
#ABCA3 SLC22A16 SLCO4C1 NCALD EIF4E3
#RPL21P19 GJA1 FEN1 ACSM5 SNORD12 INPP5B  
#SGPL1 GJA1 RAB11FIP3 GPT2

#AIC_new
#1 5 6 11 14 15 17   
#2 4 8 9 10 13 16 20   
#3 4 7 12 18 19

#ABCA3 SLC22A16 SLCO4C1 NCALD EIF4E3 PDLIM5 MUC15
#RPL21P19 GJA1 FEN1 ACSM5 SNORD12  INPP5B NLRC4  RHOT1          
#SGPL1 GJA1 RAB11FIP3  GPT2 MACF1 EIF6

transformed_data = data_transform(data$time_day, data$status, res_forward$X, res_forward$label)
result <- Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-6, frailty = F)
result$param_new[[1]]$beta
result$param_new[[2]]$beta
result$param_new[[3]]$beta

cov <- cov_theta_profile(fit_result=result, transformed_data, frailty = F, eps = 1/sqrt(length(data$time_day)), tol = 1e-3, fast = F, max_iter = 100)
cov$se_by_group

compute_iAUC(transformed_data,result = result)
compute_iBS(transformed_data,result = result)

## Cox selection
cov_names <- paste0("`", colnames(X_candidate), "`")

formula_str <- paste("Surv(data$time_day, data$status) ~", paste(cov_names, collapse = " + "))


# Fit Cox model
cox_full <- coxph(as.formula(formula_str), data = as.data.frame(X_candidate))
summary(cox_full)

cox_aic_result <- cox_full
cox_aic_result <- stepAIC(cox_full, upper = as.formula(formula_str),
                          direction = "both",
                          trace = FALSE, k=log(n)
)
cox_aic_result <- stepAIC(cox_full, upper = as.formula(formula_str),
                          direction = "both",
                          trace = FALSE, k=2
)
summary(cox_aic_result)



hazard_cox <- generate_hazard_matrix(cox_aic_result, data, quantile(data$time_day, probs = seq(0.01, 0.99, by = 0.05))
                                     , "cox")

time_AUC <- Score(
  object   = list(Cox = hazard_cox),
  formula  = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$time_day,
    status = data$status
  ),
  times    = quantile(data$time_day, probs = seq(0.01, 0.99, by = 0.05)),
  metrics  = "AUC",
  cens.model = "km"
)

mean(time_AUC$AUC$score$AUC)

# Calculate Brier Score
bs_res <- Score(
  object      = list(Cox = 1 - exp(-hazard_cox)),
  formula     = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = data$time_day,
    status = data$status
  ),
  times    = quantile(data$time_day, probs = seq(0.01, 0.99, by = 0.05)),
  metrics     = "Brier",
  cens.model  = "km"
)

mean(bs_res$Brier$score$Brier[which(bs_res$Brier$score$model=="Cox")])

#


## Cox-Snell Residual
eval_cumhaz <- function(t, bh_times, cum_lambda) {
  idx <- findInterval(t, bh_times)
  out <- numeric(length(t))
  keep <- idx > 0
  out[keep] <- cum_lambda[idx[keep]]
  out
}

n <- length(data$time_day)
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


jpeg("CS_comp_BIC.jpeg", quality = 100, units = "in", width = 10, height = 8, res = 300)

plot(cs_fit, xlab="Cox-Snell Residual", ylab="Survival",
     main="Competing Cox model",xlim = c(0,3.2),col="blue", lwd=2)

# Theoretical Exp(1): S(r) = exp(-r)
curve(exp(-x), col="black", lwd=3, add=TRUE)

legend("topright", legend=c("KM of residuals", "Exp(1)"),
       col=c("blue","black"), lwd=2)

dev.off()


