library(survival)
library(casebase)
data(support)
attach(support)
source("Cox_competing_correct_cluster.R")
source("Cox_competing_algorithm.R")
library(riskRegression)
library(MASS)
data <- support[support$dzgroup == "Coma", ]
data$sex <- ifelse(data$sex=="male",1,0)
data$ID <- seq(1,nrow(data),1)
data$status <- data$death
data$time <- data$d.time
kfold_cox(data, formula,id_var = "ID", k = 5, seed = 1234) 


kfold_cox <- function(data, formula, id_var = "ID", k = 5, seed = 123) {
  set.seed(seed)
  n <- nrow(data)
  
  fold_ids <- sample(rep(1:k, length.out = n))
  
  predictions <- numeric(n)
  
  fold_iAUC  <- numeric(k)
  fold_brier <- numeric(k)
  
  for (fold in 1:k) {
    
    cat(sprintf("Processing fold %d/%d...\n", fold, k))
    
    test_idx <- which(fold_ids == fold)
    
    train <- data[-test_idx, , drop = FALSE]
    test  <- data[test_idx, , drop = FALSE]
    
    formula_local <- as.formula(
      deparse(formula),
      env = environment()
    )
    
    fit <- suppressWarnings(
      tryCatch(
        survival::coxph(
          formula = formula_local,
          data = train
        ),
        error = function(e) {
          cat(sprintf("Model fit failed in fold %d: %s\n",
                      fold, e$message))
          return(NULL)
        }
      )
    )
    

    baseline_haz <- basehaz(fit, centered = T)  
    
    if (is.null(fit)) {
      fold_iAUC[fold]  <- NA
      fold_brier[fold] <- NA
      next
    }
    
    
    ############################################################
    # Test prediction
    ############################################################
    
    test_preds <- predict(fit, newdata=test, type = "lp")
    
    predictions[test_idx] <- test_preds
    
    ############################################################
    # Evaluation times
    ############################################################
    
    if (sum(test$status) == 0) {
      
      fold_iAUC[fold]  <- NA
      fold_brier[fold] <- NA
      
      next
    }
    
    time_points <- unique(
      quantile(
        test$time,
        probs = seq(0.01,0.99,by=0.1)
      )
    )
    
    ############################################################
    # Construct cumulative hazard matrix
    ############################################################
    
    hazard_matrix <- matrix(
      0,
      nrow = nrow(test),
      ncol = length(time_points)
    )
    
    for (m in seq_along(time_points)) {
      
      tp <- time_points[m]
      
      idx <- max(
        c(0, which(baseline_haz$time <= tp))
      )
      
      H0_tp <- if (idx == 0) {
        0
      } else {
        baseline_haz$hazard[idx]
      }
      
      hazard_matrix[, m] <-
        H0_tp * exp(test_preds)
    }
    
    colnames(hazard_matrix) <-
      paste0("t", seq_along(time_points))
    
    ############################################################
    # AUC
    ############################################################
    
    auc_res <- Score(
      object = list(
        Cox = hazard_matrix
      ),
      formula = Surv(time, status) ~ 1,
      data = test,
      times = time_points,
      metrics = "AUC",
      cens.model = "km"
    )
    
    fold_iAUC[fold] <- mean(
      auc_res$AUC$score$AUC,
      na.rm = TRUE
    )
    
    ############################################################
    # Brier Score
    ############################################################
    
    risk_matrix <- 1 - exp(-hazard_matrix)
    
    bs_res <- Score(
      object = list(
        Cox = risk_matrix
      ),
      formula = Surv(time, status) ~ 1,
      data = test,
      times = time_points,
      metrics = "Brier",
      cens.model = "km"
    )
    
    fold_brier[fold] <- mean(
      bs_res$Brier$score$Brier[
        bs_res$Brier$score$model == "Cox"
      ],
      na.rm = TRUE
    )
    
    cat(
      sprintf(
        " Fold %d: iAUC = %.4f, BS = %.4f\n",
        fold,
        fold_iAUC[fold],
        fold_brier[fold]
      )
    )
  }
  
  list(
    fold_assignments = fold_ids,
    fold_iAUC        = fold_iAUC,
    fold_brier       = fold_brier,
    mean_fold_iAUC   = mean(fold_iAUC, na.rm = TRUE),
    mean_fold_brier  = mean(fold_brier, na.rm = TRUE)
  )
}





repeat_kfold_cox <- function(data, formula, id_var = "ID", k = 5, seeds = NULL) {
  
  if (is.null(seeds)) {
    seeds <- sample.int(1e6, 100)
  }
  
  n_rep <- length(seeds)
  
  res_mat <- matrix(
    NA,
    nrow = n_rep,
    ncol = 2
  )
  
  colnames(res_mat) <- c("iAUC", "Brier")
  
  ## store fold assignments
  all_fold_ids <- vector("list", n_rep)
  
  for (i in seq_along(seeds)) {
    
    cat(sprintf(
      "\nRunning CV repetition %d/%d (seed = %d)\n",
      i, n_rep, seeds[i]
    ))
    
    fit_i <- kfold_cox(
      data    = data,
      formula = as.formula(formula),
      id_var  = id_var,
      k       = k,
      seed    = seeds[i]
    )
    
    ## save fold ids
    all_fold_ids[[i]] <- fit_i$fold_assignments
    
    res_mat[i, ] <- c(
      fit_i$mean_fold_iAUC,
      fit_i$mean_fold_brier
    )
  }
  
  res_df <- as.data.frame(res_mat)
  
  summary_df <- data.frame(
    metric = c("iAUC", "Brier"),
    mean = c(
      mean(res_df$iAUC, na.rm = TRUE),
      mean(res_df$Brier, na.rm = TRUE)
    ),
    sd = c(
      sd(res_df$iAUC, na.rm = TRUE),
      sd(res_df$Brier, na.rm = TRUE)
    )
  )
  
  summary_df$se <-
    summary_df$sd /
    sqrt(colSums(!is.na(res_df)))
  
  list(
    raw_results      = res_df,
    summary          = summary_df,
    fold_assignments = all_fold_ids,
    seeds            = seeds
  )
}


set.seed(1234)
seeds_100 <- sample(1:10000, 10)

formula = as.formula(Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                       resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                       bun + urine + adlp + adlsc+ sps + aps)

cv100_res_full <- repeat_kfold_cox(
  data    = data,
  formula = formula,
  id_var  = "ID",
  #k       = 5,
  seeds   = seeds_100
)


cv100_res_full$summary

formula = as.formula(Surv(d.time, death) ~ age + scoma + hrt + temp + adlsc + sps + 
                       aps)

cv100_res_AIC <- repeat_kfold_cox(
  data    = data,
  formula = formula,
  id_var  = "ID",
  k       = 5,
  seeds   = seeds_100
)

cv100_res_AIC$summary


#FULL
#1   iAUC 0.7959759 0.004440129
#2  Brier 0.1478048 0.001684098




kfold_weibull <- function(data, formula, id_var = "ID", k = 5, seed = 123) {
  
  set.seed(seed)
  n <- nrow(data)
  
  fold_ids <- sample(rep(1:k, length.out = n))
  
  fold_iAUC  <- numeric(k)
  fold_brier <- numeric(k)
  
  
  for (fold in 1:k) {
    
    cat(sprintf("Processing fold %d/%d...\n", fold, k))
    
    test_idx <- which(fold_ids == fold)
    
    train <- data[-test_idx, , drop = FALSE]
    test  <- data[test_idx, , drop = FALSE]
    
    
    ############################################################
    # Fit Weibull AFT model
    ############################################################
    
    fit <- suppressWarnings(
      tryCatch(
        survival::survreg(
          formula,
          data = train,
          dist = "weibull"
        ),
        error = function(e) {
          cat(sprintf("Model fit failed in fold %d\n", fold))
          return(NULL)
        }
      )
    )
    
    
    if (is.null(fit)) {
      fold_iAUC[fold]  <- NA
      fold_brier[fold] <- NA
      next
    }
    
    
    ############################################################
    # Prediction
    ############################################################
    
    lp <- predict(
      fit,
      newdata = test,
      type = "lp"
    )
    
    sigma <- fit$scale
    
    
    ############################################################
    # Evaluation times
    ############################################################
    
    if (sum(test$status) == 0) {
      
      fold_iAUC[fold]  <- NA
      fold_brier[fold] <- NA
      
      next
    }
    
    
    time_points <- unique(
      quantile(
        test$time,
        probs = seq(0.01, 0.99, by = 0.1)
      )
    )
    
    
    ############################################################
    # Construct cumulative hazard matrix
    #
    # Weibull:
    # H(t|X) = (t / exp(lp))^(1/sigma)
    ############################################################
    
    hazard_matrix <- matrix(
      0,
      nrow = nrow(test),
      ncol = length(time_points)
    )
    
    
    for (m in seq_along(time_points)) {
      
      tp <- time_points[m]
      
      hazard_matrix[, m] <-
        (tp / exp(lp))^(1 / sigma)
      
    }
    
    
    colnames(hazard_matrix) <-
      paste0("t", seq_along(time_points))
    
    
    ############################################################
    # Survival probability matrix
    ############################################################
    
    survival_matrix <- 1-exp(-hazard_matrix)
    
    
    ############################################################
    # AUC
    ############################################################
    
    auc_res <- Score(
      object = list(
        Weibull = hazard_matrix
      ),
      formula = Surv(time, status) ~ 1,
      data = test,
      times = time_points,
      metrics = "AUC",
      cens.model = "km"
    )
    
    
    fold_iAUC[fold] <- mean(
      auc_res$AUC$score$AUC,
      na.rm = TRUE
    )
    
    
    ############################################################
    # Brier score
    ############################################################
    
    bs_res <- Score(
      object = list(
        Weibull = survival_matrix
      ),
      formula = Surv(time, status) ~ 1,
      data = test,
      times = time_points,
      metrics = "Brier",
      cens.model = "km"
    )
    
    
    fold_brier[fold] <- mean(
      bs_res$Brier$score$Brier[
        bs_res$Brier$score$model == "Weibull"
      ],
      na.rm = TRUE
    )
    
    
    cat(
      sprintf(
        " Fold %d: iAUC = %.4f, BS = %.4f\n",
        fold,
        fold_iAUC[fold],
        fold_brier[fold]
      )
    )
    
  }
  
  
  list(
    fold_assignments = fold_ids,
    fold_iAUC        = fold_iAUC,
    fold_brier       = fold_brier,
    mean_fold_iAUC   = mean(fold_iAUC, na.rm = TRUE),
    mean_fold_brier  = mean(fold_brier, na.rm = TRUE)
  )
}



repeat_kfold_Weibull <- function(data, formula, k = 5, seeds = NULL) {
  
  if (is.null(seeds)) {
    seeds <- sample.int(1e6, 100)
  }
  
  n_rep <- length(seeds)
  
  res_mat <- matrix(
    NA,
    nrow = n_rep,
    ncol = 2
  )
  
  colnames(res_mat) <- c("iAUC", "Brier")
  
  ## store fold assignments
  all_fold_ids <- vector("list", n_rep)
  
  for (i in seq_along(seeds)) {
    
    cat(sprintf(
      "\nRunning CV repetition %d/%d (seed = %d)\n",
      i, n_rep, seeds[i]
    ))
    
    fit_i <- kfold_weibull(
      data    = data,
      formula = formula,
      k       = k,
      seed    = seeds[i]
    )
    
    ## save fold ids
    all_fold_ids[[i]] <- fit_i$fold_assignments
    
    res_mat[i, ] <- c(
      fit_i$mean_fold_iAUC,
      fit_i$mean_fold_brier
    )
  }
  
  res_df <- as.data.frame(res_mat)
  
  summary_df <- data.frame(
    metric = c("iAUC", "Brier"),
    mean = c(
      mean(res_df$iAUC, na.rm = TRUE),
      mean(res_df$Brier, na.rm = TRUE)
    ),
    sd = c(
      sd(res_df$iAUC, na.rm = TRUE),
      sd(res_df$Brier, na.rm = TRUE)
    )
  )
  
  summary_df$se <-
    summary_df$sd /
    sqrt(colSums(!is.na(res_df)))
  
  list(
    raw_results      = res_df,
    summary          = summary_df,
    fold_assignments = all_fold_ids,
    seeds            = seeds
  )
}





set.seed(1234)
seeds_100 <- sample(1:10000, 100)

formula = as.formula(Surv(d.time,death)~ age+ sex+ scoma + meanbp + wblc + hrt+
                       resp + temp + pafi + alb + bili + crea + sod + ph + glucose+
                       bun + urine + adlp + adlsc+ sps + aps)

cv100_res_full_W <- repeat_kfold_Weibull(
  data    = data,
  formula = formula,
  #k       = 5,
  seeds   = seeds_100
)


cv100_res_full_W$summary

##1   iAUC 0.7851412 0.004523243 0.0004523243
#2  Brier 0.1665397 0.001647146 0.0001647146


formula = as.formula(Surv(d.time, death) ~ age + scoma + hrt + temp + adlsc + sps + 
                       aps)

cv100_res_W_AIC <- repeat_kfold_Weibull(
  data    = data,
  formula = formula,
  seeds   = seeds_100
)

cv100_res_W_AIC$summary
##  metric      mean          sd           se
#1   iAUC 0.7960926 0.003870993 0.0003870993
#2  Brier 0.1627102 0.001455274 0.0001455274












kfold_competing <- function(data, covariate_names, label, fold_ids){
  
  k <- length(unique(fold_ids))
  
  fold_iAUC <- numeric(k)
  fold_iBS  <- numeric(k)
  
  for(fold in sort(unique(fold_ids))){
    
    cat(sprintf(
      "Processing fold %d/%d\n",
      fold,
      k
    ))
    
    test_idx <- which(fold_ids == fold)
    
    train <- data[-test_idx,,drop=FALSE]
    test  <- data[test_idx,,drop=FALSE]
    
    ################################################
    # Transform
    ################################################
    
    transformed_train <- data_transform(
      train$time,
      train$status,
      as.matrix(train[,covariate_names]),
      label,      cluster=as.numeric(train$ID)

    )
    
    transformed_test <- data_transform(
      test$time,
      test$status,
      as.matrix(test[,covariate_names]),
      label,
      cluster=as.numeric(test$ID)
    )
   
    ################################################
    # Fit
    ################################################
    result <- Cox_com_fit(transformed_train, maxit = 1000,tolerance = 1e-4, frailty = T)
    result$cluster <- transformed_train$cluster
    ################################################
    # Evaluation times
    ################################################
    
    time_points <- sort(unique(as.numeric(quantile(test$time,probs = seq(0.01,0.99, by = 0.1)))))
    
    ################################################
    # Predict on TEST
    ################################################
    
    pred <- calculate_total_hazard(time_points,result,transformed_test )
    
    ################################################
    # AUC
    ################################################
    
    auc_res <- riskRegression::Score(
      object = list(
        Competing = pred$hazard_matrix
      ),
      formula = Surv(time,status)~1,
      data = data.frame(
        time = transformed_test$time,
        status = transformed_test$status
      ),
      times = time_points,
      metrics = "AUC",
      cens.model = "km"
    )
    
    fold_iAUC[fold] <- mean(
      auc_res$AUC$score$AUC,
      na.rm = TRUE
    )
    
    ################################################
    # Brier
    ################################################
    
    bs_res <- riskRegression::Score(
      object = list(
        Competing = 1 - pred$S_matrix
      ),
      formula = Surv(time,status)~1,
      data = data.frame(
        time = transformed_test$time,
        status = transformed_test$status
      ),
      times = time_points,
      metrics = "Brier",
      cens.model = "km"
    )
    
    fold_iBS[fold] <- mean(
      bs_res$Brier$score$Brier[
        bs_res$Brier$score$model=="Competing"
      ],
      na.rm = TRUE
    )
    
    cat(sprintf(
      "iAUC = %.4f, iBS = %.4f\n",
      fold_iAUC[fold],
      fold_iBS[fold]
    ))
  }
  
  list(
    fold_iAUC = fold_iAUC,
    fold_iBS = fold_iBS,
    mean_iAUC = mean(
      fold_iAUC,
      na.rm=TRUE
    ),
    mean_iBS = mean(
      fold_iBS,
      na.rm=TRUE
    )
  )
}

#1 3 4 8 13 15
#2 5 8 9 10 11 14
#3 4 6 7 9 10 12

#scoma age sps adlsc sod alb
#aps hrt adlsc bili urine bun meanbp
#age sps temp pafi bili urine ph


covariate_names = c("scoma", "aps","age", "sps",  "hrt",
                    "temp", "pafi", "adlsc", "bili", "urine", 
                    "bun", "ph","sod","meanbp", "alb"
                    )

load("label_can.RData")

covariate <- covariate_names




kfold_competing(data, covariate, label=label_can, fold_ids=cv100_res_full$fold_assignments[[1]])


n_rep <- length(cv100_res_AIC$fold_assignments)

res_comp <- matrix(NA, nrow = n_rep, ncol = 2)
colnames(res_comp) <- c("iAUC", "iBS")

for (i in seq_len(n_rep)) {
  
  cat(sprintf("\nRunning competing CV repetition %d/%d\n", i, n_rep))
  
  fold_ids <- cv100_res_AIC$fold_assignments[[i]]
  
  fit_i <- kfold_competing(
    data = data,
    covariate_names = covariate_names,
    label = label_can,
    fold_ids = fold_ids
  )
  
  res_comp[i, ] <- c(
    fit_i$mean_iAUC,
    fit_i$mean_iBS
  )
}

res_comp <- as.data.frame(res_comp)

summary_comp <- data.frame(
  metric = c("iAUC", "iBS"),
  mean = c(
    mean(res_comp$iAUC, na.rm = TRUE),
    mean(res_comp$iBS, na.rm = TRUE)
  ),
  sd = c(
    sd(res_comp$iAUC, na.rm = TRUE),
    sd(res_comp$iBS, na.rm = TRUE)
  )
)

summary_comp$se <- summary_comp$sd / sqrt(n_rep)

summary_comp
##BIC
#metric      mean          sd           se
#1   iAUC 0.8135458 0.004190018 0.0010995759
#2    iBS 0.1434609 0.001514079 0.0002732985





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

  
source("Functions_new.txt")

kfold_competing_W <- function(data, covariate_names, label, fold_ids){
  
  k <- length(unique(fold_ids))
  
  fold_iAUC <- numeric(k)
  fold_iBS  <- numeric(k)
  
  for(fold in sort(unique(fold_ids))){
    
    cat(sprintf(
      "Processing fold %d/%d\n",
      fold,
      k
    ))
    
    test_idx <- which(fold_ids == fold)
    
    train <- data[-test_idx,,drop=FALSE]
    test  <- data[test_idx,,drop=FALSE]
    
    ################################################
    # Transform
    ################################################
    
    transformed_train <- data_transform(
      log(train$time),
      train$status,
      as.matrix(train[,covariate_names]),
      label
    )
    
    transformed_test <- data_transform(
      log(test$time),
      test$status,
      as.matrix(test[,covariate_names]),
      label
    )
    
    ################################################
    # Fit
    ################################################
    result <- tryCatch(
      {
        EM_fit(
          transformed_train,
          init_param = param_init(
            train[, covariate_names],
            label
          ),
          lambda1 = 0,
          lambda2 = 0,
          maxit = 10^3,
          tolerance = 1e-4
        )
      },
      error = function(e) {
        NULL
      }
    )
    if (is.null(result)) {
      cat(sprintf(
        "Model fitting failed in fold %d\n",
        fold
      ))
      
      fold_iAUC[fold] <- NA
      fold_iBS[fold] <- NA
      
      next
    }
    ################################################
    # Evaluation times
    ################################################
    
    time_points <- sort(unique(as.numeric(quantile(test$time,probs = seq(0.01,0.99, by = 0.1)))))
    
    ################################################
    # Predict on TEST
    ################################################
    
    pred <- calculate_total_hazard_weibull(time_points,result,transformed_test )
    
    ################################################
    # AUC
    ################################################
    
    auc_res <- riskRegression::Score(
      object = list(
        Competing = pred$hazard_matrix
      ),
      formula = Surv(time,status)~1,
      data = data.frame(
        time = exp(transformed_test$Y),
        status = transformed_test$delta
      ),
      times = time_points,
      metrics = "AUC",
      cens.model = "km"
    )
    
    fold_iAUC[fold] <- mean(
      auc_res$AUC$score$AUC,
      na.rm = TRUE
    )
    
    ################################################
    # Brier
    ################################################
    
    bs_res <- riskRegression::Score(
      object = list(
        Competing = 1 - pred$S_matrix
      ),
      formula = Surv(time,status)~1,
      data = data.frame(
        time = exp(transformed_test$Y),
        status = transformed_test$delta
      ),
      times = time_points,
      metrics = "Brier",
      cens.model = "km"
    )
    
    fold_iBS[fold] <- mean(
      bs_res$Brier$score$Brier[
        bs_res$Brier$score$model=="Competing"
      ],
      na.rm = TRUE
    )
    
    cat(sprintf(
      "iAUC = %.4f, iBS = %.4f\n",
      fold_iAUC[fold],
      fold_iBS[fold]
    ))
  }
  
  list(
    fold_iAUC = fold_iAUC,
    fold_iBS = fold_iBS,
    mean_iAUC = mean(
      fold_iAUC,
      na.rm=TRUE
    ),
    mean_iBS = mean(
      fold_iBS,
      na.rm=TRUE
    )
  )
}


kfold_competing_W(data, covariate, label=label_can, fold_ids=cv100_res_full_W$fold_assignments[[12]])


n_rep <- length(cv100_res_full_W$fold_assignments)

res_comp <- matrix(NA, nrow = n_rep, ncol = 2)
colnames(res_comp) <- c("iAUC", "iBS")

for (i in seq_len(n_rep)) {
  
  cat(sprintf("\nRunning competing CV repetition %d/%d\n", i, n_rep))
  
  fold_ids <- cv100_res_full_W$fold_assignments[[i]]
  
  fit_i <- kfold_competing_W(
    data = data,
    covariate_names = covariate_names,
    label = label_can,
    fold_ids = fold_ids
  )
  
  res_comp[i, ] <- c(
    fit_i$mean_iAUC,
    fit_i$mean_iBS
  )
}

res_comp <- as.data.frame(res_comp)

summary_comp <- data.frame(
  metric = c("iAUC", "iBS"),
  mean = c(
    mean(res_comp$iAUC, na.rm = TRUE),
    mean(res_comp$iBS, na.rm = TRUE)
  ),
  sd = c(
    sd(res_comp$iAUC, na.rm = TRUE),
    sd(res_comp$iBS, na.rm = TRUE)
  )
)

summary_comp$se <- summary_comp$sd / sqrt(n_rep)

summary_comp

#  metric      mean         sd          se
#1   iAUC 0.7476223 0.02542069 0.002542069
#2    iBS 0.2607665 0.06403423 0.006403423


save(cv100_res_full_W, file = "cv100_res_full_W.RData")
save(cv100_res_full, file = "cv100_res_full.RData")
#save(res_comp, file = "res_comp.RData")
save(res_comp, file = "res_comp_W.RData")
