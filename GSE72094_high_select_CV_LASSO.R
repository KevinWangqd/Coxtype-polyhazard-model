source("Cox_competing_correct_cluster_LASSO.R")
source("Cox_competing_algorithm_CV.R")

load("LUAD_process.RData")
library(riskRegression)
library(MASS)
library(glmnet)
library(survival)

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




############################################################
# Recover covariate names after forward selection
############################################################

get_selected_names <- function(X_selected, X_candidate){
  
  selected_names <- character(ncol(X_selected))
  
  for(i in seq_len(ncol(X_selected))){
    
    target <- X_selected[,i]
    
    idx <- which(
      apply(
        X_candidate,
        2,
        function(x)
          all.equal(
            as.numeric(x),
            as.numeric(target)
          ) == TRUE
      )
    )
    
    if(length(idx)!=1){
      stop("Cannot uniquely identify selected covariate")
    }
    
    selected_names[i] <- colnames(X_candidate)[idx]
  }
  
  selected_names
}



############################################################
# Generate Cox prediction matrix
############################################################

generate_cox_prediction <- function(model, newdata, times){
  
  base <- basehaz(model, centered=FALSE)
  
  lp <- predict(
    model,
    newdata=newdata,
    type="lp"
  )
  
  H <- matrix(
    0,
    nrow=nrow(newdata),
    ncol=length(times)
  )
  
  
  for(j in seq_along(times)){
    
    idx <- max(
      which(base$time <= times[j]),
      0
    )
    
    H[,j] <- 
      if(idx==0) 0 else 
        base$hazard[idx] * exp(lp)
  }
  
  S <- exp(-H)
  
  list(
    hazard=H,
    survival=S
  )
}


############################################################
# One split train-test analysis
############################################################



n <- nrow(data)

train_idx <- sample(
  1:n,
  size = floor(0.8*n),
  replace = FALSE
)

test_idx <- setdiff(
  1:n,
  train_idx
)


train <- data[train_idx, ]
test  <- data[test_idx, ]

X_train <- X_filtered[train_idx,,drop=FALSE]
X_test  <- X_filtered[test_idx,,drop=FALSE]

x <- as.matrix(X_train)

#cv_lasso_cox <- cv.glmnet(x = as.matrix(X_train), y = survival::Surv(train$time_days, train$status), family = "cox",
#          alpha = 1, nfolds = 5)
#cv_lasso_cox$lambda.min 0.11

lasso_fit <- glmnet::glmnet(
  x = as.matrix(X_train),
  y = survival::Surv(train$time_days, train$status),
  family = "cox",
  alpha = 1, #LASSO
  lambda = 0.11
)

coef_lasso <- coef(lasso_fit)

informative_vars <- coef_lasso[coef_lasso[, 1] != 0, , drop = FALSE]

# View results
selected_vars <- rownames(coef(lasso_fit))[coef(lasso_fit)[, 1] != 0]

X_train_selected <- X_train[, selected_vars, drop = FALSE]

cox_classic <- coxph(Surv(train$time_days, train$status) ~ .,
  data = as.data.frame(X_train_selected))

# View results
summary(cox_classic)

cox_classic <- stepAIC(cox_classic,trace = FALSE, k = log(nrow(train)))

####################################################
# Prediction time points
####################################################

time_points <- quantile(test$time_days, probs = seq(0.1, 0.9, by = 0.05))


####################################################
# Classical Cox prediction
####################################################

pred_cox_train <- generate_cox_prediction(
  cox_classic,
  cbind(train,X_train[,selected_vars]),
  time_points
)

pred_cox <- generate_cox_prediction(
  cox_classic,
  cbind(test,X_test[, selected_vars]),
  time_points
)


test_eval <- data.frame(time=test$time_days,status=test$status)


####################################################
# iAUC Cox
####################################################


auc_cox <- Score(
  object=list(
    Cox=pred_cox$hazard
  ),
  formula=Surv(time,status)~1,
  data=data.frame(
    time=test$time_days,
    status=test$status
  ),
  times=time_points,
  metrics="AUC",
  cens.model="km"
)


mean(auc_cox$AUC$score$AUC,na.rm=TRUE)

auc_cox_train <- Score(object=list( Cox=pred_cox_train$hazard),
  formula=Surv(time,status)~1,
  data=data.frame(
    time=train$time_days,
    status=train$status
  ),
  times=time_points,
  metrics="AUC",
  cens.model="km"
)

 mean(auc_cox_train$AUC$score$AUC, na.rm=TRUE)


####################################################
# iBS Cox
####################################################


bs_cox <- Score(
  object=list(
    Cox=pred_cox$survival
  ),
  formula=Surv(time,status)~1,
  data=data.frame(
    time=test$time_days,
    status=test$status
  ),
  times=time_points,
  metrics="Brier",
  cens.model="km"
)


mean(bs_cox$Brier$score$Brier[bs_cox$Brier$score$model=="Cox"], na.rm=TRUE)

bs_cox_train <- Score(
  object=list(
    Cox=pred_cox_train$survival
  ),
  formula=Surv(time,status)~1,
  data=data.frame(
    time=train$time_days,
    status=train$status
  ),
  times=time_points,
  metrics="Brier",
  cens.model="km"
)


mean(bs_cox_train$Brier$score$Brier[bs_cox_train$Brier$score$model=="Cox"], na.rm=TRUE)






## 
forward_res <- forward_select(train$time_days, train$status,
  X_candidate = X_train[, selected_vars, drop = FALSE], max_add = 20,
  lambda1 = log(nrow(train)),
  lambda2 = 10*log(nrow(train)))

cov_names <- get_selected_names(forward_res$X, X_train[, selected_vars, drop = FALSE])


train_comp <- data_transform(train$time_days, train$status, as.matrix(X_train[,cov_names]),
  forward_res$label)


test_comp <- data_transform(test$time_days,test$status, as.matrix(X_test[,cov_names]), forward_res$label)

comp_fit <- Cox_com_fit(train_comp, maxit=1000, tolerance=1e-3, frailty=FALSE)


train_hazard <- calculate_total_hazard(time_points, comp_fit,train_comp)

test_hazard <- calculate_total_hazard(time_points, comp_fit, test_comp)


time_AUC_train <- riskRegression::Score(
  object=list(
    competing_Cox=train_hazard$hazard_matrix
  ),
  formula=Surv(time,status)~1,
  data=data.frame(
    time=train_comp$time,
    status=train_comp$status
  ),
  times=time_points,
  metrics="AUC",
  cens.model="km"
)

mean(time_AUC_train$AUC$score$AUC, na.rm=TRUE)

time_AUC <- riskRegression::Score(
  object   = list(
    competing_Cox = test_hazard$hazard_matrix 
  ),                 # n x K risk matrix
  formula  = Surv(time, status) ~ 1,
  data     = data.frame(
    time   = test_comp$time,    # T_i
    status = test_comp$status   # delta_i (1=event, 0=censored)
  ),
  times    = time_points,
  metrics  = "AUC",
  cens.model = "km"
)

mean(time_AUC$AUC$score$AUC,na.rm = T)

####################################################
# iBS competing
####################################################

bs_train <- riskRegression::Score(
  object=list(
    competingCox=1-train_hazard$S_matrix
  ),
  formula=Surv(time,status)~1,
  data=data.frame(
    time=train_comp$time,
    status=train_comp$status
  ),
  times=time_points,
  metrics="Brier",
  cens.model="km"
)


mean(bs_train$Brier$score$Brier[ bs_train$Brier$score$model=="competingCox"],na.rm=TRUE)

bs_res <- riskRegression::Score(
  object      = list(competingCox = 1 - test_hazard$S_matrix),
  formula     = Surv(time, status) ~ 1,
  data        = data.frame(
    time   = test_comp$time,
    status = test_comp$status
  ),
  times       = time_points,
  metrics     = "Brier",
  cens.model  = "km"
)
mean(bs_res$Brier$score$Brier[which(bs_res$Brier$score$model=="competingCox")],na.rm = T)




## competing LASSO
# number of column folds
K <- 200

p <- ncol(X_train)

# randomly split columns into 11 folds
set.seed(1235)
col_folds <- sample(rep(1:K, length.out = p))

# store results
lasso_results <- vector("list", K)

for(k in 1:10){
  
  cat("Running fold", k, "/", K, "\n")
  # selected columns in this fold
  idx_k <- which(col_folds == k)
  
  X_k <- X_train[, idx_k, drop = FALSE]
  
  
  comp_lasso <- tryCatch({
    
    lasso_competing_select(
      time = train$time_day,
      status = train$status,
      X = X_k,
      n_group = 3,
      penalty = 0.1/4,
      n_draw = 1
    )
    
  }, error = function(e){
    message("Fold ", k, " failed: ", e$message)
    NULL
  })
  
  
  if(!is.null(comp_lasso)){
    
    # map back to original feature names
    selected_local <- comp_lasso$selected_variables

    
    lasso_results[[k]] <- comp_lasso$selected_variables
    
  } else {
    
    lasso_results[[k]] <- character(0)
    
  }
}


# combine all selected variables
selected_all <- unique(unlist(lasso_results))



comp_lasso_final <- lasso_competing_select(
    time = train$time_day,
    status = train$status,
    X = X_train[,selected_all],
    n_group = 3,
    penalty = 0.1/4,
    n_draw = 20
  )


comp_lasso_final$selection_frequency
selected_vars <- comp_lasso_final$selected_variables[1:30]





##
run_cv_rep <- function(data, X_filtered, n_rep = 100){
  
  res <- vector("list", n_rep)
  
  n <- nrow(data)
  
  for(r in 1:n_rep){
    
    message("Replication ", r, "/", n_rep)
    
    tryCatch({
      
      ############################################
      # Train-test split
      ############################################
      
      train_idx <- sample(
        1:n,
        size = floor(0.8*n),
        replace = FALSE
      )
      
      test_idx <- setdiff(1:n, train_idx)
      
      train <- data[train_idx,]
      test  <- data[test_idx,]
      
      X_train <- X_filtered[train_idx,,drop=FALSE]
      X_test  <- X_filtered[test_idx,,drop=FALSE]
      
      
      ############################################
      # LASSO Cox screening
      ############################################
      
      #lasso_fit <- glmnet::glmnet(x = as.matrix(X_train),
      #  y = survival::Surv(train$time_days, train$status),
      #  family = "cox",
      #  alpha = 1,
      #  lambda = 0.11)
      
      
      #selected_vars <- rownames(coef(lasso_fit))[
      #  coef(lasso_fit)[,1] != 0]
      
      
      ############################################
      # Competing LASSO screening
      ############################################
      
      K <- 50
      p <- ncol(X_train)
      
      # random split of features into K blocks
      col_folds <- sample(rep(1:K, length.out = p))
      
      lasso_results <- vector("list", K)
      
      for(k in 1:K){
        
        message("Competing LASSO block ", k, "/", K)
        
        idx_k <- which(col_folds == k)
        
        X_k <- X_train[, idx_k, drop = FALSE]
        
        comp_lasso <- tryCatch({
        
          lasso_competing_select(
            time = train$time_day,
            status = train$status,
            X = X_k,
            n_group = 2,
            penalty = 0.11/4,
            n_draw = 1)
          
        }, error=function(e){
          
          message(
            "Block ", k,
            " failed: ",
            e$message
          )
          
          NULL
        })
        
        
        if(!is.null(comp_lasso)){
  
          #selected_local <- comp_lasso$selected_variables
          
          
          lasso_results[[k]] <- comp_lasso$selected_variables
          
          
        }else{
          
          lasso_results[[k]] <- character(0)
          
        }
        
      }
      
      
      ############################################
      # Merge block-level selections
      ############################################
      
      selected_all <- unique(
        unlist(lasso_results)
      )
      
      
      
      ############################################
      # Final competing LASSO selection
      ############################################
      
      comp_lasso_final <- lasso_competing_select(
        time = train$time_day,
        status = train$status,
        X = X_train[, selected_all, drop = FALSE],
        n_group = 2,
        penalty = 0.03,
        n_draw = 10
      )
      
      
      # selection frequency
      freq <- comp_lasso_final$selection_frequency
      
      
      # choose top variables
      selected_vars <- names(
        sort(
          freq,
          decreasing = TRUE
        )[1:min(20,length(freq))]
      )
      
      X_train_selected <- X_train[,selected_vars,drop=FALSE]
      
      ############################################
      # Classical Cox + BIC
      ############################################
      
      cox_classic <- survival::coxph(
        survival::Surv(time_days,status)~.,
        data = as.data.frame(
          cbind(
            train[,c("time_days","status")],
            X_train_selected
          )
        )
      )
      
      
      cox_classic <- MASS::stepAIC(
        cox_classic,
        trace = FALSE,
        k = log(nrow(train))
      )
      
      
      time_points <- quantile(
        test$time_days,
        probs=seq(0.1,0.9,by=0.05)
      )
      
      
      ############################################
      # Cox prediction
      ############################################
      
      pred_cox <- generate_cox_prediction(
        cox_classic,
        cbind(test,X_test[,selected_vars,drop=FALSE]),
        time_points
      )
      
      
      pred_cox_train <- generate_cox_prediction(
        cox_classic,
        cbind(train,X_train[,selected_vars,drop=FALSE]),
        time_points
      )
      
      
      ############################################
      # Cox performance
      ############################################
      
      auc_cox <- Score(
        object=list(Cox=pred_cox$hazard),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=test$time_days,
          status=test$status
        ),
        times=time_points,
        metrics="AUC",
        cens.model="km"
      )
      
      
      bs_cox <- Score(
        object=list(Cox=pred_cox$survival),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=test$time_days,
          status=test$status
        ),
        times=time_points,
        metrics="Brier",
        cens.model="km"
      )
      
      
      cox_iAUC <- mean(
        auc_cox$AUC$score$AUC,
        na.rm=TRUE
      )
      
      cox_iBS <- mean(
        bs_cox$Brier$score$Brier[
          bs_cox$Brier$score$model=="Cox"
        ],
        na.rm=TRUE
      )
      
      auc_cox_train <- Score(
        object=list(Cox=pred_cox_train$hazard),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=train$time_days,
          status=train$status
        ),
        times=time_points,
        metrics="AUC",
        cens.model="km"
      )
      
      cox_train_iAUC <- mean(
        auc_cox_train$AUC$score$AUC,
        na.rm=TRUE
      )
      
      
      bs_cox_train <- Score(
        object=list(Cox=pred_cox_train$survival),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=train$time_days,
          status=train$status
        ),
        times=time_points,
        metrics="Brier",
        cens.model="km"
      )
      
      cox_train_iBS <- mean(
        bs_cox_train$Brier$score$Brier[
          bs_cox_train$Brier$score$model=="Cox"
        ],
        na.rm=TRUE
      )
      ############################################
      # Forward competing Cox
      ############################################
      
      forward_res <- forward_select(
        train$time_days,
        train$status,
        X_candidate=X_train[,selected_vars,drop=FALSE],
        max_add=15,
        lambda1=log(nrow(train)),
        lambda2=10*log(nrow(train))
      )
      
      
      cov_names <- get_selected_names(
        forward_res$X,
        X_train[,selected_vars,drop=FALSE]
      )
      
      
      train_comp <- data_transform(
        train$time_days,
        train$status,
        as.matrix(X_train[,cov_names]),
        forward_res$label
      )
      
      test_comp <- data_transform(
        test$time_days,
        test$status,
        as.matrix(X_test[,cov_names]),
        forward_res$label
      )
      
      
      comp_fit <- Cox_com_fit(
        train_comp,
        maxit=1000,
        tolerance=1e-3,
        frailty=FALSE
      )
      
      
      train_hazard <- calculate_total_hazard(
        time_points,
        comp_fit,
        train_comp
      )
      
      test_hazard <- calculate_total_hazard(
        time_points,
        comp_fit,
        test_comp
      )
      
      
      ############################################
      # Competing Cox performance
      ############################################
      
      auc_comp <- Score(
        object=list(
          competingCox=test_hazard$hazard_matrix
        ),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=test_comp$time,
          status=test_comp$status
        ),
        times=time_points,
        metrics="AUC",
        cens.model="km"
      )
      
      
      bs_comp <- Score(
        object=list(
          competingCox=1-test_hazard$S_matrix
        ),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=test_comp$time,
          status=test_comp$status
        ),
        times=time_points,
        metrics="Brier",
        cens.model="km"
      )
      
      
      comp_iAUC <- mean(
        auc_comp$AUC$score$AUC,
        na.rm=TRUE
      )
      
      comp_iBS <- mean(
        bs_comp$Brier$score$Brier[
          bs_comp$Brier$score$model=="competingCox"
        ],
        na.rm=TRUE
      )
      
      auc_comp_train <- Score(
        object=list(
          competingCox=train_hazard$hazard_matrix
        ),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=train_comp$time,
          status=train_comp$status
        ),
        times=time_points,
        metrics="AUC",
        cens.model="km"
      )
      
      comp_train_iAUC <- mean(
        auc_comp_train$AUC$score$AUC,
        na.rm=TRUE
      )
      
      
      bs_comp_train <- Score(
        object=list(
          competingCox=1-train_hazard$S_matrix
        ),
        formula=Surv(time,status)~1,
        data=data.frame(
          time=train_comp$time,
          status=train_comp$status
        ),
        times=time_points,
        metrics="Brier",
        cens.model="km"
      )
      
      comp_train_iBS <- mean(
        bs_comp_train$Brier$score$Brier[
          bs_comp_train$Brier$score$model=="competingCox"
        ],
        na.rm=TRUE
      )
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      ############################################
      # Save results
      ############################################
      
      res[[r]] <- data.frame(
        rep=r,
        
        n_selected_lasso=length(selected_vars),
        n_selected_CF=length(cov_names),
        
        Cox_train_iAUC=cox_train_iAUC,
        Cox_train_iBS=cox_train_iBS,
        Cox_test_iAUC=cox_iAUC,
        Cox_test_iBS=cox_iBS,
        
        Comp_train_iAUC=comp_train_iAUC,
        Comp_train_iBS=comp_train_iBS,
        Comp_test_iAUC=comp_iAUC,
        Comp_test_iBS=comp_iBS
      )
      
      
    }, error=function(e){
      
      message(
        "Replication ",r,
        " failed: ",
        e$message
      )
      
      res[[r]] <- data.frame(
        rep=r,
        n_selected_lasso=NA,
        n_selected_CF=NA,
        Cox_iAUC=NA,
        Cox_iBS=NA,
        Comp_iAUC=NA,
        Comp_iBS=NA
      )
    })
    
  }
  
  
  return(
    do.call(
      rbind,
      res
    )
  )
}


set.seed(43139)

final_res3 <- run_cv_rep( 
  data = data,
  X_filtered = X_filtered,
  n_rep = 5 
)



final_res <- rbind(final_res,final_res3)
final_res <- final_res3
final_res1 <- rbind(final_res,final_res1)

summary(final_res$Cox_test_iAUC)
summary(final_res$Comp_test_iAUC)


boxplot(final_res$Cox_train_iAUC,final_res$Comp_train_iAUC,final_res$Cox_test_iAUC, final_res$Comp_test_iAUC)
boxplot(final_res$Cox_train_iBS,final_res$Comp_train_iBS,final_res$Cox_test_iBS, final_res$Comp_test_iBS)


#save(final_res, file = "LASSO_L5_50.RData")
