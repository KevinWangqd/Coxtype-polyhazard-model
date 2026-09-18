source("Cox_competing_correct_cluster.R")

forward_select <- function(time, status, X_candidate, max_add = 10, lambda1 = 1, lambda2 = 1){
  
  # Step 1: Initialize with best single variable
  AIC_vec <- numeric(ncol(X_candidate))
  X_list <- vector("list", ncol(X_candidate))
  label_list <- vector("list", ncol(X_candidate))
  
  for (i in 1:ncol(X_candidate)) {
    X_list[[i]] <- as.matrix(X_candidate[, i])  
    label_list[[i]] <- list(1)
    transformed_data <- data_transform(time, status, X = X_list[[i]], label = list(1))

    model_res <- suppressMessages(
      Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-4)
    )
    AIC_vec[i] <- -2*loglikelihood(model_res$param_new, transformed_data)+ lambda1*1 + lambda2*1
  }
  
  AIC_old <- min(AIC_vec)
  X_old <- X_list[[which.min(AIC_vec)]]
  label_old <- label_list[[which.min(AIC_vec)]]
  
  
  # Step 2: Iteratively add variables
  for (step in 2:max_add) {
    message("Adding covariate in step ", step)
    
    AIC_vec <- numeric(ncol(X_candidate))
    X_list <- vector("list", ncol(X_candidate))
    label_list <- vector("list", ncol(X_candidate))
    
    for (i in 1:ncol(X_candidate)) {
      
      AIC_in_vec <- numeric(length(max(unlist(label_old))) + 1)
      X_in_list <- vector("list", length(AIC_in_vec))
      label_in_list <- vector("list", length(AIC_in_vec))
      
      for (j in 1:(length(AIC_in_vec) + 1)) {
        
        X_try <- Add_var(X = X_old, X_new = as.matrix(X_candidate[, i]), group_ind = j, label_old = label_old)
        
        if (X_try$added == FALSE) {
          AIC_in_vec[j] <- 1e8
          X_in_list[j] <- NULL
          label_in_list[j] <- NULL
          
        } else {
          transformed_data <- data_transform(time, status, X = X_try$X, label = X_try$label)
          model_res <- tryCatch({
            suppressMessages(
              Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-3)
            )
          }, error = function(e) {
            NULL
          })
          if (is.null(model_res)) {
            AIC_in_vec[j] <- 1e8
            X_in_list[[j]] <- NULL  
            label_in_list[[j]] <- NULL
          } else {
            AIC_in_vec[j] <- -2*loglikelihood(model_res$param_new, transformed_data) + lambda1*length(unlist(X_try$label))+ lambda2*max(unlist(X_try$label))
            X_in_list[[j]] <- X_try$X  
            label_in_list[[j]] <- X_try$label
          }
          
          
         
        }
      }
      
      AIC_vec[i] <- min(AIC_in_vec)
      X_list[[i]] <- X_in_list[[which.min(AIC_in_vec)]]
      label_list[[i]] <- label_in_list[[which.min(AIC_in_vec)]]
      
    }
    
    if (min(AIC_vec) < AIC_old - 1e-3) {
      AIC_old <- min(AIC_vec)
      X_old <- X_list[[which.min(AIC_vec)]]
      label_old <- label_list[[which.min(AIC_vec)]]
    } else {
      message("No improvement in forward selection at variable count = ", step)
      break  # Stop early if no improvement
    }
  }
  
  return(list(AIC = AIC_old, X = X_old, label = label_old))
}





Add_var <- function(X, X_new, group_ind, label_old) {
  X_copy <- as.matrix(X)
  success <- TRUE
  
  # Check if X_new was not added: No duplication will cause by addition
  dup_col <- which(sapply(1:ncol(X_copy), function(j) identical(as.vector(X_copy[, j]), as.vector(X_new))))
  
  if (length(dup_col) == 0) {
    # New covariate: add to X and create new label entry
    X_copy <- cbind(X_copy, X_new)
    label <- c(label_old, list(group_ind))
  } else {
    # Existing covariate: check within-group duplication
    if (group_ind %in% label_old[[dup_col]]) {
      success <- FALSE  # Duplicate: same covariate in same group; End
    } else {

      label_temp <- label_old
      label_temp[[dup_col]] <- c(label_temp[[dup_col]], group_ind)
      
      group_sets <- lapply(unique(unlist(label_temp)), function(g) {which(sapply(label_temp, function(x) g %in% x))})
      
      if (any(duplicated(lapply(group_sets, sort)))) {success <- FALSE}  ## across group duplication
      else { 
        label <- label_temp
      }}}
  
  return(list(X = if (success) X_copy else X, added = success, label = if (success) label else label_old))
}





high_select <- function(time, status, X, n_group = 3, size_group = 3, 
                        init_iauc = 0.5, in_iauc = 0.05, tar_P = 15, n_draw = 100){
  
  n_covariates <- ncol(X)
  m <- 0
  zeta <- init_iauc
  candidate_pool <- 1:n_covariates
  
  message("Starting high-dimensional selection with ", n_covariates, " covariates")
  message("Target size: ", tar_P, ", Initial threshold: ", zeta)
  repeat {
    message("\n--- Iteration ", m, " ---")
    message("Current candidate pool size: ", length(candidate_pool))
    message("Current candidate pool: ", paste(sort(candidate_pool), collapse = ", "))
    
    message("Current iAUC threshold: ", round(zeta, 4))

    iauc_vec <- numeric(n_draw)
    retained_covariates <- list()  # Store covariates from good draws
    
    for (r in 1:n_draw) {
      
      group_sizes <- sample(2:size_group, n_group, replace = TRUE)
      n_selected <- sum(group_sizes)
      
      if (n_selected > length(candidate_pool)) {
        selected_idx <- sample(candidate_pool, n_selected, replace = TRUE)
      } else {
        selected_idx <- sample(candidate_pool, n_selected, replace = FALSE) ## or FALSE to speed up
      }
      
      X_draw <- as.matrix(X[, selected_idx, drop = FALSE])
      
      label_draw <- vector("list", n_selected)
      cov_position <- 1
      for (l in 1:n_group) {
        for (k in 1:group_sizes[l]) {
          label_draw[[cov_position]] <- l
          cov_position <- cov_position + 1
        }
      }
      
      transformed_data <- data_transform(time, status, X = X_draw, label = label_draw)
      
      model_res <- tryCatch({
        suppressMessages(
          Cox_com_fit(transformed_data, maxit = 10^3, tolerance = 1e-3)
        )
      }, error = function(e) {
        NULL
      })
      
      # Compute iAUC: if model fit fails, set to 0.5 (below initial threshold)
      if (is.null(model_res)) {
        iauc_vec[r] <- 0.5
      } else {
        iauc_vec[r] <- compute_iAUC(transformed_data, result = model_res)
      }

      
      # If iAUC meets threshold, store the selected indices
      if (iauc_vec[r] >= zeta) {
        retained_covariates <- c(retained_covariates, list(selected_idx))
      }
    }
    
    if (length(retained_covariates) == 0) {
      warning("No draws meet the iAUC threshold at iteration ", m)
      message("Returning current selection with ", length(candidate_pool), " covariates")
      
      
      break
    }
    
    # Simple merge: combine all retained indices and remove duplicates
    candidate_pool <- unique(unlist(retained_covariates))
    
    m <- m + 1
    zeta <- init_iauc + m * in_iauc
    
    if (length(candidate_pool) <= tar_P || zeta >= 1) {
      break
      message("Finish selection with ", length(candidate_pool), " covariates")
    }
  }

  X_final <- as.matrix(X[, candidate_pool, drop = FALSE])
  
  return(list(
    selected_indices = candidate_pool,
    X_selected = X_final,
    n_selected = length(candidate_pool),
    iterations = m,
    final_threshold = zeta
  ))
}














