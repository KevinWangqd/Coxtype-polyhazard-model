data_transform <- function(time, status, X, label, param=NULL) {
  # Check for NA values in inputs
  if (any(is.na(time))) stop("Time contains NA values.")
  if (any(is.na(status))) stop("status contains NA values.")
  if (any(is.na(X))) stop("X contains NA values.")
  if (any(sapply(label, is.null))) stop("label contains NULL values.")
  
  # Check the same length
  if (length(time) != length(status)) {
    stop("Time and status must have the same length.")
  }
  
  if (nrow(X) != length(time)) {
    stop("The number of samples in X must match the length of time.")
  }
  
  if (!is.list(label) || length(label) != ncol(X)) {
    stop("label must be a list of the same length as the number of columns in X.")
  }
  
  ## --- Order by time ---
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]
  X <- X[ord, , drop = FALSE]
 
  ## --- Group structure ---
  groups <- unique(unlist(label))
  L <- length(groups)
  
  beta <- vector("list", L)
  lambda <- vector("list", L)
  cum_lambda <- vector("list", L)
  bh_times <- vector("list", L)
  X_list <- vector("list", L)
  
  
  ## --- Construct group-specific objects ---
  for (l in seq_along(groups)) {
    # Subset the corresponding features for the current group
    X_list[[l]] <- X[, which(sapply(label, function(x) groups[l] %in% x)), drop = FALSE]
    if (is.null(param)) {
      # Initialize  beta, lambda
      beta[[l]] <- rep(0, ncol(X_list[[l]]))   
      bh_times[[l]] <- sort(unique(time))
      lambda[[l]] <- rep(1/length(bh_times[[l]]), length(bh_times[[l]]))
    }
    else {
      # Use provided param values
      beta[[l]] <- param$beta[[l]]
      bh_times[[l]] <- param$bh_times[[l]]
      lambda[[l]] <- param$lambda[[l]]
    } 
    cum_lambda[[l]] <- cumsum(lambda[[l]])
  }
  

  init_param <- list(
    beta = beta,
    bh_times = bh_times,
    cum_lambda = cum_lambda,
    lambda = lambda)
  
  # Return a list of transformed data and initialized parameters
  return(list(
    time = time,            
    status = status, 
    X = X_list, init_param = init_param
  ))
}



Cox_com_fit <- function(transformed_data, penalty=0, maxit = 1000, tolerance = 1e-6, frailty = FALSE, sigma2=1) {
  
  # Number of groups
  L <- length(transformed_data$X)
  
  # Initialize parameter lists
  param_new <- vector("list", L)
  
  for (l in 1:L) {
    param_new[[l]] <- list(
      bh_times = transformed_data$init_param$bh_times[[l]],
      lambda = transformed_data$init_param$lambda[[l]],
      cum_lambda = transformed_data$init_param$cum_lambda[[l]],
      beta   = transformed_data$init_param$beta[[l]]
    )
  }
  p_total <- sum(sapply(transformed_data$X, ncol))
  #eps = max(min(0.05, 0.25 /p_total),0.01)
  eps = ifelse(p_total<8, 0.05, ifelse(p_total<10, 0.025, 0.01))
  
  # Convergence tracker for each group
  stop_vec <- rep(FALSE, L)
  
  
  # Initialize frailty variance (global)
  if (frailty) {
    sigma2 = sigma2
  } else {
    sigma2 = NULL
  }

  # EM algorithm
  for (iter in 1:maxit) {
    
    # E-step: compute eta matrix
    if (frailty) {
      eta_matrix <- eta_matrix_fun(
        time   = transformed_data$time,
        status = transformed_data$status,
        X      = transformed_data$X,
        bh_times = lapply(param_new, function(x) x$bh_times),
        lambda = lapply(param_new, function(x) x$lambda),
        beta   = lapply(param_new, function(x) x$beta),
        eps = eps
      )
      
      b_hat_fit <- update_frailty_gaussian_sigma(
        time   = transformed_data$time,
        status = transformed_data$status,
        X      = transformed_data$X,
        bh_times = lapply(param_new, function(x) x$bh_times),
        cum_lambda   = lapply(param_new, function(x) x$cum_lambda),
        beta     = lapply(param_new, function(x) x$beta),
        sigma2    = sigma2)  
      
      b_hat <- b_hat_fit$b
      sigma2 <- b_hat_fit$sigma2

    } else {
      eta_matrix <- eta_matrix_fun(
        time   = transformed_data$time,
        status = transformed_data$status,
        X      = transformed_data$X,
        bh_times = lapply(param_new, function(x) x$bh_times),
        lambda = lapply(param_new, function(x) x$lambda),
        beta   = lapply(param_new, function(x) x$beta), eps = eps
      )
      b_hat <- NULL
    }

        # M-step: update each group's parameters
    for (l in 1:L) {
      
      param_updater <- update_l(
        time   = transformed_data$time,
        status = transformed_data$status,
        X_l    = transformed_data$X[[l]],
        eta_l  = eta_matrix[, l],
        b_hat  = b_hat
      )
      # Check that bh_times and lambda have the same length
      if (length(param_updater$bh_times) != length(param_updater$lambda_l)) {
        stop(paste0("Length mismatch in group ", l))
      }
      # Check convergence for beta
      beta_diff <- sum((param_updater$beta_l - param_new[[l]]$beta)^2)
      if (beta_diff < tolerance) stop_vec[l] <- TRUE
      
      # Update parameters
      param_new[[l]] <- list(
        bh_times = param_updater$bh_times, 
        lambda = param_updater$lambda_l,
        beta   = param_updater$beta_l,
        cum_lambda = param_updater$cum_lambda_l
      )
    }

    
    
    # If all groups have converged, stop iteration
    if (all(stop_vec)) {
      message(paste("Converged at iteration", iter))
      break
    }
  }
  
  # Return final parameters and last eta matrix
  return(list(
    param_new  = param_new,
    eta_matrix = eta_matrix,
    b_hat = b_hat
  ))
}



eta_matrix_fun <- function(time, status, X, bh_times, lambda, beta, eps = 5e-2) {
  n <- length(time)      # number of observations
  L <- length(X)      # number of groups
  
  
  eta_matrix <- matrix(NA, nrow = n, ncol = L)
  
  for (l in seq_len(L)) {
    #idx <- match(time, bh_times[[l]])
    idx <- findInterval(time, bh_times[[l]]) 
    lambda_aligned <- lambda[[l]][idx]
    
    # --- Compute η for group l ---
    eta_matrix[, l] <- lambda_aligned * exp(X[[l]] %*% beta[[l]])
  }
  eta_matrix[status == 0, ] <- 1 / L
  
  sum_eta <- rowSums(eta_matrix)
  eta_matrix <- eta_matrix / sum_eta
  
  eta_matrix <- pmax(pmin(eta_matrix, 1 - eps), eps)
  #eta_matrix <- pmax(eta_matrix, eps)
  
  return(eta_matrix)
}




update_l <- function(time, status, X_l, eta_l, b_hat){
  
  pseudo_data <- make_pseudo(time, status, X_l, eta_l, b_hat)
  
  covariate_names <-  names(pseudo_data)[!(colnames(pseudo_data) %in% c("time", "status", "eta", "b_hat"))]
  
  ## --------------------------------
  ## Offset term
  ## --------------------------------
  if (is.null(b_hat)) {
    
    # -------- No frailty --------
    cox_formula <- as.formula(
      paste(
        "survival::Surv(time, status) ~",
        paste(covariate_names, collapse = " + ")
      )
    )
    
  } else {
    # -------- With frailty (offset) --------
    cox_formula <- as.formula(
      paste(
        "survival::Surv(time, status) ~",
        paste(covariate_names, collapse = " + "),
        "+ offset(b_hat)"
      )
    )
  }
  
  fit <- survival::coxph(
    cox_formula,
    data    = pseudo_data,
    weights = pseudo_data$eta,
    ties    = "breslow",
    control = survival::coxph.control(
      eps = 1e-11,
      iter.max = 1000,
      toler.chol = 1e-14
    )
  )

  bh <- survival::basehaz(fit, centered = FALSE)  
  lambda_l <- c(bh$hazard[1], diff(bh$hazard))
  cum_lambda_l <- bh$hazard
    
  return(list(bh_times = bh$time, lambda_l = lambda_l, cum_lambda_l = cum_lambda_l, beta_l = fit$coefficients))
}




make_pseudo <- function(time, status, X, eta, b_hat) {
  # If X is a matrix, ensure column names exist
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }
  
  # Combine into data frame
  
  df <- data.frame(time = time, status = status, X, eta = eta)
  
  if (!is.null(b_hat)) {
    df$b_hat <- b_hat
  }
  
  # For censored subjects, set weight to 1
  df$eta[df$status == 0] <- 1
  
  # Create pseudo-censored copies for observed failures
  df_pseudo <- df[df$status == 1, ]          # select failures
  df_pseudo$eta <- 1 - df_pseudo$eta         # complementary weight
  df_pseudo$status <- 0
  
  # mark as censored
  # Combine original and pseudo observations
  df_combined <- rbind(df, df_pseudo)
  
  return(df_combined)
}




update_frailty_gaussian_sigma <- function(time, status, X, bh_times, cum_lambda, beta, sigma2) {
  
  n <- length(time)
  L <- length(X)
  
  # Compute cumulative hazards
  cumhaz_total <- numeric(n)
  for (l in seq_len(L)) {
    idx <- findInterval(time, bh_times[[l]])
    cumhaz_total <- cumhaz_total + cum_lambda[[l]][idx] * exp(drop(X[[l]] %*% beta[[l]]))
  }
  
  # Posterior mode of b_i
  b <- sapply(seq_len(n), function(i) {
    f <- function(bi) status[i] - cumhaz_total[i] * exp(bi) - bi / sigma2
    tryCatch(
      uniroot(f, lower = -20, upper = 20)$root,
      error = function(e) {
        max(min(log(max(status[i], 1e-6) / max(cumhaz_total[i], 1e-6)), 10), -10)
      }
    )
  })
  
  var_b <- 1 / (cumhaz_total * exp(b) + 1 / sigma2)
  
  # Update sigma^2
  sigma2_new <- mean(b^2 + var_b)
  
  return(list(b = b, sigma2 = sigma2_new))
}


update_lambdat <- function(fit_result, transformed_data, tol = 1e-4, max_iter = 100) {
  L <- length(fit_result)
  n <- length(transformed_data$time)
  
  # Extract beta and initial lambda from fit_result
  beta_list <- lapply(fit_result, function(x) x$beta)
  lambda_list <- lapply(fit_result, function(x) x$lambda)
  bh_times_list <- lapply(fit_result, function(x) x$bh_times)
  
  # Compute risk for each cause (exp(X beta))
  risk_list <- lapply(seq_len(L), function(l) {
    exp(transformed_data$X[[l]] %*% beta_list[[l]])
  })
  
  # Compute at-risk process Ys for each cause
  Ys_list <- lapply(risk_list, function(risk_vec) rev(cumsum(rev(risk_vec))))
  
  # Self-consistency loop
  for (iter in seq_len(max_iter)) {
    
    # Step 1: compute eta matrix (numerator weights)
    eta_matrix <- eta_matrix_fun(
      time     = transformed_data$time,
      status   = transformed_data$status,
      X        = transformed_data$X,
      bh_times = bh_times_list,
      lambda   = lambda_list,
      beta     = beta_list
    )
    
    # Step 2: update lambda at bh_times using self-consistency
    lambda_new <- vector("list", L)
    for (l in seq_len(L)) {
      lambda_new[[l]] <- numeric(length(bh_times_list[[l]]))
      for (j in seq_along(bh_times_list[[l]])) {
        idx <- which(findInterval(transformed_data$time, bh_times_list[[l]]) == j)
        if (length(idx) > 0) {
          lambda_new[[l]][j] <- sum(eta_matrix[idx, l] / Ys_list[[l]][idx])
        }
      }
    }
    
    # Check convergence
    diff <- max(sapply(seq_len(L), function(l) max(abs(lambda_new[[l]] - lambda_list[[l]]))))
    lambda_list <- lambda_new
    if (diff < tol) break
  }
  
  return(lambda_list)
}


loglikelihood <- function(fit_result, transformed_data) {
  n <- length(transformed_data$time)
  L <- length(fit_result)
  
  # Extract beta, lambda, and bh_times from fit_result
  beta_list <- lapply(fit_result, function(x) x$beta)
  lambda_list <- lapply(fit_result, function(x) x$lambda)
  bh_times_list <- lapply(fit_result, function(x) x$bh_times)
  
  loglik <- 0
  
  # Compute linear predictors for each group
  linpred_list <- lapply(seq_len(L), function(l) {
    as.vector(transformed_data$X[[l]] %*% beta_list[[l]])
  })
  
  # Term 1: sum over failures
  for (i in seq_len(n)) {
    if (transformed_data$status[i] == 1) {
      tmp <- 0
      for (l in seq_len(L)) {
        idx <- findInterval(transformed_data$time[i], bh_times_list[[l]])
        idx <- max(idx, 1)
        tmp <- tmp + lambda_list[[l]][idx] * exp(linpred_list[[l]][i])
      }
      loglik <- loglik + log(tmp)
    }
  }
  
  # Term 2: cumulative hazard integral
  for (l in seq_len(L)) {
    bh_times_l <- bh_times_list[[l]]
    lambda_l <- lambda_list[[l]]
    X_l <- transformed_data$X[[l]]
    linpred_l <- linpred_list[[l]]
    
    for (j in seq_along(bh_times_l)) {
      t_j <- bh_times_l[j]
      at_risk <- which(transformed_data$time >= t_j)
      loglik <- loglik - lambda_l[j] * sum(exp(linpred_l[at_risk]))
    }
  }
  
  return(loglik)
}

profile_loglikelihood <- function(beta_list_fixed, fit_result, transformed_data, tol = 1e-4, max_iter = 100) {
  L <- length(fit_result)
  
  # Step 1: Update lambda to NPMLE given fixed beta
  # We replace the beta in fit_result with the fixed beta_list
  fit_beta_fixed <- fit_result
  for (l in seq_len(L)) {
    fit_beta_fixed[[l]]$beta <- beta_list_fixed[[l]]
  }
  
  # Step 2: Update lambda using self-consistency
  lambda_npmle <- update_lambdat(fit_beta_fixed, transformed_data, tol = tol, max_iter = max_iter)
  
  # Replace the lambda in fit_result with the updated lambda
  for (l in seq_len(L)) {
    fit_beta_fixed[[l]]$lambda <- lambda_npmle[[l]]
  }
  
  # Step 3: Compute log-likelihood with updated lambda
  loglik <- loglikelihood(fit_beta_fixed, transformed_data)
  
  return(loglik)
}


cov_beta_profile <- function(fit_result, transformed_data, eps = 1e-4, tol = 1e-2, fast=TRUE, max_iter = 100) {
  L <- length(fit_result)
  
  # Extract the estimated beta
  beta_hat <- lapply(fit_result, function(x) x$beta)
  
  # Flatten beta into a single vector for convenience
  beta_vec <- unlist(beta_hat)
  p <- length(beta_vec)
  
  # Initialize covariance matrix
  cov_mat <- matrix(0, nrow = p, ncol = p)
  
  # Helper to reconstruct beta_list from vector
  vector_to_beta_list <- function(vec) {
    beta_list <- vector("list", L)
    idx <- 1
    for (l in seq_len(L)) {
      pl <- length(fit_result[[l]]$beta)
      beta_list[[l]] <- vec[idx:(idx + pl - 1)]
      idx <- idx + pl
    }
    return(beta_list)
  }
  
  # Compute the profile log-likelihood at beta_hat
  loglik_0 <- profile_loglikelihood(beta_hat, fit_result, transformed_data, tol, max_iter)
  
  # Finite-difference to approximate Hessian
  for (i in seq_len(p)) {
    
    e_i <- rep(0, p)
    e_i[i] <- eps
    
    beta_plus  <- vector_to_beta_list(beta_vec + e_i)
    beta_minus <- vector_to_beta_list(beta_vec - e_i)
    
    loglik_plus  <- profile_loglikelihood(beta_plus, fit_result, transformed_data, tol, max_iter)
    loglik_minus <- profile_loglikelihood(beta_minus, fit_result, transformed_data, tol, max_iter)
    
    # Negative second derivative along i-th coordinate
    cov_mat[i, i] <- abs (loglik_plus - 2 * loglik_0 + loglik_minus) / (eps^2)
  }
  
  if(!fast){
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        e_i <- rep(0, p); e_i[i] <- eps
        e_j <- rep(0, p); e_j[j] <- eps
        
        beta_pp <- vector_to_beta_list(beta_vec + e_i + e_j)
        beta_pm <- vector_to_beta_list(beta_vec + e_i - e_j)
        beta_mp <- vector_to_beta_list(beta_vec - e_i + e_j)
        beta_mm <- vector_to_beta_list(beta_vec - e_i - e_j)
        
        log_pp <- profile_loglikelihood(beta_pp, fit_result, transformed_data, tol, max_iter)
        log_pm <- profile_loglikelihood(beta_pm, fit_result, transformed_data, tol, max_iter)
        log_mp <- profile_loglikelihood(beta_mp, fit_result, transformed_data, tol, max_iter)
        log_mm <- profile_loglikelihood(beta_mm, fit_result, transformed_data, tol, max_iter)
        
        cov_mat[i, j] <- cov_mat[j, i] <- - (log_pp - log_pm - log_mp + log_mm) / (4 * eps^2)
      }
    }
  }
  # Invert Hessian to get covariance matrix
  cov_mat <- solve(cov_mat)
  se_beta <- sqrt(diag(cov_mat))
  
  return(list(cov_mat= cov_mat, se_beta =  se_beta))
}



calculate_total_hazard <- function(time_points, result, transformed_data) {
  # time_points: vector of time points to evaluate
  L <- length(result$param_new)  # Number of causes
  n_samples <- nrow(transformed_data$X[[1]])
  
  # Extract beta_list once
  beta_list <- lapply(result$param_new, function(x) x$beta)
  
  # Pre-compute exp(Xβ) for each cause
  exp_Xbeta <- lapply(1:L, function(l) {
    as.numeric(exp(transformed_data$X[[l]] %*% beta_list[[l]]))
  })
  
  # For each time point
  hazard_matrix <- sapply(time_points, function(t) {
    total_hazard <- numeric(n_samples)
    
    for (l in 1:L) {
      # Find baseline hazard at time t for cause l
      # Use findInterval for efficiency
      idx <- findInterval(t, result$param_new[[l]]$bh_times)
      if (idx == 0) {
        baseline_hazard_l <- 0  # Before first baseline hazard time
      } else {
        baseline_hazard_l <- result$param_new[[l]]$cum_lambda[idx]
      }
      
      # Add cause-specific hazard
      total_hazard <- total_hazard + baseline_hazard_l * exp_Xbeta[[l]]
    }
    
    return(total_hazard)
  })
  
  # Add time points as column names
  colnames(hazard_matrix) <- paste0("t_", round(time_points, 2))
  
  S_matrix <- exp(-hazard_matrix)
  
  return(list(hazard_matrix= hazard_matrix, S_matrix =  S_matrix))
}

