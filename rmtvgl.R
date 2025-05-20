#1. solve_tvgl.R
# TVGL baseline (Time-Varying Graphical Lasso) using ADMM

library(Matrix)
library(MASS)

soft_threshold <- function(X, tau) {
  return(sign(X) * pmax(abs(X) - tau, 0))
}

solve_tvgl <- function(S_list, lambda, beta, rho = 1, max_iter = 1000, tol = 1e-4, penalty_type = "l1") {
  T <- length(S_list)  # number of time points
  p <- nrow(S_list[[1]])  # number of variables
  
  # Initialize variables
  Theta_list <- lapply(1:T, function(t) diag(1, p))
  Z_list <- lapply(1:T, function(t) diag(1, p))
  U_list <- lapply(1:T, function(t) matrix(0, p, p))
  
  for (iter in 1:max_iter) {
    Theta_list_old <- Theta_list
    
    # Step 1: Theta update
    for (t in 1:T) {
      A <- S_list[[t]] - rho * (Z_list[[t]] - U_list[[t]])
      if (t > 1) A <- A + beta * Theta_list[[t-1]]
      if (t < T) A <- A + beta * Theta_list[[t+1]]
      
      eig <- eigen(A, symmetric = TRUE)
      d <- eig$values
      V <- eig$vectors
      d_new <- (d + sqrt(d^2 + 4*rho)) / (2*rho)
      Theta_list[[t]] <- V %*% diag(d_new) %*% t(V)
    }
    
    # Step 2: Z update (penalty-specific)
    for (t in 1:T) {
      V_mat <- Theta_list[[t]] + U_list[[t]]
      if (penalty_type == "l1") {
        Z_list[[t]] <- soft_threshold(V_mat, lambda / rho)
      } else if (penalty_type == "l2") {
        Z_list[[t]] <- V_mat / (1 + lambda / rho)
      } else {
        stop("Invalid penalty_type. Use 'l1' or 'l2'.")
      }
      diag(Z_list[[t]]) <- diag(V_mat)  # don't penalize diagonals
    }
    
    # Step 3: Dual update
    for (t in 1:T) {
      U_list[[t]] <- U_list[[t]] + Theta_list[[t]] - Z_list[[t]]
    }
    
    # Step 4: Check convergence
    primal_resid <- max(sapply(1:T, function(t) norm(Theta_list[[t]] - Z_list[[t]], "F")))
    dual_resid <- max(sapply(1:T, function(t) norm(Z_list[[t]] - Theta_list_old[[t]], "F")))
    
    if (primal_resid < tol && dual_resid < tol) {
      cat("TVGL converged at iteration", iter, "\n")
      break
    }
  }
  
  return(Theta_list)
}


#############
#2.# solve_rmtvgl.R

library(Matrix)

# Soft-thresholding operator
soft_threshold <- function(X, tau) {
  return(sign(X) * pmax(abs(X) - tau, 0))
}

# Huber loss function and gradient
huber_grad <- function(X, S, delta) {
  R <- X - S
  grad <- matrix(0, nrow = nrow(R), ncol = ncol(R))
  small_idx <- abs(R) <= delta
  large_idx <- !small_idx
  grad[small_idx] <- R[small_idx]
  grad[large_idx] <- delta * sign(R[large_idx])
  return(grad)
}

# Update Theta using Huber loss
update_theta_huber <- function(S, Z, U, Theta_prev, rho, beta, delta) {
  p <- nrow(S)
  max_iter <- 100
  tol <- 1e-4
  Theta <- Theta_prev
  
  for (iter in 1:max_iter) {
    G <- huber_grad(Theta, S, delta) + rho * (Theta - Z + U)
    if (!is.null(Theta_prev)) {
      G <- G + 2 * beta * (2*Theta - Theta_prev)
    }
    Theta_new <- Theta - 0.5 * G  # simple gradient step (step size=0.5)
    
    if (norm(Theta_new - Theta, type = "F") < tol) {
      break
    }
    Theta <- Theta_new
  }
  return((Theta + t(Theta)) / 2)  # ensure symmetry
}

# E-step: update empirical covariance estimates
em_estimate_S <- function(X_list) {
  S_list <- list()
  for (X in X_list) {
    # Column mean imputation
    for (j in 1:ncol(X)) {
      col_mean <- mean(X[, j], na.rm = TRUE)
      X[is.na(X[, j]), j] <- col_mean
    }
    
    # Safe covariance estimation
    S <- tryCatch(cov(X, use = "pairwise.complete.obs"),
                  error = function(e) diag(1e-4, ncol(X)))
    
    # Ensure square p x p matrix
    p <- ncol(X)
    if (!is.matrix(S) || any(dim(S) != c(p, p))) {
      S <- diag(1e-4, p)
    }
    
    # Symmetrize and ridge regularization
    S <- S + diag(1e-4, ncol(S))
    S <- (S + t(S)) / 2
    
    S_list[[length(S_list) + 1]] <- S
  }
  return(S_list)
}


# Main RMTVGL solver
solve_rmtvgl <- function(X_list, lambda, beta, delta = 1.0, rho = 1.0, 
                         em_tol = 1e-3, admm_tol = 1e-3, max_em_iter = 10, max_admm_iter = 100, 
                         penalty_type = "l1") {
  T <- length(X_list)
  p <- nrow(X_list[[1]])
  
  Theta_list <- lapply(1:T, function(t) diag(1, p))
  Z_list <- Theta_list
  U_list <- lapply(1:T, function(t) matrix(0, p, p))
  
  S_list <- em_estimate_S(X_list, Theta_list)
  
  for (em_iter in 1:max_em_iter) {
    Theta_old_list <- Theta_list
    
    # M-step: ADMM updates
    for (admm_iter in 1:max_admm_iter) {
      # Update Theta
      for (t in 1:T) {
        prev_theta <- if (t == 1) NULL else Theta_list[[t-1]]
        Theta_list[[t]] <- update_theta_huber(S_list[[t]], Z_list[[t]], U_list[[t]], prev_theta, rho, beta, delta)
      }
      
      # Update Z
      for (t in 1:T) {
        V_mat <- Theta_list[[t]] + U_list[[t]]
        if (penalty_type == "l1") {
          Z_list[[t]] <- soft_threshold(V_mat, lambda / rho)
        } else if (penalty_type == "l2") {
          Z_list[[t]] <- V_mat / (1 + lambda / rho)
        } else {
          stop("Invalid penalty_type. Use 'l1' or 'l2'.")
        }
        diag(Z_list[[t]]) <- diag(V_mat)  # don't penalize diagonals
      }
      
      # Update U
      for (t in 1:T) {
        U_list[[t]] <- U_list[[t]] + Theta_list[[t]] - Z_list[[t]]
      }
      
      # Check ADMM convergence
      primal_resid <- sum(sapply(1:T, function(t) norm(Theta_list[[t]] - Z_list[[t]], type = "F")^2))
      dual_resid <- sum(sapply(1:T, function(t) norm(Z_list[[t]] - Theta_old_list[[t]], type = "F")^2))
      if (sqrt(primal_resid) < admm_tol && sqrt(dual_resid) < admm_tol) {
        break
      }
    }
    
    # E-step: update S
    S_list <- em_estimate_S(X_list, Theta_list)
    
    # Check EM convergence
    em_change <- sum(sapply(1:T, function(t) norm(Theta_list[[t]] - Theta_old_list[[t]], type = "F")^2))
    if (sqrt(em_change) < em_tol) {
      break
    }
  }
  
  return(Theta_list)
}

#############
#
# tvgl_function <- function(data_list, lambda, beta, penalty_type = "l1", rho = 1, max_iter = 100, tol = 1e-3) {
#   T <- length(data_list)
#   p <- nrow(data_list[[1]])
#   
#   # Fallback if invalid penalty
#   if (!(penalty_type %in% c("l1", "l2"))) {
#     warning("TVGL only supports 'l1' or 'l2'. Defaulting to 'l1'.")
#     penalty_type <- "l1"
#   }
#   
#   # Initialize
#   Theta_list <- lapply(1:T, function(t) diag(p))
#   Z_list <- lapply(1:T, function(t) diag(p))
#   U_list <- lapply(1:T, function(t) matrix(0, p, p))
#   
#   # Precompute sample covariance matrices with small ridge
#   S_list <- lapply(data_list, function(X) {
#     S <- cov(X, use = "pairwise.complete.obs")
#     S + diag(1e-4, p)
#   })
#   
#   for (iter in 1:max_iter) {
#     Theta_prev <- Theta_list
#     
#     # Theta update
#     for (t in 1:T) {
#       A <- S_list[[t]] + rho * (Z_list[[t]] - U_list[[t]])
#       eig <- eigen(A, symmetric = TRUE)
#       D <- eig$values
#       V <- eig$vectors
#       D_new <- (D + sqrt(D^2 + 4 * rho)) / (2 * rho)
#       Theta_list[[t]] <- V %*% diag(D_new) %*% t(V)
#     }
#     
#     # Z update
#     for (t in 1:T) {
#       V_mat <- Theta_list[[t]] + U_list[[t]]
#       if (penalty_type == "l1") {
#         Z_soft <- sign(V_mat) * pmax(abs(V_mat) - lambda / rho, 0)
#         Z_list[[t]] <- Z_soft
#       } else if (penalty_type == "l2") {
#         Z_list[[t]] <- V_mat / (1 + lambda / rho)
#       }
#       diag(Z_list[[t]]) <- diag(V_mat)
#     }
#     
#     # U update
#     for (t in 1:T) {
#       U_list[[t]] <- U_list[[t]] + Theta_list[[t]] - Z_list[[t]]
#     }
#     
#     # Check convergence
#     primal_resid <- max(sapply(1:T, function(t) norm(Theta_list[[t]] - Z_list[[t]], type = "F")))
#     if (primal_resid < tol) break
#   }
#   
#   return(Theta_list)
# }

#new version tvgl_function:
tvgl_function <- function(data_list, lambda, beta, penalty_type = "l1", rho = 1, max_iter = 100, tol = 1e-3) {
  T <- length(data_list)
  p <- ncol(data_list[[1]])  # Use number of features
  
  # Validate penalty type
  if (!(penalty_type %in% c("l1", "l2"))) {
    warning("TVGL only supports 'l1' or 'l2'. Defaulting to 'l1'.")
    penalty_type <- "l1"
  }
  
  # Impute missing values in each time point
  data_list <- lapply(data_list, function(X) {
    for (j in 1:ncol(X)) {
      col_mean <- mean(X[, j], na.rm = TRUE)
      X[is.na(X[, j]), j] <- col_mean
    }
    return(X)
  })
  
  # Compute regularized covariance matrices
  S_list <- lapply(data_list, function(X) {
    S <- cov(X, use = "pairwise.complete.obs")
    
    # Ensure S is valid
    if (!is.matrix(S) || !all(dim(S) == c(p, p)) || any(is.na(S)) || any(is.infinite(S))) {
      S <- diag(1e-4, p)
    }
    
    # Regularize and symmetrize
    S <- S + diag(1e-4, p)
    S <- (S + t(S)) / 2
    
    # Ensure positive definite
    eig <- eigen(S)
    eig$values[eig$values < 1e-6] <- 1e-6
    S <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    
    return(S)
  })
  
  # Initialize variables
  Theta_list <- lapply(1:T, function(t) diag(p))
  Z_list <- lapply(1:T, function(t) diag(p))
  U_list <- lapply(1:T, function(t) matrix(0, p, p))
  
  # ADMM iterations
  for (iter in 1:max_iter) {
    Theta_prev <- Theta_list
    
    # Theta update
    for (t in 1:T) {
      A <- S_list[[t]] + rho * (Z_list[[t]] - U_list[[t]])
      eig <- eigen(A, symmetric = TRUE)
      D <- eig$values
      V <- eig$vectors
      D_new <- (D + sqrt(D^2 + 4 * rho)) / (2 * rho)
      Theta_list[[t]] <- V %*% diag(D_new) %*% t(V)
    }
    
    # Z update
    for (t in 1:T) {
      V_mat <- Theta_list[[t]] + U_list[[t]]
      if (penalty_type == "l1") {
        Z_list[[t]] <- sign(V_mat) * pmax(abs(V_mat) - lambda / rho, 0)
      } else if (penalty_type == "l2") {
        Z_list[[t]] <- V_mat / (1 + lambda / rho)
      }
      diag(Z_list[[t]]) <- diag(V_mat)
    }
    
    # U update
    for (t in 1:T) {
      U_list[[t]] <- U_list[[t]] + Theta_list[[t]] - Z_list[[t]]
    }
    
    # Check convergence
    primal_resid <- max(sapply(1:T, function(t) norm(Theta_list[[t]] - Z_list[[t]], "F")))
    if (primal_resid < tol) break
  }
  
  return(Theta_list)
}



##############
# rmtvgl_function: Robust Missing-data-aware TVGL


rmtvgl_function <- function(data_list, lambda, beta, delta = 1.0, rho = 1.0,
                            em_tol = 1e-3, admm_tol = 1e-3, max_em_iter = 5, max_admm_iter = 10,
                            penalty_type = "l1", alpha = 0.5) {
  T <- length(data_list)
  p <- ncol(data_list[[1]])
  
  # === Robust covariance + NA imputation ===
  em_estimate_S <- function(X_list) {
    S_list <- list()
    new_X_list <- list()
    for (i in 1:length(X_list)) {
      X <- X_list[[i]]
      for (j in 1:ncol(X)) {
        col_mean <- mean(X[, j], na.rm = TRUE)
        if (is.nan(col_mean)) col_mean <- 0  # fallback if all values are NA
        X[is.na(X[, j]), j] <- col_mean
      }
      S <- tryCatch(cov(X), error = function(e) diag(1e-4, p))
      S <- S + diag(1e-4, p)
      S <- (S + t(S)) / 2
      S_list[[i]] <- S
      new_X_list[[i]] <- X
    }
    return(list(S_list = S_list, X_list = new_X_list))
  }
  
  update_theta_huber <- function(S, Z, U, Theta_prev, rho, beta, delta) {
    Theta <- diag(nrow(S))
    for (iter in 1:10) {
      R <- Theta - S
      grad <- matrix(0, nrow = nrow(R), ncol = ncol(R))
      small_idx <- abs(R) <= delta
      large_idx <- !small_idx
      grad[small_idx] <- R[small_idx]
      grad[large_idx] <- delta * sign(R[large_idx])
      G <- grad + rho * (Theta - Z + U)
      if (!is.null(Theta_prev)) {
        G <- G + 2 * beta * (Theta - Theta_prev)
      }
      Theta_new <- Theta - 0.5 * G
      if (norm(Theta_new - Theta, "F") < 1e-4) break
      Theta <- Theta_new
    }
    return((Theta + t(Theta)) / 2)
  }
  
  # === Initial E-step ===
  em_out <- em_estimate_S(data_list)
  S_list <- em_out$S_list
  data_list <- em_out$X_list  # update with imputed data
  
  # === Initialize ADMM variables ===
  Theta_list <- lapply(1:T, function(t) diag(p))
  Z_list <- Theta_list
  U_list <- lapply(1:T, function(t) matrix(0, p, p))
  
  for (em_iter in 1:max_em_iter) {
    cat(sprintf("EM iter %d\n", em_iter))
    Theta_old_list <- Theta_list
    
    for (admm_iter in 1:max_admm_iter) {
      for (t in 1:T) {
        Theta_prev <- if (t == 1) NULL else Theta_list[[t - 1]]
        Theta_list[[t]] <- update_theta_huber(S_list[[t]], Z_list[[t]], U_list[[t]],
                                              Theta_prev, rho, beta, delta)
      }
      for (t in 1:T) {
        V <- Theta_list[[t]] + U_list[[t]]
        if (penalty_type == "l1") {
          Z_list[[t]] <- sign(V) * pmax(abs(V) - lambda / rho, 0)
        } else if (penalty_type == "l2") {
          Z_list[[t]] <- V / (1 + lambda / rho)
        } else if (penalty_type == "elastic") {
          l1_part <- alpha * lambda
          l2_part <- (1 - alpha) * lambda
          Z_temp <- V / (1 + l2_part / rho)
          Z_list[[t]] <- sign(Z_temp) * pmax(abs(Z_temp) - l1_part / rho, 0)
        }
        diag(Z_list[[t]]) <- diag(V)
      }
      for (t in 1:T) {
        U_list[[t]] <- U_list[[t]] + Theta_list[[t]] - Z_list[[t]]
      }
    }
    
    # === Update S and imputed data for next EM step ===
    em_out <- em_estimate_S(data_list)
    S_list <- em_out$S_list
    data_list <- em_out$X_list
  }
  
  return(Theta_list)
}



##################
#3.synthetic data generation
generate_synthetic_data <- function(T = 100, p = 10, type = "global", noise_sd = 0.05) {
  #set.seed(123)
  
  # Helper: Generate a random sparse positive definite matrix
  generate_sparse_precision <- function(p, sparsity = 0.05) {
    repeat {
      mat <- matrix(0, p, p)
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          if (runif(1) < sparsity) {
            value <- runif(1, min = 0.2, max = 0.6) * sample(c(-1, 1), 1)
            mat[i, j] <- value
            mat[j, i] <- value
          }
        }
      }
      diag(mat) <- abs(rowSums(mat)) + 0.5 # Make it diagonally dominant
      
      if (all(eigen(mat)$values > 0)) {
        return(mat)
      }
    }
  }
  
  Theta_list <- list()
  
  if (type == "global") {
    Theta1 <- generate_sparse_precision(p)
    Theta2 <- generate_sparse_precision(p)
    change_point <- T %/% 2
    
    for (t in 1:T) {
      if (t <= change_point) {
        Theta_list[[t]] <- Theta1
      } else {
        Theta_list[[t]] <- Theta2
      }
    }
    
  } else if (type == "local") {
    Theta_base <- generate_sparse_precision(p)
    for (t in 1:T) {
      Theta_t <- Theta_base
      if (t == (T %/% 2)) {
        # Local shift: modify only a small submatrix
        idx <- sample(1:p, 3) # randomly pick 3 nodes
        Theta_t[idx, idx] <- Theta_t[idx, idx] + diag(runif(3, 0.5, 1.0))
      }
      Theta_list[[t]] <- Theta_t
    }
    
  } else {
    stop("Unknown type. Choose 'global' or 'local'.")
  }
  
  # Generate synthetic data
  data_list <- list()
  for (t in 1:T) {
    Sigma <- solve(Theta_list[[t]])
    Sigma <- (Sigma + t(Sigma)) / 2  # symmetrize
    
    # If still not perfect, adjust slightly
    eig <- eigen(Sigma)
    eig$values[eig$values < 1e-6] <- 1e-6
    Sigma <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    
    samples <- MASS::mvrnorm(n = 10, mu = rep(0, p), Sigma = Sigma)
    samples <- samples + matrix(rnorm(length(samples), sd = noise_sd), nrow = nrow(samples))
    
    data_list[[t]] <- samples
  }
  
  return(list(data = data_list, Theta = Theta_list))
}

################
## evaluation_metrics.R

library(pROC)
library(ggplot2)

# record runtime
time_run <- function(expr) {
  start_time <- Sys.time()
  result <- eval(expr)
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  return(list(result = result, runtime = runtime))
}

# Metric functions
evaluate_f1 <- function(true_matrix, est_matrix) {
  true_edges <- (true_matrix != 0)
  est_edges <- (est_matrix != 0)
  tp <- sum(true_edges & est_edges)
  fp <- sum(!true_edges & est_edges)
  fn <- sum(true_edges & !est_edges)
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  if (precision + recall == 0) return(0)
  return(2 * precision * recall / (precision + recall))
}

# evaluate_auc <- function(true_matrix, est_matrix) {
#   true_labels <- as.vector(true_matrix != 0)
#   scores <- as.vector(abs(est_matrix))
#   roc_obj <- roc(true_labels, scores)
#   return(auc(roc_obj))
# }
evaluate_auc <- function(true_matrix, est_matrix) {
  true_labels <- as.vector(true_matrix != 0)
  scores <- as.vector(abs(est_matrix))
  
  # Only compute AUC if there are both TRUE and FALSE in true_labels
  if (length(unique(true_labels)) < 2) return(NA)
  
  roc_obj <- roc(true_labels, scores)
  return(auc(roc_obj))
}

evaluate_mac <- function(est_list) {
  changes <- sapply(2:length(est_list), function(t) {
    mean(abs(as.matrix(est_list[[t]]) - as.matrix(est_list[[t - 1]])))
  })
  mean(changes)
}

evaluate_edge_stability <- function(est_list, threshold = 1e-6) {
  T <- length(est_list)
  p <- nrow(est_list[[1]])
  binary_list <- lapply(est_list, function(mat) abs(as.matrix(mat)) > threshold)
  edge_changes <- sum(sapply(2:T, function(t) sum(binary_list[[t]] != binary_list[[t - 1]])))
  total_edges <- (T - 1) * p^2
  return(edge_changes / total_edges)
}

# Model comparison
compare_methods <- function(true_list, data_list, method_tvgl, method_rmtvgl,
                            lambda, beta, penalty_type = "l1") {
  T <- length(data_list)
  
  # Run TVGL
  tvgl_out <- time_run({ method_tvgl(data_list, lambda, beta, penalty_type) })
  est_tvgl <- tvgl_out$result
  runtime_tvgl <- tvgl_out$runtime
  
  # Run RM-TVGL
  rmtvgl_out <- time_run({ method_rmtvgl(data_list, lambda, beta, delta = 1.0, penalty_type = penalty_type) })
  est_rmtvgl <- rmtvgl_out$result
  runtime_rmtvgl <- rmtvgl_out$runtime
  
  # Metric helper function with NA-safe AUC
  compute_metrics <- function(est_list) {
    f1 <- mean(sapply(1:T, function(t) evaluate_f1(as.matrix(true_list[[t]]), as.matrix(est_list[[t]]))))
    aucs <- sapply(1:T, function(t) evaluate_auc(as.matrix(true_list[[t]]), as.matrix(est_list[[t]])))
    auc_mean <- mean(na.omit(aucs))  # Skip NAs due to 1-class true labels
    mac <- evaluate_mac(est_list)
    stab <- evaluate_edge_stability(est_list)
    return(list(F1 = f1, AUC = auc_mean, MAC = mac, Stability = stab))
  }
  
  results <- list(
    TVGL = c(compute_metrics(est_tvgl), Runtime = runtime_tvgl),
    RM_TVGL = c(compute_metrics(est_rmtvgl), Runtime = runtime_rmtvgl)
  )
  
  return(results)
}

# Wrapper to evaluate RM-TVGL model
evaluate_model <- function(true_list, data_list, my_penalty_type, alpha = 0.5, lambda = 0.5, beta = 1) {
  res <- compare_methods(
    true_list, data_list,
    method_tvgl = tvgl_function,
    method_rmtvgl = function(data_list, lambda, beta, ...) {
      rmtvgl_function(data_list = data_list, lambda = lambda, beta = beta,
                      delta = 1.0, penalty_type = my_penalty_type, alpha = alpha)
    },
    lambda = lambda, beta = beta,
    penalty_type = my_penalty_type
  )
  return(res$RM_TVGL)
}



result_l1 <- evaluate_model(true_list, data_list, "l1")
result_l2 <- evaluate_model(true_list, data_list, "l2")
result_elastic <- evaluate_model(true_list, data_list, "elastic", alpha = 0.5)

# Wrapper to run multiple times and average
run_comparison_multiple <- function(n_runs = 5, T = 100, p = 10, type = "global", 
                                    lambda = 0.5, beta = 1,
                                    method_tvgl, method_rmtvgl, penalty_type = "l1") {
  library(dplyr)
  
  metric_names <- c("F1", "AUC", "MAC", "Stability", "Runtime")
  results_list <- list(
    TVGL = lapply(metric_names, function(x) c()),
    RM_TVGL = lapply(metric_names, function(x) c())
  )
  names(results_list$TVGL) <- metric_names
  names(results_list$RM_TVGL) <- metric_names
  
  for (i in 1:n_runs) {
    cat(sprintf("Run %d/%d\n", i, n_runs))
    synthetic <- generate_synthetic_data(T = T, p = p, type = type)
    true_list <- synthetic$Theta
    data_list <- synthetic$data
    
    out <- compare_methods(true_list, data_list, method_tvgl, method_rmtvgl,
                           lambda, beta, penalty_type)
    
    for (metric in metric_names) {
      results_list$TVGL[[metric]] <- c(results_list$TVGL[[metric]], out$TVGL[[metric]])
      results_list$RM_TVGL[[metric]] <- c(results_list$RM_TVGL[[metric]], out$RM_TVGL[[metric]])
    }
  }
  
  summarize <- function(x) sprintf("%.3f ± %.3f", mean(x), sd(x))
  
  summary_df <- data.frame(
    Method = c("TVGL", "RM-TVGL"),
    F1 = c(summarize(results_list$TVGL$F1), summarize(results_list$RM_TVGL$F1)),
    AUC = c(summarize(results_list$TVGL$AUC), summarize(results_list$RM_TVGL$AUC)),
    MAC = c(summarize(results_list$TVGL$MAC), summarize(results_list$RM_TVGL$MAC)),
    Stability = c(summarize(results_list$TVGL$Stability), summarize(results_list$RM_TVGL$Stability)),
    Runtime_sec = c(summarize(results_list$TVGL$Runtime), summarize(results_list$RM_TVGL$Runtime))
  )
  
  return(summary_df)
}


############## 05/15/2025 example code to get L1 ,L2,elastic for 10d-synthetic dataset, 1 run:
# Assumes all your previous helper functions are defined:
# generate_synthetic_data, tvgl_function, rmtvgl_function, time_run, evaluate_f1, evaluate_auc, evaluate_mac, evaluate_edge_stability

# Step 1: Generate synthetic data
set.seed(123)
T <- 100
p <- 10
synthetic_data <- generate_synthetic_data(T = T, p = p, type = "global")
true_list <- synthetic_data$Theta
data_list <- synthetic_data$data

# Add 5% missing values and 5% outliers
for (t in 1:T) {
  X <- data_list[[t]]
  mask_na <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
  X[mask_na] <- NA
  mask_outlier <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
  X[mask_outlier] <- X[mask_outlier] * 10
  data_list[[t]] <- X
}

# Step 2: Define parameters
lambda <- 0.5
beta <- 1
delta <- 1.0
alpha_elastic <- 0.5

penalties <- c("l1", "l2", "elastic")
methods <- c("TVGL", "RM-TVGL")

# Step 3: Run experiments
results_list <- list()

for (penalty in penalties) {
  cat("Running penalty:", penalty, "\n")
  
  # TVGL
  out_tvgl <- time_run({
    tvgl_function(data_list, lambda = lambda, beta = beta, penalty_type = ifelse(penalty == "elastic", "l1", penalty))
  })
  est_tvgl <- out_tvgl$result
  time_tvgl <- out_tvgl$runtime
  
  # RM-TVGL
  out_rmtvgl <- time_run({
    rmtvgl_function(data_list, lambda = lambda, beta = beta, delta = delta,
                    penalty_type = penalty, alpha = alpha_elastic)
  })
  est_rmtvgl <- out_rmtvgl$result
  time_rmtvgl <- out_rmtvgl$runtime
  
  # Metrics
  for (method in methods) {
    est <- if (method == "TVGL") est_tvgl else est_rmtvgl
    time <- if (method == "TVGL") time_tvgl else time_rmtvgl
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Method = method,
      Penalty = penalty,
      F1 = mean(sapply(1:T, function(t) evaluate_f1(true_list[[t]], est[[t]]))),
      AUC = mean(sapply(1:T, function(t) evaluate_auc(true_list[[t]], est[[t]]))),
      MAC = evaluate_mac(est),
      Stability = evaluate_edge_stability(est),
      Runtime_sec = time
    )
  }
}

# Combine results
library(dplyr)
summary_df <- bind_rows(results_list)

# Print results
print(summary_df)

