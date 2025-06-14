#### here are the all example codes of RMTVGL:

##### 1. RM-TVGL 05/15/2025  5 runs, example code, synthetic dataset(10dimentional), all penalties(Table 1):
library(Matrix)
library(MASS)
library(pROC)
library(dplyr)

# Function to summarize with mean ± sd
summarize_metric <- function(x) sprintf("%.3f ± %.3f", mean(x), sd(x))

# some Parameters
T <- 100
p <- 10
lambda <- 0.5
beta <- 1
delta <- 1.0
alpha_elastic <- 0.5
penalties <- c("l1", "l2", "elastic")
methods <- c("TVGL", "RM-TVGL")
n_repeats <- 5

# Initialize result storage
all_results <- list()

for (rep in 1:n_repeats) {
  cat("Running repetition", rep, "\n")
  
  # Generate synthetic data
  set.seed(100 + rep)
  synthetic_data <- generate_synthetic_data(T = T, p = p, type = "global")
  true_list <- synthetic_data$Theta
  data_list <- synthetic_data$data
  
  # Add missing + outlier values (5% each)
  for (t in 1:T) {
    X <- data_list[[t]]
    mask_na <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
    X[mask_na] <- NA
    mask_outlier <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
    X[mask_outlier] <- X[mask_outlier] * 10
    data_list[[t]] <- X
  }
  
  # Run across all penalties
  for (penalty in penalties) {
    cat("  Penalty:", penalty, "\n")
    
    # TVGL
    out_tvgl <- time_run({
      tvgl_function(data_list, lambda = lambda, beta = beta,
                    penalty_type = ifelse(penalty == "elastic", "l1", penalty))
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
    
    # Store metrics
    for (method in methods) {
      est <- if (method == "TVGL") est_tvgl else est_rmtvgl
      time <- if (method == "TVGL") time_tvgl else time_rmtvgl
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Method = method,
        Penalty = penalty,
        Rep = rep,
        F1 = mean(sapply(1:T, function(t) evaluate_f1(true_list[[t]], est[[t]]))),
        AUC = mean(sapply(1:T, function(t) evaluate_auc(true_list[[t]], est[[t]]))),
        MAC = evaluate_mac(est),
        Stability = evaluate_edge_stability(est),
        Runtime_sec = time
      )
    }
  }
}

# Combine all repetitions
df_all <- bind_rows(all_results)

# Summarize with mean ± sd
summary_df <- df_all %>%
  group_by(Method, Penalty) %>%
  summarise(
    F1 = summarize_metric(F1),
    AUC = summarize_metric(AUC),
    MAC = summarize_metric(MAC),
    Stability = summarize_metric(Stability),
    Runtime_sec = summarize_metric(Runtime_sec),
    .groups = "drop"
  )

# Print final result
print(summary_df)
write.csv(summary_df, "RM_TVGL_10gene synthetic_summary.csv", row.names = FALSE)

################
#2.10gene synthetic data the plot for different delta value of huber loss(Fig1):
#curves of auc ,f1,mac,stability fr different delta:
# Load required libraries
library(Matrix)
library(MASS)
library(pROC)
library(dplyr)
library(tidyr)
library(ggplot2)

# Parameters
T <- 100
p <- 10
lambda <- 0.5
beta <- 1
alpha_elastic <- 0.5
penalties <- c("elastic", "l1", "l2")
delta_list <- c(0.5, 1, 2, 5)
n_runs <- 5

# Placeholder for results
all_results <- list()

# Function to run one RM-TVGL setting
run_one_setting <- function(data_list, true_list, lambda, beta, delta, penalty) {
  result <- rmtvgl_function(data_list,
                            lambda = lambda,
                            beta = beta,
                            delta = delta,
                            penalty_type = penalty,
                            alpha = alpha_elastic)
  list(
    F1 = mean(sapply(1:T, function(t) evaluate_f1(true_list[[t]], result[[t]]))),
    AUC = mean(sapply(1:T, function(t) evaluate_auc(true_list[[t]], result[[t]]))),
    MAC = evaluate_mac(result),
    Stability = evaluate_edge_stability(result)
  )
}

# Run all settings 
for (penalty in penalties) {
  for (delta in delta_list) {
    cat("Running penalty =", penalty, ", delta =", delta, "\n")
    
    f1_vec <- auc_vec <- mac_vec <- stab_vec <- numeric(n_runs)
    
    for (run in 1:n_runs) {
      set.seed(100 + run)
      
      synthetic_data <- generate_synthetic_data(T = T, p = p, type = "global")
      true_list <- synthetic_data$Theta
      data_list <- synthetic_data$data
      
      for (t in 1:T) {
        X <- data_list[[t]]
        mask_na <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
        X[mask_na] <- NA
        mask_outlier <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
        X[mask_outlier] <- X[mask_outlier] * 10
        data_list[[t]] <- X
      }
      
      res <- run_one_setting(data_list, true_list, lambda, beta, delta, penalty)
      f1_vec[run] <- res$F1
      auc_vec[run] <- res$AUC
      mac_vec[run] <- res$MAC
      stab_vec[run] <- res$Stability
    }
    
    all_results[[length(all_results) + 1]] <- data.frame(
      Penalty = penalty,
      Delta = delta,
      F1 = mean(f1_vec),
      AUC = mean(auc_vec),
      MAC = mean(mac_vec),
      Stability = mean(stab_vec)
    )
  }
}

# Prepare data for plotting
df_all <- bind_rows(all_results)
df_long <- df_all %>%
  pivot_longer(cols = c(F1, AUC, MAC, Stability), names_to = "Metric", values_to = "Value")

# --- Plot: one file per metric, no title/caption ---
color_map <- c("elastic" = "#E64B35", "l1" = "#4DBBD5", "l2" = "#000000")
line_style <- c("elastic" = "dashed", "l1" = "dashed", "l2" = "dashed")

metrics <- c("F1", "AUC", "MAC", "Stability")

for (metric in metrics) {
  p <- ggplot(df_long %>% filter(Metric == metric),
              aes(x = Delta, y = Value, color = Penalty, group = Penalty)) +
    geom_line(size = 1.2, aes(linetype = Penalty)) +
    geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
    scale_x_continuous(trans = "log2", breaks = delta_list) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = line_style) +
    theme_minimal(base_size = 16) +
    labs(x = expression(delta), y = metric) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
  
  print(p)
  # Optional: save for Overleaf
  # ggsave(paste0("Figure_", metric, "_vs_delta.pdf"), p, width = 6, height = 4)
}


########## 3. 20 gene synthetic dataset and output(I didn't put this result in the paper, since I think one synthetic dataset is enough):
library(Matrix)
library(MASS)
library(pROC)
library(dplyr)

# --- AUC FUNCTION: Original version you used ---
evaluate_auc <- function(true_mat, est_mat) {
  true_vec <- as.vector(true_mat)
  est_vec <- as.vector(est_mat)
  
  idx <- which(row(true_mat) != col(true_mat))
  true_edges <- true_vec[idx]
  est_scores <- est_vec[idx]
  
  roc_obj <- pROC::roc(response = true_edges, predictor = est_scores, quiet = TRUE)
  return(as.numeric(pROC::auc(roc_obj)))
}

# Summarizer 
summarize_metric <- function(x) sprintf("%.3f ± %.3f", mean(x), sd(x))

# Parameters 
T <- 100
p <- 20
lambda <- 0.5
beta <- 1
delta <- 1.0
alpha_elastic <- 0.5
penalties <- c("l1", "l2", "elastic")
methods <- c("TVGL", "RM-TVGL")
n_repeats <- 5
all_results <- list()

# loop
for (rep in 1:n_repeats) {
  cat("Running repetition", rep, "\n")
  set.seed(100 + rep)
  
  # Generate synthetic data
  synthetic_data <- generate_synthetic_data(T = T, p = p, type = "global")
  true_list <- synthetic_data$Theta
  data_list <- synthetic_data$data
  
  # Add 5% missing and 5% outliers
  for (t in 1:T) {
    X <- data_list[[t]]
    mask_na <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
    X[mask_na] <- NA
    mask_outlier <- matrix(runif(length(X)) < 0.05, nrow = nrow(X))
    X[mask_outlier] <- X[mask_outlier] * 10
    data_list[[t]] <- X
  }
  
  for (penalty in penalties) {
    cat("  Penalty:", penalty, "\n")
    
    #  TVGL 
    out_tvgl <- time_run({
      tvgl_function(data_list, lambda = lambda, beta = beta,
                    penalty_type = ifelse(penalty == "elastic", "l1", penalty))
    })
    est_tvgl <- out_tvgl$result
    time_tvgl <- out_tvgl$runtime
    
    #  RM-TVGL 
    out_rmtvgl <- time_run({
      rmtvgl_function(data_list, lambda = lambda, beta = beta, delta = delta,
                      penalty_type = penalty, alpha = alpha_elastic)
    })
    est_rmtvgl <- out_rmtvgl$result
    time_rmtvgl <- out_rmtvgl$runtime
    
    #  Store metrics 
    for (method in methods) {
      est <- if (method == "TVGL") est_tvgl else est_rmtvgl
      time <- if (method == "TVGL") time_tvgl else time_rmtvgl
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Method = method,
        Penalty = penalty,
        Rep = rep,
        F1 = mean(sapply(1:T, function(t) evaluate_f1(true_list[[t]], est[[t]]))),
        AUC = mean(sapply(1:T, function(t) evaluate_auc(true_list[[t]], est[[t]]))),
        MAC = evaluate_mac(est),
        Stability = evaluate_edge_stability(est),
        Runtime_sec = time
      )
    }
  }
}

#  Final summary 
df_all <- bind_rows(all_results)

summary_df <- df_all %>%
  group_by(Method, Penalty) %>%
  summarise(
    F1 = summarize_metric(F1),
    AUC = summarize_metric(AUC),
    MAC = summarize_metric(MAC),
    Stability = summarize_metric(Stability),
    Runtime_sec = summarize_metric(Runtime_sec),
    .groups = "drop"
  )

# Output results
print(summary_df)
write.csv(summary_df, "Synthetic20_TVGL_RMTVGL_summary.csv", row.names = FALSE)







#######4. real dataset part: 05/15/2025 collin cancer real dataset :
# Required packages
install.packages("e1071")
install.packages("Matrix")

# Load necessary packages
library(e1071)
library(Matrix)

# Read LIBSVM file
data <- read.matrix.csr("colon-cancer")  # or provide full path

# Convert sparse matrix to dense
X <- as.matrix(data$x)  # This works because data$x is a 'matrix.csr' (SparseM)

# Response vector
y <- data$y

# Check dimensions
cat("Feature matrix shape:", dim(X), "\n")
cat("Label vector length:", length(y), "\n")

# Split into T time points
T <- 10  # e.g., 10 time steps
samples_per_time <- floor(nrow(X) / T)

data_list <- list()
for (t in 1:T) {
  idx <- ((t - 1) * samples_per_time + 1):(t * samples_per_time)
  data_list[[t]] <- X[idx, ]
}

# Handle the remainder (if 62 not divisible by T)
if ((T * samples_per_time) < nrow(X)) {
  data_list[[T]] <- rbind(data_list[[T]], X[(T * samples_per_time + 1):nrow(X), ])
}

# Check shape
sapply(data_list, dim)

theta_tvgl <- tvgl_function(data_list, lambda = 0.5, beta = 1, penalty_type = "l1")
theta_rmtvgl <- rmtvgl_function(data_list, lambda = 0.5, beta = 1, delta = 1.0, penalty_type = "l1")
f1_tvgl <- mean(sapply(1:length(data_list), function(t) evaluate_f1(theta_tvgl[[t]], theta_tvgl[[t]])))  # placeholder
auc_tvgl <- mean(sapply(1:length(data_list), function(t) evaluate_auc(theta_tvgl[[t]], theta_tvgl[[t]])))




#### real data colon cancer all output using 3 penalties, since the dataset is too large, we can pick the top 100 most variable genes(Table 2):
library(Matrix)
library(MASS)
library(pROC)
library(dplyr)
library(e1071)

# Helper functions 
summarize_metric <- function(x) sprintf("%.3f ± %.3f", mean(x), sd(x))

#  Parameters 
T <- 10
lambda <- 0.5
beta <- 1
delta <- 1.0
alpha_elastic <- 0.5
penalties <- c("l1", "l2", "elastic")
methods <- c("TVGL", "RM-TVGL")
n_repeats <- 5

#  Load Colon Cancer Data 
data <- read.matrix.csr("colon-cancer")  # LIBSVM format
X_full <- as.matrix(data$x)              # 62 × 2000
y <- data$y

# Select Top 100 Most Variable Genes
var_rank <- order(apply(X_full, 2, var), decreasing = TRUE)
X <- X_full[, var_rank[1:100]]

#  Result storage
all_results <- list()

for (rep in 1:n_repeats) {
  cat("Running repetition", rep, "\n")
  set.seed(100 + rep)
  
  # Split into T time points
  samples_per_time <- floor(nrow(X) / T)
  data_list <- list()
  for (t in 1:T) {
    idx <- ((t - 1) * samples_per_time + 1):(t * samples_per_time)
    data_list[[t]] <- X[idx, , drop = FALSE]
  }
  # Handle leftovers
  if ((T * samples_per_time) < nrow(X)) {
    data_list[[T]] <- rbind(data_list[[T]], X[(T * samples_per_time + 1):nrow(X), , drop = FALSE])
  }
  
  # Add 5% missing values + outliers
  for (t in 1:T) {
    X_t <- data_list[[t]]
    mask_na <- matrix(runif(length(X_t)) < 0.05, nrow = nrow(X_t))
    X_t[mask_na] <- NA
    mask_outlier <- matrix(runif(length(X_t)) < 0.05, nrow = nrow(X_t))
    X_t[mask_outlier] <- X_t[mask_outlier] * 10
    data_list[[t]] <- X_t
  }
  
  # Ground truth using baseline TVGL (for F1, AUC)
  true_list <- tvgl_function(data_list, lambda = 5, beta = 0.5, penalty_type = "l1")
  true_list_sparse <- lapply(true_list, function(mat) {
    mat[abs(mat) < 0.05] <- 0
    mat
  })
  
  # Loop through penalties
  for (penalty in penalties) {
    cat("  Penalty:", penalty, "\n")
    
    # TVGL
    out_tvgl <- time_run({
      tvgl_function(data_list, lambda = lambda, beta = beta,
                    penalty_type = ifelse(penalty == "elastic", "l1", penalty))
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
    
    # Store metrics
    for (method in methods) {
      est <- if (method == "TVGL") est_tvgl else est_rmtvgl
      time <- if (method == "TVGL") time_tvgl else time_rmtvgl
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Method = method,
        Penalty = penalty,
        Rep = rep,
        F1 = mean(sapply(1:T, function(t) evaluate_f1(true_list_sparse[[t]], est[[t]]))),
        AUC = mean(sapply(1:T, function(t) evaluate_auc(true_list_sparse[[t]], est[[t]]))),
        MAC = evaluate_mac(est),
        Stability = evaluate_edge_stability(est),
        Runtime_sec = time
      )
    }
  }
}

# Combine and summarize
df_all <- bind_rows(all_results)

summary_df <- df_all %>%
  group_by(Method, Penalty) %>%
  summarise(
    F1 = summarize_metric(F1),
    AUC = summarize_metric(AUC),
    MAC = summarize_metric(MAC),
    Stability = summarize_metric(Stability),
    Runtime_sec = summarize_metric(Runtime_sec),
    .groups = "drop"
  )

# Output
print(summary_df)
write.csv(summary_df, "colon_summary_results.csv", row.names = FALSE)

###### 5.leu real dataset all penalties output, we used top 200 genes(Table 3 ):
# Required Libraries
library(Matrix)
library(MASS)
library(dplyr)
library(pROC)

# Parameters
T <- 10
p <- 200
lambda <- 0.5
beta <- 1.0
delta <- 1.0
alpha_elastic <- 0.5
penalties <- c("l1", "l2", "elastic")
methods <- c("TVGL", "RM-TVGL")
n_runs <- 5

# Initialize result storage
all_results <- list()

for (rep in 1:n_runs) {
  cat("Running repetition", rep, "\n")
  
  # Load and preprocess Leukemia
  data <- read.matrix.csr("leu")
  X_full <- as.matrix(data$x)
  var_rank <- order(apply(X_full, 2, var), decreasing = TRUE)
  X <- X_full[, var_rank[1:p]]
  
  # Split into T time points
  samples_per_time <- floor(nrow(X) / T)
  data_list <- list()
  for (t in 1:T) {
    idx <- ((t - 1) * samples_per_time + 1):(t * samples_per_time)
    data_list[[t]] <- X[idx, , drop = FALSE]
  }
  if ((T * samples_per_time) < nrow(X)) {
    data_list[[T]] <- rbind(data_list[[T]], X[(T * samples_per_time + 1):nrow(X), , drop = FALSE])
  }
  
  # Add 5% missing and outliers
  set.seed(100 + rep)
  for (t in 1:T) {
    X_t <- data_list[[t]]
    mask_na <- matrix(runif(length(X_t)) < 0.05, nrow = nrow(X_t))
    X_t[mask_na] <- NA
    mask_outlier <- matrix(runif(length(X_t)) < 0.05, nrow = nrow(X_t))
    X_t[mask_outlier] <- X_t[mask_outlier] * 10
    data_list[[t]] <- X_t
  }
  
  # Generate surrogate ground truth using TVGL
  true_list <- tvgl_function(data_list, lambda = 5, beta = 0.5, penalty_type = "l1")
  true_list_sparse <- lapply(true_list, function(mat) {
    mat[abs(mat) < 0.05] <- 0
    mat
  })
  
  for (penalty in penalties) {
    cat("  Penalty:", penalty, "\n")
    
    # TVGL
    out_tvgl <- time_run({
      tvgl_function(data_list, lambda = lambda, beta = beta,
                    penalty_type = ifelse(penalty == "elastic", "l1", penalty))
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
    
    # Store metrics
    for (method in methods) {
      est <- if (method == "TVGL") est_tvgl else est_rmtvgl
      time <- if (method == "TVGL") time_tvgl else time_rmtvgl
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Method = method,
        Penalty = penalty,
        Rep = rep,
        F1 = mean(sapply(1:T, function(t) evaluate_f1(true_list_sparse[[t]], est[[t]]))),
        AUC = mean(sapply(1:T, function(t) evaluate_auc(true_list_sparse[[t]], est[[t]]))),
        MAC = evaluate_mac(est),
        Stability = evaluate_edge_stability(est),
        Runtime_sec = time
      )
    }
  }
}

# Combine and summarize
df_all <- bind_rows(all_results)
summarize_metric <- function(x) sprintf("%.3f ± %.3f", mean(x), sd(x))

summary_df <- df_all %>%
  group_by(Method, Penalty) %>%
  summarise(
    F1 = summarize_metric(F1),
    AUC = summarize_metric(AUC),
    MAC = summarize_metric(MAC),
    Stability = summarize_metric(Stability),
    Runtime_sec = summarize_metric(Runtime_sec),
    .groups = "drop"
  )

# Print and save
print(summary_df)
write.csv(summary_df, "Leukemia_TVGL_RMTVGL_summary.csv", row.names = FALSE)

