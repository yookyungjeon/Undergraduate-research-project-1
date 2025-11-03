# Title: Graphical Lasso–guided MICE Imputation (p = 20, initial imputation: mean, diag(wi)=0 in voting)


# 1) Setup ------------------------------------------------
library(MASS)
library(glasso)
library(mice)
library(ggplot2)
library(pracma)
library(ROCR)
library(gridExtra)
library(grid)

n <- 100
p <- 20

time_results <- list()
roc_dfs <- list(p1 = list(), p2 = list(), p3 = list(), p4 = list())
auc_values <- data.frame(Method = character(), Run = integer(), AUC = numeric(), stringsAsFactors = FALSE)


# 2) Simulation runs --------------------------------------
for (run in 1:10) {

  print(paste("Starting run", run))
  set.seed(run * 100)

  # 2-1) Generate data
  # (i) Precision matrix (Omega): 1 on diagonal, 0.5 on immediate off-diagonals
  Omega <- diag(1, p)
  for (i in 1:(p - 1)) {
    Omega[i, i + 1] <- 0.5
    Omega[i + 1, i] <- 0.5
  }
  Sigma <- solve(Omega)

  # (ii) Generate multivariate normal data
  complete_data <- mvrnorm(n, rep(0, p), Sigma)

  # (iii) Inject MCAR missingness at rate 0.1
  missing_pattern <- matrix(rbinom(n * p, 1, 0.1), n, p)
  data_na <- complete_data
  data_na[missing_pattern == 1] <- NA


  # 2-2) Initial mean imputation 
  start_time <- Sys.time()

  initial_imputation <- function(data) {
    data_imputed <- data
    for (col in 1:ncol(data_imputed)) {
      if (is.numeric(data_imputed[, col])) {
        mean_value <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[is.na(data_imputed[, col]), col] <- mean_value
      }
    }
    data_imputed
  }

  data <- initial_imputation(data_na)

  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "Initial Imputation",
                                                time = difftime(end_time, start_time, units = "mins"))))


  # 2-3) Iterative MICE guided by glasso 
  start_time <- Sys.time()

  for (data_set_number in 1:5) {
    for (iter in 1:5) {
      # (i) Fit glasso path over rho grid
      rho_values <- seq(0, 200, by = 1)
      glasso_results <- glassopath(cov(data), rholist = rho_values, trace = 0)

      # (ii) Select rho via BIC-like criterion
      BIC_values <- sapply(seq_along(rho_values), function(j) {
        Omega_j <- glasso_results$wi[, , j]
        non_zero_count <- sum(abs(Omega_j) > 0)
        log(det(Omega_j)) - sum(diag(cov(data) %*% Omega_j)) - non_zero_count * log(n)
      })
      optimal_rho_index <- which.min(BIC_values)

      # (iii) Build predictorMatrix from nonzeros of selected precision
      glasso_result_wi <- glasso_results$wi[, , optimal_rho_index]
      predictorMatrix <- ifelse(abs(glasso_result_wi) > 0, 1, 0)
      diag(predictorMatrix) <- 0

      # (iv) One-step MICE using predictorMatrix; only impute where missing
      imputed_data <- mice(
        data,
        m = 1,
        maxit = 1,
        predictorMatrix = predictorMatrix,
        where = (missing_pattern == 1),
        seed = 500,
        printFlag = FALSE
      )
      completed_data <- complete(imputed_data, action = 1)

      # (v) Update data
      data <- completed_data
    }

    assign(paste("pseudo_complete", data_set_number, sep = ""), data, envir = .GlobalEnv)
  }

  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "Imputation",
                                                time = difftime(end_time, start_time, units = "mins"))))


  # 2-4) ROC evaluation vs. true Omega 
  Omega_bin <- ifelse(abs(Omega) > 0, 1, 0)
  diag(Omega_bin) <- 0
  rho_values <- seq(0, 200, by = 1)

  # Helper: ROC curve and AUC calculation
  roc_curve <- function(final_Omega_est_bin, Omega_bin, subtitle) {
    roc_df <- data.frame(FPR = numeric(), TPR = numeric(), rho = numeric())
    for (rho_idx in 1:length(rho_values)) {
      Omega_est_bin <- final_Omega_est_bin[, , rho_idx]
      pred <- prediction(as.vector(Omega_est_bin), as.vector(Omega_bin))
      perf <- performance(pred, "tpr", "fpr")
      roc_df <- rbind(roc_df,
                      data.frame(FPR = unlist(perf@x.values),
                                 TPR = unlist(perf@y.values),
                                 rho = rho_values[rho_idx]))
    }
    roc_df <- roc_df[order(roc_df$FPR), ]
    auc_value <- trapz(roc_df$FPR, roc_df$TPR)
    list(roc_df = roc_df, auc = auc_value)
  }

  # Helper: majority vote across datasets
  calculate_majority_vote <- function(datasets, rho_values) {
    num_datasets <- length(datasets)
    combined_Omega_est_bin <- array(0, dim = c(p, p, length(rho_values)))
    for (rho_idx in seq_along(rho_values)) {
      rho <- rho_values[rho_idx]
      Omega_est_list <- lapply(datasets, function(dat) {
        gl <- glasso(cov(dat), rho = rho)
        diag(gl$wi) <- 0
        ifelse(abs(gl$wi) > 0, 1, 0)
      })
      combined_Omega_est_bin[, , rho_idx] <- apply(
        array(unlist(Omega_est_list), c(p, p, num_datasets)),
        c(1, 2),
        function(x) ifelse(sum(x) > (num_datasets / 2), 1, 0)
      )
    }
    combined_Omega_est_bin
  }


  # (i) Proposed Method
  start_time <- Sys.time()
  dataset1 <- mget(paste0("pseudo_complete", 1:5))
  Omega_imputation_maj <- calculate_majority_vote(dataset1, rho_values)
  roc_result1 <- roc_curve(Omega_imputation_maj, Omega_bin, "Proposed Method")
  roc_dfs$p1 <- append(roc_dfs$p1, list(roc_result1$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "Proposed Method", Run = run, AUC = roc_result1$auc))
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "ROC by Proposed Method",
                                                time = difftime(end_time, start_time, units = "mins"))))


  # (ii) Complete Data (oracle)
  start_time <- Sys.time()
  dataset2 <- list(complete_data)
  Omega_oracle_bin <- calculate_majority_vote(dataset2, rho_values)
  roc_result2 <- roc_curve(Omega_oracle_bin, Omega_bin, "Complete Data")
  roc_dfs$p2 <- append(roc_dfs$p2, list(roc_result2$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "Complete Data", Run = run, AUC = roc_result2$auc))
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "ROC by Complete Data",
                                                time = difftime(end_time, start_time, units = "mins"))))


   # (iii) MICE default
  start_time <- Sys.time()
  imputed_data1 <- mice(data_na, m = 5, seed = 123, printFlag = FALSE)
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "MICE Default",
                                                time = difftime(end_time, start_time, units = "mins"))))

  start_time <- Sys.time()
  for (i in 1:5) assign(paste0("imputed_data1_", i), complete(imputed_data1, action = i))
  dataset3 <- mget(paste0("imputed_data1_", 1:5))
  Omega_mice_default_maj <- calculate_majority_vote(dataset3, rho_values)
  roc_result3 <- roc_curve(Omega_mice_default_maj, Omega_bin, "MICE Default")
  roc_dfs$p3 <- append(roc_dfs$p3, list(roc_result3$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "MICE Default", Run = run, AUC = roc_result3$auc))
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "ROC by MICE Default",
                                                time = difftime(end_time, start_time, units = "mins"))))


  # (iv) MICE with lasso.select.norm
  start_time <- Sys.time()
  imputed_data2 <- mice(data_na, method = "lasso.select.norm", m = 5, seed = 123, printFlag = FALSE)
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "MICE Lasso",
                                                time = difftime(end_time, start_time, units = "mins"))))

  start_time <- Sys.time()
  for (i in 1:5) assign(paste0("imputed_data2_", i), complete(imputed_data2, action = i))
  dataset4 <- mget(paste0("imputed_data2_", 1:5))
  Omega_mice_lasso_maj <- calculate_majority_vote(dataset4, rho_values)
  roc_result4 <- roc_curve(Omega_mice_lasso_maj, Omega_bin, "MICE Lasso")
  roc_dfs$p4 <- append(roc_dfs$p4, list(roc_result4$roc_df))
  auc_values <- rbind(auc_values, data.frame(Method = "MICE Lasso", Run = run, AUC = roc_result4$auc))
  end_time <- Sys.time()
  time_results <- append(time_results, list(list(run = run, step = "ROC by MICE Lasso",
                                                time = difftime(end_time, start_time, units = "mins"))))

  print(paste("Finished run", run))
}


# 3) Average ROC curves -----------------------------------
# (i) Average ROC across 10 runs for each method
average_roc_df <- function(roc_dfs) {
  fpr_vals <- seq(0, 1, length.out = 100)
  combined_roc_df <- data.frame(FPR = numeric(), TPR = numeric())
  for (roc_df in roc_dfs) {
    interp_tpr <- approx(roc_df$FPR, roc_df$TPR, xout = fpr_vals)$y
    combined_roc_df <- rbind(combined_roc_df, data.frame(FPR = fpr_vals, TPR = interp_tpr))
  }
  aggregate(TPR ~ FPR, data = combined_roc_df, mean)
}

# (ii) Plot helper with AUC in title
plot_roc_curve <- function(avg_roc_df, title) {
  auc_value <- trapz(avg_roc_df$FPR, avg_roc_df$TPR)
  ggplot(avg_roc_df, aes(x = FPR, y = TPR)) +
    geom_point(shape = 20) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = sprintf("ROC Curve (AUC = %.5f)", auc_value),
         subtitle = title, x = "1 - SP", y = "SE") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
}

average_roc_dfs <- lapply(roc_dfs, average_roc_df)
avg_plots <- list(
  plot_roc_curve(average_roc_dfs$p1, "Proposed Method"),
  plot_roc_curve(average_roc_dfs$p2, "Complete Data"),
  plot_roc_curve(average_roc_dfs$p3, "MICE Default"),
  plot_roc_curve(average_roc_dfs$p4, "MICE Lasso")
)


# 4) Final outputs ----------------------------------------
# (i) Show four average ROC plots
grid.arrange(grobs = avg_plots, ncol = 2)

# (ii) Print AUC table (sorted by Method)
auc_values <- auc_values[order(auc_values$Method), ]
print(auc_values)

# (iii) Print timing summary
time_df <- do.call(rbind, lapply(time_results, as.data.frame))
time_df$time <- format(time_df$time, scientific = FALSE)
print(time_df)
