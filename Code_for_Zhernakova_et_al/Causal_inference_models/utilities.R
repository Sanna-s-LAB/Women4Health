library(broom)
library(tidyr)
library(purrr)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(reshape2)
library(rlang)
library(mgcv)
library(parallel)
# Function to extract p-value from a linear model for a given term
extract_p_value <- function(model, term) {
  if (term %in% rownames(coef(summary(model)))) {
    return(coef(summary(model))[term, "Pr(>|t|)"])
  } else {
    return(NA)
  }
}

# --- Utility function: extract p-value from the first smooth term of a GAM if it exists ---
extract_p_value_gam <- function(model) {
  st <- summary(model)$s.table
  if (is.null(st) || nrow(st) == 0) return(NA_real_)  # Avoid errors if s.table is missing
  as.numeric(st[1, 4])
}


# Function to fit a linear model in either long or wide format
fit_model_adaptive <- function(data, formula,
                               id_var = "id", time_var = "t",
                               mode = c("long", "wide")) {
  mode <- match.arg(mode)

  if (mode == "long") {
    # Create lag variables if they do not exist
    vars_in_formula <- all.vars(formula)
    predictors <- vars_in_formula[-1]
    for (pred in predictors) {
      if (grepl("_lag$", pred)) {
        base <- sub("_lag$", "", pred)
        if (!pred %in% names(data)) {
          data[[pred]] <- ave(data[[base]], data[[id_var]],
                              FUN = function(x) c(NA, head(x, -1)))
        }
      }
    }
  } else if (mode == "wide") {
    # Convert long → wide if needed
    vars_in_formula <- all.vars(formula)
    needs_wide <- any(grepl("_[0-9]+$", vars_in_formula))
    if (needs_wide) {
      long_vars <- unique(gsub("_[0-9]+$", "", vars_in_formula))
      # Pivot only if at least one wide column is missing
      missing <- !vars_in_formula %in% names(data)
      if (any(missing)) {
        data <- data %>%
          pivot_wider(names_from = !!sym(time_var),
                      values_from = all_of(long_vars))
      }
    }
  }

  fit <- lm(formula, data = data)
  return(fit)
}

# Generate multivariate time series data with optional missingness
generate_multivariate_time_series <- function(
    n = 1000,
    T_points = 4,
    formulas = list(),
    trend_fun_list = list(),
    sd = 0.1,
    seed = NULL,
    missing_value = FALSE, perc = 0.1
) {

  if(!is.null(seed)){
    set.seed(seed)}
  var_names <- sapply(formulas, function(f) as.character(f[[2]]))
  id <- rep(1:n, each = T_points)
  t <- rep(1:T_points, times = n)
  data <- data.frame(id = id, t = t)
  for (v in var_names) {
    data[[v]] <- NA
  }
  for (v in var_names) {
    init <- rnorm(n, 0, 1)
    data[data$t == 1, v] <- init
  }
  for (j in 1:n) {
    idx <- which(data$id == j)
    for (i in 2:T_points) {
      row <- idx[i]
      for (f in formulas) {
        lhs <- as.character(f[[2]])
        rhs <- f[[3]]
        env <- new_environment(parent = baseenv())
        env$trend <- if (lhs %in% names(trend_fun_list)) {
          trend_fun_list[[lhs]](i)
        } else 0
        env$error <- rnorm(1, 0, sd)
        for (v in var_names) {
          env[[v]] <- data[row, v]
        }
        # ⚙️ lag() based on ave()
        env$lag <- function(var) {
          var_name <- deparse(substitute(var))
          lags <- ave(data[[var_name]], data$id, FUN = function(x) c(NA, head(x, -1)))
          return(lags[row])
        }
        data[row, lhs] <- eval(rhs, env)
      }
    }
  }

  if(missing_value == 'MNAR'){
    n_miss <- round(perc * n * T_points)  # total proportion missing

    # Distribution of missing values across time points
    n_miss_t3 <- round(n_miss * 0.30)
    n_miss_t4 <- round(n_miss * 0.40)
    n_miss_12 <- n_miss - n_miss_t4 - n_miss_t3  # remainder for t1–t2

    # Indices for each time point
    idx_t1 <- which(data$t == 1)
    idx_t2 <- which(data$t == 2)
    idx_t3 <- which(data$t == 3)
    idx_t4 <- which(data$t == 4)

    # Sample indices to be set to NA (common for X and Y)
    missing_idx <- c(
      sample(c(idx_t1, idx_t2), n_miss_12, replace = FALSE),
      sample(idx_t3, n_miss_t3, replace = FALSE),
      sample(idx_t4, n_miss_t4, replace = FALSE)
    )

    # Apply NA to all variables except id and t
    for(col in setdiff(colnames(data),c('id','t'))){
      data[missing_idx,col] <- NA
    }

  }else if(missing_value == "MAR"){
    missing_idx <- sample(1:(n * T), size = round(perc * n * T))
    for(col in setdiff(colnames(data),c('id','t'))){
      data[missing_idx,col] <- NA
    }
  } else {
    data <- data
  }
  return(data)
}

# --- Fit GAM adaptively for long or wide data ---
fit_gam_adaptive <- function(data, formula, id_var = "id", time_var = "t", mode = c("long","wide")) {
  mode <- match.arg(mode)

  # Long mode: create lagged variables if predictor ends with "_lag"
  if(mode == "long") {
    vars_in_formula <- all.vars(formula)
    predictors <- vars_in_formula[-1]
    for(pred in predictors) {
      if(grepl("_lag$", pred)) {
        base <- sub("_lag$", "", pred)
        if(!pred %in% names(data)) {
          data[[pred]] <- ave(data[[base]], data[[id_var]], FUN = function(x) c(NA, head(x, -1)))
        }
      }
    }
  }

  # Wide mode: pivot if necessary (here base R reshape is used)
  if(mode == "wide") {
    vars_in_formula <- all.vars(formula)
    needs_wide <- any(grepl("_[0-9]+$", vars_in_formula))
    if(needs_wide) {
      long_vars <- unique(gsub("_[0-9]+$", "", vars_in_formula))
      missing <- !vars_in_formula %in% names(data)
      if(any(missing)) {
        data <- data %>%
          pivot_wider(names_from = !!sym(time_var), values_from = all_of(long_vars))
      }
    }
  }

  # Remove missing values and fit GAM
  data <- na.omit(data)
  fit <- gam(formula, data = data, method = "REML")
  return(fit)
}

# Function to combine plots for each sample size
create_combined_plots_per_n <- function(result_list, n) {
  plots_pval <- list()
  plots_perf <- list()
  csv_to_save <- data.frame(TP = numeric(3), FN = numeric(3), FP = numeric(3),
                            TN = numeric(3), Power = numeric(3), T1E = numeric(3))

  for (i in 1:3) {
    res <- result_list[[paste0("model", i)]]
    p_long <- res$p_long
    p_long$b1 <- p_long$b1_true

    # Confusion matrix
    verdict <- p_long %>%
      dplyr::filter(metric == "p_value", !is.na(decision)) %>%
      dplyr::count(decision) %>%
      tibble::deframe()

    TP <- verdict[["TP"]]; FN <- verdict[["FN"]]
    FP <- verdict[["FP"]]; TN <- verdict[["TN"]]
    power <- TP / (TP + FN); type1 <- FP / (FP + TN)

    perf_df <- data.frame(metric = c("Power", "Type I error"), value = c(power, type1))
    csv_to_save[i, ] <- c(TP, FN, FP, TN, power, type1)

    # Plot p-values
    pval_plot <- ggplot(subset(p_long, metric == "p_value"), aes(x = b1, y = -log10(value), color = direction)) +
      geom_point(size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      facet_wrap(~ direction) +
      labs(x = "True effect (b1)", y = "-log10(p)", title = paste0("Model ", i, " (n = ", n, ")")) +
      theme_minimal()

    # Plot performance
    perf_plot <- ggplot(perf_df, aes(x = metric, y = value, fill = metric)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_manual(values = c("Power" = "orange", "Type I error" = "steelblue")) +
      ylim(0, 1) +
      labs(y = "Rate", title = paste0("Model ", i, " – Performance")) +
      theme_minimal() +
      theme(legend.position = "none")

    plots_pval[[i]] <- pval_plot
    plots_perf[[i]] <- perf_plot
  }

  combined_pval <- plots_pval[[1]] + plots_pval[[2]] + plots_pval[[3]] + plot_layout(nrow = 1)
  combined_perf <- plots_perf[[1]] + plots_perf[[2]] + plots_perf[[3]] + plot_layout(nrow = 1)
  final_plot <- (combined_pval / combined_perf) + plot_annotation(title = paste("Simulation Results for n =", n))

  return(list(plot = final_plot, data = csv_to_save))
}


generate_multivariate_time_series_nl <- function(
    n = 100,
    T_points = 4,
    seed = NULL,
    relation_types = c("quadratic", "cubic", "sin", "exp", "log","linear"),
    trends = function(x) {0.2 * x}, sd = 0.1
) {
  if(!is.null(seed)) set.seed(seed)

  # Choose a random relationship type for the whole dataset
  rel_type <- sample(relation_types, 1)

  id <- rep(1:n, each = T_points)
  t <- rep(1:T_points, times = n)
  data <- data.frame(id = id, t = t, X = NA_real_, Y = NA_real_, relation_type = rel_type)

  # Initialize first timepoint
  data$X[data$t == 1] <- rnorm(n)
  data$Y[data$t == 1] <- rnorm(n)

  for (i in 1:n) {
    idx <- which(data$id == i)
    for (j in 2:T_points) {
      row <- idx[j]
      prev <- idx[j-1]

      # X autocorrelated
      data$X[row] <- trends(j) + data$X[prev] + rnorm(1, 0, sd)

      # Y depends on lagged X with chosen relationship
      data$Y[row] <- switch(
        rel_type,
        "quadratic" =  data$X[prev]^2 + data$Y[prev] + rnorm(1, 0, sd),
        "cubic"     =  data$X[prev]^3 + data$Y[prev] + rnorm(1, 0, sd),
        "sin"       =  sin(data$X[prev]) + data$Y[prev] + rnorm(1, 0, sd),
        "exp"       =  exp(data$X[prev]/3) + data$Y[prev] + rnorm(1, 0, sd),
        "log"       =  log(abs(data$X[prev]) + 1) + data$Y[prev] + rnorm(1, 0, sd),
        "linear"    =  2* data$X[prev]  + data$Y[prev] + rnorm(1, 0, sd),
        # default
        sin(data$X[prev]) + data$Y[prev] + rnorm(1, 0, sd)
      )
    }
  }

  data$id <- as.factor(data$id)  # Useful for s(id, bs="re") in long format
  return(data)
}



create_combined_plots_nl_per_n <- function(result_list, n) {

  plots_pval <- list()
  plots_perf <- list()
  csv_to_save <- data.frame(
    TP = numeric(3),
    FN = numeric(3),
    FP = numeric(3),
    TN = numeric(3),
    Power = numeric(3),
    T1E = numeric(3)
  )
  for (i in 1:3) {

    res <- result_list[[paste0("model", i)]]
    p_long <- res$p_long

    # Ensure simulation index is correct
    p_long$b1 <-p_long$sim

    # Compute confusion matrix (only p-values)
    verdict <- table(
      factor(p_long$decision,
             levels = c("TP","FN","FP","TN")),
      useNA = "no"
    )

    TP <- verdict[["TP"]]
    FN <- verdict[["FN"]]
    FP <- verdict[["FP"]]
    TN <- verdict[["TN"]]

    power <- TP / (TP + FN)
    type1 <- FP / (FP + TN)

    perf_df <- data.frame(
      metric = c("Power", "Type I error"),
      value  = c(power, type1)
    )

    csv_to_save[i, ] <- c(TP, FN, FP, TN, power, type1)

    # -----------------
    # Plot p-values
    # -----------------
    pval_plot <- ggplot(
      subset(p_long, metric == "p_value"),
      aes(x = b1, y = -log10(value), color = direction)
    ) +
      geom_point(size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      facet_wrap(~ direction) +
      labs(
        x = "simulation index",
        y = "-log10(p)",
        title = paste0("Model ", i, " (n = ", n, ")")
      ) +
      theme_minimal()

    # -----------------
    # Plot Power & Type I
    # -----------------
    perf_plot <- ggplot(perf_df, aes(x = metric, y = value, fill = metric)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_manual(values = c(
        "Power" = "orange",
        "Type I error" = "steelblue"
      )) +
      ylim(0, 1) +
      labs(
        y = "Rate",
        title = paste0("Model ", i, " – Performance")
      ) +
      theme_minimal() +
      theme(legend.position = "none")

    plots_pval[[i]]  <- pval_plot
    plots_perf[[i]]  <- perf_plot
  }

  # Combinations
  combined_pval <- plots_pval[[1]] + plots_pval[[2]] + plots_pval[[3]] +
    plot_layout(nrow = 1)

  combined_perf <- plots_perf[[1]] + plots_perf[[2]] + plots_perf[[3]] +
    plot_layout(nrow = 1)

  # Final layout
  final_plot <- (combined_pval / combined_perf) +
    plot_annotation(title = paste("Simulation Results for n =", n))

  return(list(plot=final_plot, data= csv_to_save))
}

