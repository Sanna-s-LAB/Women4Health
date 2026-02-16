source("../utilities.R")

rel_nonlinear = as.character(commandArgs(trailingOnly = TRUE)[1])

# Create directory if it does not exist
if (!dir.exists(paste0("../Final_lagged/linear_on_non_linear_with_time/",rel_nonlinear))){
  dir.create(paste0("../Final_lagged/linear_on_non_linear_with_time/",rel_nonlinear))
}
setwd(paste0("../Final_lagged/linear_on_non_linear_with_time/",rel_nonlinear))


run_complex_sim <- function(n = 1000, T_points = 4, seed = 53, a = 0, sd = 0.1,
                            relation_types = c("quadratic", "cubic", "sin", "exp", "log","linear"),
                            ncores = detectCores() - 1) {

  single_sim <- function(l) {
    # Generate time series data
    data <- generate_multivariate_time_series_nl(
      n = n, T_points = T_points,
      relation_types = relation_types,
      trends= function(x) {0.2 * x},
      sd = sd
    )

    # Fit lagged linear models in long format
    fit_lm_YtoX <- fit_model_adaptive(data = data, X ~ Y_lag + X_lag, mode = "long")
    fit_lm_XtoY <- fit_model_adaptive(data = data, Y ~ X_lag + Y_lag, mode = "long")
    p_y2x <- extract_p_value(fit_lm_YtoX, "Y_lag")
    b_y2x <- fit_lm_YtoX$coefficients[[2]]
    p_x2y <- extract_p_value(fit_lm_XtoY, "X_lag")
    b_x2y <- fit_lm_XtoY$coefficients[[2]]

    # Convert data to wide format
    wide_data <- reshape(
      data,
      timevar = "t",
      idvar = "id",
      direction = "wide"
    )

    # Fit wide models
    fit_lm_YtoX_1 <- fit_model_adaptive(data = wide_data, X.2 ~ Y.1 + X.1, mode = "wide")
    fit_lm_XtoY_1 <- fit_model_adaptive(data = wide_data, Y.2 ~ X.1 + Y.1, mode = "wide")
    p_y2x_1 <- extract_p_value(fit_lm_YtoX_1, "Y.1")
    p_x2y_1 <- extract_p_value(fit_lm_XtoY_1, "X.1")

    fit_lm_YtoX_2 <- fit_model_adaptive(data = wide_data, X.3 ~ Y.2 + X.2, mode = "wide")
    fit_lm_XtoY_2 <- fit_model_adaptive(data = wide_data, Y.3 ~ X.2 + Y.2, mode = "wide")
    p_y2x_2 <- extract_p_value(fit_lm_YtoX_2, "Y.2")
    p_x2y_2 <- extract_p_value(fit_lm_XtoY_2, "X.2")

    fit_lm_YtoX_3 <- fit_model_adaptive(data = wide_data, X.4 ~ Y.3 + X.3, mode = "wide")
    fit_lm_XtoY_3 <- fit_model_adaptive(data = wide_data, Y.4 ~ X.3 + Y.3, mode = "wide")
    p_y2x_3 <- extract_p_value(fit_lm_YtoX_3, "Y.3")
    p_x2y_3 <- extract_p_value(fit_lm_XtoY_3, "X.3")

    list(
      p_values = data.frame(
        sim = l,
        p_X_to_Y = p_x2y,
        p_Y_to_X = p_y2x
      ),
      p_values_1 = data.frame(
        sim = l,
        p_X_to_Y = p_x2y_1,
        p_Y_to_X = p_y2x_1
      ),
      p_values_2 = data.frame(
        sim = l,
        p_X_to_Y = p_x2y_2,
        p_Y_to_X = p_y2x_2
      ),
      p_values_3 = data.frame(
        sim = l,
        p_X_to_Y = p_x2y_3,
        p_Y_to_X = p_y2x_3
      )
    )
  }

  # Run simulations in parallel
  results <- mclapply(1:200, single_sim, mc.cores = ncores)

  # Combine results
  p_values <- do.call(rbind, lapply(results, function(x) x$p_values))
  p_values_1 <- do.call(rbind, lapply(results, function(x) x$p_values_1))
  p_values_2 <- do.call(rbind, lapply(results, function(x) x$p_values_2))
  p_values_3 <- do.call(rbind, lapply(results, function(x) x$p_values_3))

  # Function to create plots for p-values
  create_plots <- function(p_values, model_name) {
    p_long <- melt(
      p_values,
      id.vars = c("sim"),
      measure.vars = c("p_X_to_Y", "p_Y_to_X"),
      variable.name = "stat", value.name = "value"
    )
    p_long$direction <- ifelse(grepl("X_to_Y", p_long$stat), "X → Y", "Y → X")
    p_long$metric <- "p_value"
    # True direction (X → Y)
    p_long$truth <- ifelse(p_long$direction == "X → Y", "true", "null")
    p_long$decision <- ifelse(
      p_long$metric == "p_value" & p_long$truth == "true" & p_long$value < 0.05, "TP",
      ifelse(p_long$metric == "p_value" & p_long$truth == "true" & p_long$value >= 0.05, "FN",
             ifelse(p_long$metric == "p_value" & p_long$truth == "null" & p_long$value < 0.05, "FP",
                    ifelse(p_long$metric == "p_value" & p_long$truth == "null" & p_long$value >= 0.05, "TN", NA))))

    plot_pvalue <- ggplot(subset(p_long, metric == "p_value"), aes(x = sim, y = -log10(value), color = direction)) +
      geom_point(size = 1.4, alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      facet_wrap(~direction, scales = "free_y") +
      theme_minimal() +
      labs(title = paste(model_name, "- GAM p-values"), x = "simulation index", y = "-log10(p-value)")

    list(p_long = p_long, plot_pvalue = plot_pvalue)
  }

  # Create plots for models 1–3
  list_1 <- create_plots(p_values_1, "Model 1")
  list_2 <- create_plots(p_values_2, "Model 2")
  list_3 <- create_plots(p_values_3, "Model 3")

  # Compute Power and Type I error for lag models
  p_long <- reshape2::melt(
    p_values,
    id.vars = c("sim"),
    measure.vars = c("p_X_to_Y", "p_Y_to_X"),
    variable.name = "stat", value.name = "value"
  )
  p_long$direction <- ifelse(grepl("X_to_Y", p_long$stat), "X → Y", "Y → X")
  p_long$metric <- "p_value"
  p_long$truth <- ifelse(p_long$direction == "X → Y", "true", "null")
  p_long$decision <- with(p_long, ifelse(
    metric == "p_value" & truth == "true" & value < 0.05, "TP",
    ifelse(metric == "p_value" & truth == "true" & value >= 0.05, "FN",
           ifelse(metric == "p_value" & truth == "null" & value < 0.05, "FP",
                  ifelse(metric == "p_value" & truth == "null" & value >= 0.05, "TN", NA)))
  ))

  levels_decision <- c("TP", "FP", "TN", "FN")

  # Create confusion matrix including all levels
  verdict <- table(factor(p_long$decision, levels = levels_decision))

  # Replace any NA with 0
  verdict[is.na(verdict)] <- 0

  TP <- verdict["TP"]
  FN <- verdict["FN"]
  FP <- verdict["FP"]
  TN <- verdict["TN"]

  # Compute Power and Type I error
  power <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  type1 <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)

  verdict_df <- data.frame(
    metric = c("Power", "Type I error"),
    value  = c(power, type1)
  )

  csv_to_save <- data.frame(TP = TP,
                            FN = FN,
                            FP = FP,
                            TN = TN,Power = power,T1E=type1)
  rownames(csv_to_save) <- NULL

  # -----------------
  # Plots
  # -----------------

  # p-values vs true effect
  p1 <- ggplot(subset(p_long, metric == "p_value"), aes(x = sim, y = -log10(value), color = direction)) +
    geom_point(size = 1.4, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    facet_wrap(~ direction, scales = "free_y") +
    labs(x = "simulation index", y = "-log10(p)", title = paste0("n=", n, ": p-value (lag models)")) +
    theme_minimal()

  # Power & Type I error plot
  p3 <- ggplot(verdict_df, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("Power" = "orange", "Type I error" = "steelblue")) +
    ylim(0, 1) +
    labs(
      title = paste0("n = ", n, ": Power and Type I Error Rate"),
      x = NULL,
      y = "Rate"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  return(list(
    relation_sim = relation_types,
    model_lag = list(p_values, p1, p3, csv_to_save),
    model1 = list(p_values = p_values_1, p_long = list_1$p_long, plot_pvalue = list_1$plot_pvalue),
    model2 = list(p_values = p_values_2, p_long = list_2$p_long, plot_pvalue = list_2$plot_pvalue),
    model3 = list(p_values = p_values_3, p_long = list_3$p_long, plot_pvalue = list_3$plot_pvalue)
  ))
}

# Sample sizes to simulate
n_values <- c(100, 300, 500)
results_list <- list()

# Run simulations for each sample size
for (nn in n_values) {

  cat("Simulating for n =", nn, "\n")
  result <- run_complex_sim(n = nn, rel =rel_nonlinear, seed=101)
  results_list[[as.character(nn)]] <- result

  plot_combined <- create_combined_plots_per_n(result, nn)[["plot"]]
  ggsave(paste0("Wide_models_n", nn, ".png"), plot_combined, width = 16, height = 10)
  write.csv(create_combined_plots_nl_per_n(result, nn)[["data"]],paste0("Wide_models_n", nn,"_", ".csv"))


  plots <- results_list[[as.character(nn)]][["model_lag"]]
  combined_plot <- (plots[[2]] |  plots[[3]]) +
    plot_annotation(title = paste("Simulation results for n =", as.character(nn)))
  ggsave(paste0("Long_model_n", nn, ".png"), combined_plot, width = 12, height = 8)
  write.csv(plots[[4]],paste0("Long_model_n", nn,"_", ".csv"))
}
