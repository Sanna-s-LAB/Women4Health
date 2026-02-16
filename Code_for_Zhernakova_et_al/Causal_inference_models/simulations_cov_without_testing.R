source("../utilities.R")

# a = as.numeric(commandArgs(trailingOnly = TRUE)[1])
a = 0.2

# Create output directory if it does not exist
if (!dir.exists("../Final_lagged/Covariates/Without_Z_simulation")) {
  dir.create("../Final_lagged/Covariates/Without_Z_simulation")
}
setwd("../Final_lagged/Covariates/Without_Z_simulation")

run_complex_sim <- function(n = 1000, T_points = 4, seed = 101, a = 0, ncores = 70) {

  b1_values <- seq(-2, 2, length.out = 200)
  se_values <- rnorm(200, 2, 0.5)

  single_sim <- function(l) {
    b1 <- b1_values[l]
    se <- se_values[l]

    f_list <- list(
      X ~ lag(X) + trend + error,
      Z ~ lag(Z) + error + trend,
      as.formula(paste("Y ~ lag(Y) +", b1, "* lag(X) + Z + error"))
    )

    trends <- list(
      X = function(t) a * t,
      Z = function(t) a * t
    )

    data <- generate_multivariate_time_series(
      n = n, T_points = T_points,
      formulas = f_list,
      trend_fun_list = trends, sd = se
    )

    fit_lm_YtoX <- fit_model_adaptive(data = data, X ~ Y_lag + X_lag, mode = "long")
    fit_lm_XtoY <- fit_model_adaptive(data = data, Y ~ X_lag + Y_lag, mode = "long")
    p_y2x <- extract_p_value(fit_lm_YtoX, "Y_lag")
    b_y2x <- fit_lm_YtoX$coefficients[[2]]
    p_x2y <- extract_p_value(fit_lm_XtoY, "X_lag")
    b_x2y <- fit_lm_XtoY$coefficients[[2]]

    # Base R wide reshape
    wide_data <- reshape(
      data,
      timevar = "t",
      idvar = "id",
      direction = "wide"
    )

    fit_lm_YtoX_1 <- fit_model_adaptive(data = wide_data, X.2 ~ Y.1 + X.1, mode = "wide")
    fit_lm_XtoY_1 <- fit_model_adaptive(data = wide_data, Y.2 ~ X.1 + Y.1, mode = "wide")
    p_y2x_1 <- extract_p_value(fit_lm_YtoX_1, "Y.1")
    b_y2x_1 <- fit_lm_YtoX_1$coefficients[[2]]
    p_x2y_1 <- extract_p_value(fit_lm_XtoY_1, "X.1")
    b_x2y_1 <- fit_lm_XtoY_1$coefficients[[2]]

    fit_lm_YtoX_2 <- fit_model_adaptive(data = wide_data, X.3 ~ Y.2 + X.2, mode = "wide")
    fit_lm_XtoY_2 <- fit_model_adaptive(data = wide_data, Y.3 ~ X.2 + Y.2, mode = "wide")
    p_y2x_2 <- extract_p_value(fit_lm_YtoX_2, "Y.2")
    b_y2x_2 <- fit_lm_YtoX_2$coefficients[[2]]
    p_x2y_2 <- extract_p_value(fit_lm_XtoY_2, "X.2")
    b_x2y_2 <- fit_lm_XtoY_2$coefficients[[2]]

    fit_lm_YtoX_3 <- fit_model_adaptive(data = wide_data, X.4 ~ Y.3 + X.3, mode = "wide")
    fit_lm_XtoY_3 <- fit_model_adaptive(data = wide_data, Y.4 ~ X.3 + Y.3, mode = "wide")
    p_y2x_3 <- extract_p_value(fit_lm_YtoX_3, "Y.3")
    b_y2x_3 <- fit_lm_YtoX_3$coefficients[[2]]
    p_x2y_3 <- extract_p_value(fit_lm_XtoY_3, "X.3")
    b_x2y_3 <- fit_lm_XtoY_3$coefficients[[2]]

    list(
      p_values = data.frame(
        b1 = b1, se = se,
        p_sem_X_to_Y = p_x2y, b_sem_X_to_Y = b_x2y,
        p_sem_Y_to_X = p_y2x, b_sem_Y_to_X = b_y2x
      ),
      p_values_1 = data.frame(
        sim = l, b1_true = b1, se = se,
        p_sem_X_to_Y = p_x2y_1, b_sem_X_to_Y = b_x2y_1,
        p_sem_Y_to_X = p_y2x_1, b_sem_Y_to_X = b_y2x_1
      ),
      p_values_2 = data.frame(
        sim = l, b1_true = b1, se = se,
        p_sem_X_to_Y = p_x2y_2, b_sem_X_to_Y = b_x2y_2,
        p_sem_Y_to_X = p_y2x_2, b_sem_Y_to_X = b_y2x_2
      ),
      p_values_3 = data.frame(
        sim = l, b1_true = b1, se = se,
        p_sem_X_to_Y = p_x2y_3, b_sem_X_to_Y = b_x2y_3,
        p_sem_Y_to_X = p_y2x_3, b_sem_Y_to_X = b_y2x_3
      )
    )
  }

  results <- mclapply(seq_along(b1_values), single_sim, mc.cores = ncores)

  p_values <- do.call(rbind, lapply(results, function(x) x$p_values))
  p_values_1 <- do.call(rbind, lapply(results, function(x) x$p_values_1))
  p_values_2 <- do.call(rbind, lapply(results, function(x) x$p_values_2))
  p_values_3 <- do.call(rbind, lapply(results, function(x) x$p_values_3))

  # Create plots
  create_plots <- function(p_values, model_name) {
    p_long <- melt(
      p_values,
      id.vars = c("sim", "b1_true", "se"),
      measure.vars = c("p_sem_X_to_Y", "b_sem_X_to_Y",
                       "p_sem_Y_to_X", "b_sem_Y_to_X"),
      variable.name = "stat", value.name = "value"
    )
    p_long$direction <- ifelse(grepl("X_to_Y", p_long$stat), "X → Y", "Y → X")
    p_long$metric <- ifelse(grepl("^p_", p_long$stat), "p_value", "Beta")
    p_long$truth <- ifelse(p_long$direction == "X → Y" & p_long$b1_true != 0, "true",
                           ifelse(p_long$metric == "p_value", "null", NA))
    p_long$decision <- ifelse(
      p_long$metric == "p_value" & p_long$truth == "true" & p_long$value < 0.05, "TP",
      ifelse(p_long$metric == "p_value" & p_long$truth == "true" & p_long$value >= 0.05, "FN",
             ifelse(p_long$metric == "p_value" & p_long$truth == "null" & p_long$value < 0.05, "FP",
                    ifelse(p_long$metric == "p_value" & p_long$truth == "null" & p_long$value >= 0.05, "TN", NA))))
    plot1 <- ggplot(subset(p_long, metric == "p_value"),
                    aes(x = b1_true, y = -log10(value), color = direction)) +
      geom_point(size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      facet_wrap(~direction) +
      theme_minimal() +
      labs(title = paste(model_name, "- SEM p-values"),
           x = "True b1", y = "-log10(p-value)")
    plot2 <- ggplot(subset(p_long, metric == "Beta" & direction == "X → Y"),
                    aes(x = b1_true, y = value)) +
      geom_point(color = "forestgreen", alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      theme_minimal() +
      labs(title = paste(model_name, "- Estimated vs True β"),
           x = "True b1", y = "Estimated β")
    list(plot_pvalue = plot1, plot_beta = plot2, p_long = p_long)
  }

  list_1 <- create_plots(p_values_1, "Model 1")
  list_2 <- create_plots(p_values_2, "Model 2")
  list_3 <- create_plots(p_values_3, "Model 3")

  p_long <- reshape2::melt(
    p_values,
    id.vars = c("b1", "se"),
    measure.vars = c("p_sem_X_to_Y", "b_sem_X_to_Y",
                     "p_sem_Y_to_X", "b_sem_Y_to_X"),
    variable.name = "stat", value.name = "value"
  )
  p_long$direction <- ifelse(grepl("X_to_Y", p_long$stat), "X → Y", "Y → X")
  p_long$metric <- ifelse(grepl("^p_", p_long$stat), "p_value", "Beta")
  p_long$truth <- ifelse(p_long$direction == "X → Y" & p_long$b1 != 0, "true",
                         ifelse(p_long$metric == "p_value", "null", NA))
  p_long$decision <- with(p_long, ifelse(
    metric == "p_value" & truth == "true" & value < 0.05, "TP",
    ifelse(metric == "p_value" & truth == "true" & value >= 0.05, "FN",
           ifelse(metric == "p_value" & truth == "null" & value < 0.05, "FP",
                  ifelse(metric == "p_value" & truth == "null" & value >= 0.05, "TN", NA)))))
  verdict <- table(p_long$decision, useNA = "no")

  TP <- verdict[["TP"]]
  FN <- verdict[["FN"]]
  FP <- verdict[["FP"]]
  TN <- verdict[["TN"]]

  power <- TP / (TP + FN)
  type1 <- FP / (FP + TN)

  verdict_df <- data.frame(
    metric = c("Power", "Type I error"),
    value  = c(power, type1)
  )
  csv_to_save <- data.frame(TP = TP,
                            FN = FN,
                            FP = FP,
                            TN = TN,
                            Power = power,
                            T1E = type1)

  # Plotting
  p1 <- ggplot(subset(p_long, metric == "p_value"),
               aes(x = b1, y = -log10(value), color = direction)) +
    geom_point(size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    facet_wrap(~ direction) +
    labs(x = "True effect (b1)", y = "-log10(p)",
         title = paste0("n=", n, ": p-value vs true effect")) +
    theme_minimal()
  p2 <- ggplot(subset(p_long, metric == "Beta" & direction == "X → Y"),
               aes(x = b1, y = value)) +
    geom_point(alpha = 0.6, color = "forestgreen") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "True effect", y = "Estimated effect",
         title = paste0("n=", n, ": Estimated vs true")) +
    theme_minimal()
  p3 <- ggplot(verdict_df, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c(
      "Power" = "orange",
      "Type I error" = "steelblue"
    )) +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.position = "none")

  return(list(
    model_lag = list(p_values, p1, p2, p3, data = csv_to_save),
    model1 = list(p_values = p_values_1, p_long = list_1$p_long,
                  plot_pvalue = list_1$plot_pvalue,
                  plot_beta = list_1$plot_beta),
    model2 = list(p_values = p_values_2, p_long = list_2$p_long,
                  plot_pvalue = list_2$plot_pvalue,
                  plot_beta = list_2$plot_beta),
    model3 = list(p_values = p_values_3, p_long = list_3$p_long,
                  plot_pvalue = list_3$plot_pvalue,
                  plot_beta = list_3$plot_beta)
  ))
}

n_values <- c(100, 300, 500)
results_list <- list()

for (nn in n_values) {
  cat("Simulating for n =", nn, "\n")

  # Run simulation
  result <- run_complex_sim(n = nn, a = a)
  results_list[[as.character(nn)]] <- result

  # Save wide-format combined plots and CSV
  plot_combined <- create_combined_plots_per_n(result, nn)[["plot"]]
  ggsave(paste0("Wide_models_n", nn, ".png"),
         plot_combined, width = 16, height = 10)
  write.csv(create_combined_plots_per_n(result, nn)[["data"]],
            paste0("Wide_models_n", nn, ".csv"))

  plots <- results_list[[as.character(nn)]][["model_lag"]]

  # Combine long-format plots
  combined_plot <- (plots[[2]] | plots[[3]]) / plots[[4]] +
    plot_annotation(title = paste("Simulation results for n =", as.character(nn)))

  # Save long-format results
  ggsave(paste0("Long_model_n", nn, ".png"), combined_plot, width = 12, height = 8)
  write.csv(plots[[5]], paste0("Long_model_n", nn, ".csv"))
}
