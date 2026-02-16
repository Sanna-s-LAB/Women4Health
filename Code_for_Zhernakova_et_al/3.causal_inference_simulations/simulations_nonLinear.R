source("../utilities.R")


rel_nonlinear = as.character(commandArgs(trailingOnly = TRUE)[1])
out_dir <- "../Final_lagged/Nonlinear/Nonlinear"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}
setwd(out_dir)

# --- Run full simulation with GAM fits in parallel ---
run_complex_sim <- function(n = 1000, T_points = 4, seed = NULL, rel = c("quadratic", "cubic", "sin", "exp", "log", "linear"), ncores = 30) {
  if(!is.null(seed)) set.seed(seed)

  # Choose one relation type for the simulation
  rel <- sample(rel,1)

  # Function for a single simulation
  single_sim <- function(l) {
    tryCatch({
      data_sim <- generate_multivariate_time_series_nl(n = n, T_points = T_points, relation_types = rel)

      # --- LONG GAM fits ---
      fit_lm_XtoY <- fit_gam_adaptive(
        data_sim,
        Y ~ s(X_lag, k=3, bs="tp", sp=1.5) +
          s(Y_lag, k=3, bs="tp", sp=1.5) +
          s(id, bs="re"),
        mode = "long"
      )

      fit_lm_YtoX <- fit_gam_adaptive(
        data_sim,
        X ~ s(Y_lag, k=3, bs="tp", sp=1.5) +
          s(X_lag, k=3, bs="tp", sp=1.5) +
          s(id, bs="re"),
        mode = "long"
      )

      p_y2x <- extract_p_value_gam(fit_lm_YtoX)
      p_x2y <- extract_p_value_gam(fit_lm_XtoY)

      # --- WIDE reshape ---
      wide_data <- reshape(
        data_sim[,c("t","id","X","Y")],
        timevar = "t",
        idvar = "id",
        direction = "wide"
      )

      # WIDE GAM fits for all lags (no random effect)
      fit_lm_XtoY_1 <- fit_gam_adaptive(wide_data, Y.2 ~ s(X.1, k=5) +  s(Y.1, k=5), mode = "wide")
      fit_lm_YtoX_1 <- fit_gam_adaptive(wide_data, X.2 ~ s(Y.1, k=5) +  s(X.1, k=5), mode = "wide")
      p_y2x_1 <- extract_p_value_gam(fit_lm_YtoX_1)
      p_x2y_1 <- extract_p_value_gam(fit_lm_XtoY_1)

      fit_lm_XtoY_2 <- fit_gam_adaptive(wide_data, Y.3 ~ s(X.2, k=5) +  s(Y.2, k=5), mode = "wide")
      fit_lm_YtoX_2 <- fit_gam_adaptive(wide_data, X.3 ~ s(Y.2, k=5) +  s(X.2, k=5), mode = "wide")
      p_y2x_2 <- extract_p_value_gam(fit_lm_YtoX_2)
      p_x2y_2 <- extract_p_value_gam(fit_lm_XtoY_2)

      fit_lm_XtoY_3 <- fit_gam_adaptive(wide_data, Y.4 ~ s(X.3, k=5) + s(Y.3, k=5), mode = "wide")
      fit_lm_YtoX_3 <- fit_gam_adaptive(wide_data, X.4 ~ s(Y.3, k=5) +  s(X.3, k=5), mode = "wide")
      p_y2x_3 <- extract_p_value_gam(fit_lm_YtoX_3)
      p_x2y_3 <- extract_p_value_gam(fit_lm_XtoY_3)

      # Return list of p-values
      list(
        p_values = data.frame(sim = l, p_X_to_Y = p_x2y, p_Y_to_X = p_y2x),
        p_values_1 = data.frame(sim = l, p_X_to_Y = p_x2y_1, p_Y_to_X = p_y2x_1),
        p_values_2 = data.frame(sim = l, p_X_to_Y = p_x2y_2, p_Y_to_X = p_y2x_2),
        p_values_3 = data.frame(sim = l, p_X_to_Y = p_x2y_3, p_Y_to_X = p_y2x_3)
      )
    }, error = function(e) {
      # Return NAs if a single simulation fails
      na_row <- data.frame(sim = l, p_X_to_Y = NA_real_, p_Y_to_X = NA_real_)
      list(p_values = na_row, p_values_1 = na_row, p_values_2 = na_row, p_values_3 = na_row)
    })
  }

  # Run 200 simulations in parallel
  results <- parallel::mclapply(1:200, single_sim, mc.cores = ncores)

  # Combine p-values
  p_values   <- do.call(rbind, lapply(results, function(x) x$p_values))
  p_values_1 <- do.call(rbind, lapply(results, function(x) x$p_values_1))
  p_values_2 <- do.call(rbind, lapply(results, function(x) x$p_values_2))
  p_values_3 <- do.call(rbind, lapply(results, function(x) x$p_values_3))

  # --- Create plots for WIDE models ---
  create_plots <- function(p_values, model_name) {
    p_long <- melt(
      p_values,
      id.vars = c("sim"),
      measure.vars = c("p_X_to_Y", "p_Y_to_X"),
      variable.name = "stat", value.name = "value"
    )
    p_long$direction <- ifelse(grepl("X_to_Y", p_long$stat), "X → Y", "Y → X")
    p_long$metric <- "p_value"
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

  list_1 <- create_plots(p_values_1, "Model 1")
  list_2 <- create_plots(p_values_2, "Model 2")
  list_3 <- create_plots(p_values_3, "Model 3")

  # --- LONG model lag logic ---
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
                  ifelse(metric == "p_value" & truth == "null" & value >= 0.05, "TN", NA)))))

  # Confusion matrix counts
  verdict <- table(factor(p_long$decision, levels = c("TP", "FN", "FP", "TN")), useNA = "no")
  TP <- verdict[["TP"]]
  FN <- verdict[["FN"]]
  FP <- verdict[["FP"]]
  TN <- verdict[["TN"]]

  # Performance metrics
  power  <- TP / (TP + FN)
  type1  <- FP / (FP + TN)

  verdict_df <- data.frame(metric = c("Power", "Type I error"), value = c(power, type1))
  csv_to_save <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN, Power = power, T1E = type1)

  # --- PLOTS ---
  p1 <- ggplot(subset(p_long, metric == "p_value"), aes(x = sim, y = -log10(value), color = direction)) +
    geom_point(size = 1.4, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    facet_wrap(~ direction, scales = "free_y") +
    labs(x = "simulation index", y = "-log10(p)", title = paste0("n=", n, ": p-value (lag models)")) +
    theme_minimal()

  p3 <- ggplot(verdict_df, aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("Power" = "orange", "Type I error" = "steelblue")) +
    ylim(0, 1) +
    labs(title = paste0("n = ", n, ": Power and Type I Error Rate"), x = NULL, y = "Rate") +
    theme_minimal() +
    theme(legend.position = "none")

  return(list(
    relation_sim = rel,
    model_lag = list(p_values, p1, p3, csv_to_save),
    model1 = list(p_values = p_values_1, p_long = list_1$p_long, plot_pvalue = list_1$plot_pvalue),
    model2 = list(p_values = p_values_2, p_long = list_2$p_long, plot_pvalue = list_2$plot_pvalue),
    model3 = list(p_values = p_values_3, p_long = list_3$p_long, plot_pvalue = list_3$plot_pvalue)
  ))
}

                                      
# --- Batch simulations across multiple sample sizes ---
n_values <- c(100, 300, 500)
results_list <- list()

for (nn in n_values) {
  cat("Simulating for n =", nn, "\n")
  result <- run_complex_sim(n = nn, rel = rel_nonlinear, seed = 101)
  results_list[[as.character(nn)]] <- result

  plot_combined <- create_combined_plots_nl_per_n(result, nn)[["plot"]]
  ggsave(paste0("Wide_models_n", nn, "_", result$relation_sim, ".png"), plot_combined, width = 16, height = 10)
  write.csv(create_combined_plots_per_n(result, nn)[["data"]], paste0("Wide_models_n", nn, "_", result$relation_sim, ".csv"))

  plots <- results_list[[as.character(nn)]][["model_lag"]]
  combined_plot <- (plots[[2]] | plots[[3]]) + plot_annotation(title = paste("Simulation results for n =", nn))
  ggsave(paste0("Long_model_n", nn, "_", result$relation_sim, ".png"), combined_plot, width = 12, height = 8)
  write.csv(plots[[4]], paste0("Long_model_n", nn, "_", result$relation_sim, ".csv"))
}
