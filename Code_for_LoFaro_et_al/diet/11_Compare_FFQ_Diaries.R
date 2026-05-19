# PAOLA FORABOSCO 28/04/2026

# ============================================================================ #
#  DESCRIPTION
# ============================================================================ #

# This script performs two complementary analyses on nutritional traits
# measured via FFQ (Food Frequency Questionnaire) and Daily Diaries :
#
# 1. Paired comparisons of FFQ vs Diaries:
#    - Paired Wilcoxon Signed-Rank tests are applied
#    - Boxplots are generated with significance annotations
#    - Raw p-values are stored and Bonferroni correction is applied
#    - Combined boxplots are exported as a JPEG image
#
# 2. Spearman correlation between FFQ and Diaries measurements:
#    - Spearman correlation coefficient (Rs) and p-value (Pr) are computed
#    - Scatterplots with regression line and confidence interval are generated
#    - Marginal histograms are added to scatterplots
#    - Results are stored for export and Bonferroni correction is applied
#    - Combined correlation plots are exported as a JPEG image
#
# Input:
#   - Wide-format CSV with FFQ and Diaries (Total Means) measurements for each trait
# Output:
#   - CSV files with Wilcoxon and Spearman p-values (raw and Bonferroni-adjusted)
#   - JPEG files with combined boxplots and correlation plots


# ============================================================================ #
#  INITIALIZATION
# ============================================================================ #
# Load required libraries for plotting, data manipulation, and statistics
library(ggplot2)    # plotting
library(dplyr)      # data manipulation
library(ggsignif)   # significance bars for boxplots
library(gridExtra)  # arrange multiple plots in a grid
library(grid)       # textGrob and gpar for titles
library(psych)      # corr.test() for Spearman p-values
library(ggExtra)    # marginal histograms on scatterplots
library(tidyr)      # for drop_na()

# Read input data (wide format with FFQ and Diaries)
dat <- read.csv("FFQ_DiariesTM_wide.csv", header = TRUE)

# Label used for output file names
label = "FFQ_Diaries"               

# List of nutritional traits to analyze
traits <- c("Calories", "Carbohydrates", "Fats", "Proteins",
            "Cholesterol", "Sodium", "Sugar", "Fiber")

# Initialize dataframe to store all pairwise test results
all_pvalues <- data.frame()

# Initialize list to store boxplots
plots <- list()          

# ============================================================================ #
#  BOXPLOT GENERATION & PAIRWISE TESTS
# ============================================================================ #

# Loop through each nutritional trait
for (trait in traits) {
  
  # Subset FFQ and Diaries columns corresponding to the current trait
  dat.trait <- data.frame(dat[ , grepl(trait, names(dat))])
  colnames(dat.trait ) <- c("FFQ", "Diaries")
  
  # Add subject ID
  dat.trait$ID <- dat$ID
  
  # Create long-format dataframe with FFQ vs Diaries as categories
  data <- data.frame(
    category = factor(rep(c("FFQ", "Diaries"), each = nrow(dat))),
    value    = c(dat.trait$FFQ, dat.trait$Diaries),
    label    = dat.trait$ID
  )
  
  # Count number of valid (non-missing) observations per category
  n_counts <- data %>%
    drop_na(value) %>%           # Remove missing values
    group_by(category) %>%       # Group by FFQ / Diaries
    summarise(N = n()) %>%       # Count observations
    ungroup()
  
  # Define category labels (optionally including sample size)
  # category_labels <- paste0(n_counts$category, " (N = ", n_counts$N, ")")
  category_labels <- n_counts$category
  
  # Update factor levels to ensure correct ordering and labeling
  data$category <- factor(
    data$category,
    levels = n_counts$category,
    labels = category_labels
  )
  
  # ------- #
  # Boxplot #
  # ------- #
  
  # Generate boxplot for the current trait
  bp <- ggplot(data, aes(x = category, y = value)) +
    geom_boxplot(outlier.shape = 16, fill = "royalblue1", alpha = 0.7) +
    stat_boxplot(geom = "errorbar", width = 0.1) +
    labs(
      title = "",
      x = NULL,
      y = trait
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, color = "royalblue4", face = "bold")
    )
  
  # ----------------------------------- #
  # Pairwise Wilcoxon Test (PAIRED 2×2) #
  # ----------------------------------- #
  
  # Extract category levels (FFQ and Diaries)
  groups <- levels(data$category)
  
  # Generate all pairwise comparisons
  comparisons <- combn(groups, 2, simplify = FALSE)
  n_comparisons <- length(comparisons)
  
  # Function to perform paired Wilcoxon Signed-Rank test
  paired_wilcoxon_pairwise <- function(comparison) {
    group1 <- comparison[1]
    group2 <- comparison[2]
    
    # Extract data for each group
    data1 <- data[data$category == group1, c("label", "value")]
    data2 <- data[data$category == group2, c("label", "value")]
    
    # Identify subjects with non-missing paired values
    common_ids <- intersect(
      data1$label[!is.na(data1$value)], 
      data2$label[!is.na(data2$value)]
    )
    
    # Align paired values by subject ID
    paired_data1 <- data1$value[match(common_ids, data1$label)]
    paired_data2 <- data2$value[match(common_ids, data2$label)]
    
    # Perform paired Wilcoxon Signed-Rank test
    test_result <- wilcox.test(paired_data1, paired_data2,
                               paired = TRUE, exact = FALSE)
    
    # Return p-value and number of paired observations
    return(list(p_value = test_result$p.value,
                n_paired = length(common_ids)))
  }
  
  # Apply paired Wilcoxon test to all comparisons
  test_results <- lapply(comparisons, paired_wilcoxon_pairwise)
  
  # Extract p-values and paired sample sizes
  p_values  <- sapply(test_results, function(x) x$p_value)
  pairwise_n <- sapply(test_results, function(x) x$n_paired)
  
  # -------------------------------- #
  # Add significance bars to boxplot #
  # -------------------------------- #
  
  # Initialize annotated boxplot
  bpp <- bp
  
  # Define vertical positions for significance bars
  max_value <- max(data$value, na.rm = TRUE)
  y_positions <- max_value + (1:n_comparisons) * max_value / 10
  
  # Add significance annotations for each comparison
  for (i in seq_along(comparisons)) {
    bpp <- bpp +
      geom_signif(
        comparisons = list(comparisons[[i]]),
        annotations = ifelse(p_values[i] < 0.001, "***",
                             ifelse(p_values[i] < 0.01, "**",
                                    ifelse(p_values[i] < 0.05, "*", "ns"))),
        y_position = y_positions[i],
        tip_length = 0.01
      )
  }
  
  # ------------------------------ #
  # Create a results summary table #
  # ------------------------------ #
  
  results_table <- data.frame(
    trait,
    Group1 = sapply(comparisons, function(x) x[1]),
    Group2 = sapply(comparisons, function(x) x[2]),
    Paired_N = pairwise_n,
    pvalue = p_values
  )
  
  # Append results for the current trait
  all_pvalues <- rbind(all_pvalues, results_table)
  
  # Store plot
  plots[[trait]] <- bpp
}

# ============================================================================ #
#  EXPORT ALL P-VALUES AND ADD BONFERRONI CORRECTION
# ============================================================================ #

# Apply Bonferroni correction across all pairwise tests
all_pvalues$p.adj.Bonf <- p.adjust(all_pvalues$pvalue, method = "bonferroni")

# Export results table
write.csv(all_pvalues, paste0(label, "_Wilcoxcon_pvalues.csv"),
          row.names = FALSE)

# NOTE:
# height = 1400 for 8 plots
# height = 800 for 4 plots

# Save all boxplots into a single JPEG file
jpeg(file = paste0(label, "_means.jpeg"),
     width = 1000, height = 1400, units = "px", res = 100)

grid.arrange(
  grobs = plots,
  ncol = 2,
  top = textGrob(
    paste0("Paired Wilcoxon Signed-Rank test", "\n"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)

dev.off()


# ============================================================================ #
#  SPEARMAN CORRELATIONS 
# ============================================================================ #

# Initialize dataframe to store all p-values
all_pvalues <- data.frame()

# Initialize list to store plots
plots <- list()          

# Loop through each nutritional trait
for (trait in traits) {
  
  # Subset columns corresponding to the current trait
  dat.trait <- data.frame(dat[, grepl(trait, names(dat))])
  colnames(dat.trait) <- c("FFQ", "Diaries")
  
  # Compute Spearman correlation coefficient (Rs)
  Rs <- round(cor(dat.trait[,1], dat.trait[,2],
                  use = "pairwise.complete.obs", method = "spearman"), 2)
  
  # Compute Spearman p-value using corr.test from psych package
  pvalue <- round(corr.test(dat.trait[,1], dat.trait[,2],
                        method = "spearman",
                        use = "pairwise.complete.obs",
                        adjust = "none")$p, 7)
  
  # Combine trait name, correlation, and p-value into a dataframe
  pw_df <- data.frame(trait, Rs, pvalue)
  
  # Append results to cumulative p-value dataframe
  all_pvalues <- rbind(all_pvalues, pw_df)
  
  # Determine plot limits based on observed data
  lim <- unlist(na.omit(dat.trait))
  lim1 <- min(lim)
  lim2 <- max(lim)
  
  # Create scatterplot with regression line and confidence interval
  
  plot <- ggplot(dat.trait, aes(FFQ, Diaries)) +
    geom_point(colour = "royalblue1") +               # points
    xlim(lim1, lim2) + ylim(lim1, lim2) +            # consistent axis limits
    geom_smooth(method = "lm", se = TRUE) +          # regression line with CI
    labs(
      title = trait,                                  # main title: trait name
      subtitle = paste0("Spearman r = ", Rs, "\n", "p-value = ", pvalue)  # subtitle with stats
    ) +
    theme(
      plot.title = element_text(color = "royalblue4",   # change title color
                                face = "bold", 
                                size = 14),
      plot.subtitle = element_text(color = "black",  # keep subtitle black
                                   size = 12)
    )
  
  # Add marginal histograms to scatterplot
  plot_w_hist <- ggMarginal(plot, type = "histogram", fill = "orange")
  
  # Store plot in the list
  plots[[trait]] <- plot_w_hist 
}

# ============================================================================ #
#  EXPORT RESULTS
# ============================================================================ #

# Apply Bonferroni correction for multiple comparisons
all_pvalues$p.adj.Bonf <- p.adjust(all_pvalues$pvalue, method = "bonferroni")

# Export correlation coefficients and p-values to CSV
write.csv(all_pvalues, paste0(label, "_Spearman_pvalues.csv"), row.names = FALSE)

# Save all plots into a single JPEG file
jpeg(file = paste0(label, "_correlations.jpeg"),
     width = 700, height = 1600, units = "px", res = 100)

grid.arrange(
  grobs = plots,
  ncol = 2,
  top = textGrob(
    paste0("Spearman correlations", "\n"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)

# Close JPEG device
dev.off()

