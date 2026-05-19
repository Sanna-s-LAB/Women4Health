# PAOLA FORABOSCO 28/04/2026

# ============================================================================ #
#  DESCRIPTION
# ============================================================================ #

# Generates boxplots for nutritional traits / diet scores across cohorts,
# correcting for Age and BMI, and performs pairwise Wilcoxon rank-sum tests.
# Also performs Kruskal–Wallis tests for global differences.
# 
# All pairwise p-values (raw and Bonferroni-corrected) are exported to CSV.
# Boxplots with significance annotations are saved as JPEG.

# Note, for FFQ analysis a variable “first30” is added in the regression model

# INPUT FILES:
# Diaries_TM / FFQ / Diet_Scores
# containing nutritional traits from Daily Diaries/FFQ or Diet scores
# The columns sample/Age/BMI must be present

# ============================================================================ #
#  INITIALIZATION
# ============================================================================ #

# Load required libraries

library(ggplot2)    # For creating boxplots
library(dplyr)      # For data manipulation (group_by, mutate, summarise)
library(ggpubr)     # For stat_compare_means() to add significance annotations
library(gridExtra)  # For arranging multiple plots into a single grid

# Data label - change to match dataset ("ALL_Diaries_TM", "FFQ", "Diet_Scores")
# label <- "ALL_Diaries_TM" 
# label <- "FFQ" 
label <- "Diet_Scores" 

# Read data
dat <- read.csv(paste0(label, ".csv"), header=TRUE)

# Traits to analyze  - change to match dataset 
#traits <- c("Calories", "Carbohydrates", "Fats", "Proteins", "Cholesterol",
#             "Sodium", "Sugar", "Fiber")

traits <- c("aMED_total","HEI_total","DASH_total","AHEI_total")


# Column name for cohort/category variable
Category.Name <- "sample"  

# Initialize dataframes to store p-values and list to store plots
all_pvalues <- data.frame()
kw_pvalues <- data.frame()
plots <- list()          

# Quick data summaries
summary(dat)
by(dat, dat$sample, summary)

# ============================================================================ #
#  BOXPLOT GENERATION & PAIRWISE TESTS
# ============================================================================ #

# Loop over each trait to generate boxplots and perform tests
for (trait in traits) {
  
  # --------------------- #
  # Correct for Age & BMI #
  # --------------------- #
  
  lm_model <- lm(dat[,trait] ~ Age + BMI , data = dat, na.action = na.exclude)
  
# For FFQ correction of the first 30 participants
# lm_model <- lm(dat[,trait] ~ Age + BMI + first30, data = dat, na.action = na.exclude)
  
  res <- data.frame(resid(lm_model))
  colnames(res) <- trait
  res$category <- dat[,Category.Name]
  res$ID <- dat$ID
  
  # Prepare dataframe for plotting and testing
  data <- data.frame(
    category = factor(res$category),
    value = res[, trait],
    label = res$ID
  )
  
  # Remove missing values
  data <- na.omit(data)
  
  # ---------------------------- #
  # Detect outliers per category #
  # ---------------------------- #
  
  outliers <- data %>%
    group_by(category) %>%
    mutate(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      is_outlier = value < (Q1 - 1.5 * IQR) | value > (Q3 + 1.5 * IQR)
    ) %>%
    ungroup() %>%
    select(category, value, is_outlier)
  
  max_val <- max(data$value, na.rm = TRUE)
  
  # ------- #
  # Boxplot #
  # ------- #
  
  n_labels <- data %>%
    group_by(category) %>%
    summarise(n = sum(!is.na(value))) %>%
    mutate(label = paste0(category, " (N=", n, ")"))
  
  
  # Boxplot generation
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
  
  
  # --------------------------------- #
  # Pairwise Wilcoxon tests (raw p)  #
  # --------------------------------- #
  
  pw_test <- pairwise.wilcox.test(data$value, data$category, 
                                  p.adjust.method = "none", paired = FALSE, exact = FALSE)
  
  # Convert pairwise p-values matrix to dataframe for export
  pw_df <- as.data.frame(as.table(pw_test$p.value))
  pw_df <- pw_df[-3,]  # Remove self comparisons
  colnames(pw_df) <- c("sample1", "sample2", "pvalue")
  pw_df$trait <- trait
  pw_df <- pw_df[, c("trait", "sample1", "sample2", "pvalue")]
  
  all_pvalues <- rbind(all_pvalues, pw_df)
  
  # Prepare pairwise comparisons for significance bars on plot
  categories <- unique(as.character(data$category))
  comparisons <- combn(categories, 2, simplify = FALSE)
  
  # Add significance bars to boxplot
  bpp <- bp + stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.signif",
    label.y = seq(from = max_val + max_val/5, by = max_val/5, length.out = length(comparisons))
  )
  
  plots[[trait]] <- bpp
  
  # ---------------------- #
  # Kruskal–Wallis test    #
  # ---------------------- #
  
  kw_test <- kruskal.test(data$value ~ data$category)
  
  kw_df <- data.frame(
    trait = trait,
    statistic = kw_test$statistic,
    pvalue = kw_test$p.value
  )
  
  kw_pvalues <- rbind(kw_pvalues, kw_df)
}

# ============================================================================ #
#  EXPORT ALL P-VALUES WITH MULTIPLE TESTING CORRECTION
# ============================================================================ #

all_pvalues$p.adj.Bonf <- p.adjust(all_pvalues$pvalue, method = "bonferroni")
kw_pvalues$p.adj.Bonf <- p.adjust(kw_pvalues$pvalue, method = "bonferroni")

write.csv(all_pvalues, paste0(label,"_pairwise_pvalues_all.csv"), row.names = FALSE)
write.csv(kw_pvalues, paste0(label,"_kruskal_wallis_pvalues_all.csv"), row.names = FALSE)

# ============================================================================ #
#  SAVE PLOTS TO FILE
# ============================================================================ #

# NOTE: height = 1400 for 8 plots, height = 800 for 4 plots 
jpeg(file = paste0(label,"_samples.jpeg"), width=1000, height=800, units="px", res=100)
print(grid.arrange(grobs=plots, ncol=2))
dev.off()



