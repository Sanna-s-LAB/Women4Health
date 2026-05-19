# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Computes correlations between PCoA axes (beta diversity structure)
## and FSH levels stratified by phase, and generates multi-panel plots ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
abundanceFILE <- "species.MGS_FECAL_merged_abundance_table_clr.rds"

inDF <- readRDS(abundanceFILE)

# Remove samples without phase
inDF <- inDF[!is.na(inDF$phase), ]

# ------------------------------------------------------------
### DEFINE VARIABLES
# ------------------------------------------------------------
# PCoA axes to test
pcoa_vars <- c("PCOA1", "PCOA2", "PCOA3")

# ------------------------------------------------------------
### COMPUTE CORRELATIONS + GENERATE PLOTS
# ------------------------------------------------------------
plot_list <- list()
k <- 1

for (ph in sort(unique(inDF$phase))) {

  df <- inDF[inDF$phase == ph, ]

  for (pc in pcoa_vars) {

    x <- df[[pc]]
    y <- df$FSH

    idx <- complete.cases(x, y)

    # Skip if insufficient data
    if (sum(idx) < 5) next

    # Spearman correlation
    cor_res <- cor.test(
      x[idx],
      y[idx],
      method = "spearman",
      exact = FALSE
    )

    # Plot
    plot_list[[k]] <- ggplot(df, aes(x = .data[[pc]], y = FSH)) +
      geom_point(color = "dodgerblue", alpha = 0.7) +
      geom_smooth(method = "lm", color = "red", linewidth = 0.6) +
      labs(
        x = pc,
        y = "FSH",
        title = paste0("Phase ", ph),
        subtitle = paste0(
          "Spearman rho = ", round(cor_res$estimate, 3),
          " (p = ", signif(cor_res$p.value, 3), ")"
        )
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9)
      )

    k <- k + 1
  }
}

# ------------------------------------------------------------
### OUTPUT PLOTS
# ------------------------------------------------------------
# pdf("PCOA_FSH_correlations.pdf", width = 8, height = 12)

for (i in seq(1, length(plot_list), by = 6)) {

  grid.arrange(
    grobs = plot_list[i:min(i + 5, length(plot_list))],
    ncol = 2,
    nrow = 3
  )
}

# dev.off()





