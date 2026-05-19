# Microbiome scripts
# by: Valeria Lo Faro, IRGB_CA
# ==============================================


## Visualizes microbiome compositional differences across phases
## using mean relative abundance and pie charts ##

# ------------------------------------------------------------
### LOAD LIBRARIES
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(data.table)

# ------------------------------------------------------------
### INPUT DATA
# ------------------------------------------------------------
abundanceFILE <- "MGS_FECAL_metaphlan.rds"

# Read data
inDF <- readRDS(abundanceFILE)

# ------------------------------------------------------------
### DEFINE GROUPS (PHASES)
# ------------------------------------------------------------
group_definitions <- data.frame(
  Name    = c("Follicular", "Ovulatory", "Early luteal", "Late luteal"),
  Pattern = c(".*_1", ".*_2", ".*_3", ".*_4"),
  Color   = c("#f38120", "#c72249", "#a6c1d6", "#a32d93")
)

# ------------------------------------------------------------
### PREPARE DATA
# ------------------------------------------------------------
# Select species features
colsFeatures <- grep("s__", colnames(inDF), value = TRUE)

# Keep relevant columns
inDF <- inDF[, c("ID", "phase", colsFeatures)]

# Create unique sample ID
inDF$sample_ID <- paste0(inDF$ID, "_", inDF$phase)

# Average per sample_ID
inDF <- inDF %>%
  group_by(sample_ID) %>%
  summarise(across(all_of(colsFeatures), mean, na.rm = TRUE)) %>%
  ungroup()

# ------------------------------------------------------------
### TRANSPOSE TO TAXA x SAMPLES
# ------------------------------------------------------------
tdata <- as.data.frame(t(inDF))
colnames(tdata) <- tdata[1, ]
tdata <- tdata[-1, ]

tdata <- tdata %>%
  mutate(across(everything(), as.numeric))

tdata$Lineage <- rownames(tdata)

# ------------------------------------------------------------
### EXTRACT SPECIES NAMES
# ------------------------------------------------------------
tdata$Taxon <- sapply(strsplit(tdata$Lineage, "\\|"), function(x) {
  sp <- grep("s__", x, value = TRUE)
  if (length(sp) > 0) sp else NA
})

tdata <- tdata[!is.na(tdata$Taxon), ]

# ------------------------------------------------------------
### DEFINE SAMPLE GROUPS
# ------------------------------------------------------------
samples_list <- lapply(group_definitions$Pattern, function(p) {
  grep(p, colnames(tdata), value = TRUE)
})

# ------------------------------------------------------------
### COMPUTE GROUP MEAN ABUNDANCE
# ------------------------------------------------------------
group_means <- tdata %>% select(Taxon)

for (i in seq_len(nrow(group_definitions))) {

  group_name <- group_definitions$Name[i]
  sample_cols <- samples_list[[i]]

  if (length(sample_cols) > 0) {
    group_means[[group_name]] <- rowMeans(
      tdata[, sample_cols, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    group_means[[group_name]] <- NA
  }
}

# ------------------------------------------------------------
### RESHAPE DATA
# ------------------------------------------------------------
long_data <- group_means %>%
  pivot_longer(-Taxon, names_to = "Group", values_to = "Mean") %>%
  filter(Mean > 0)

# Total abundance per group
group_totals <- long_data %>%
  group_by(Group) %>%
  summarise(TotalAbundance = sum(Mean), .groups = "drop")

# ------------------------------------------------------------
### SELECT TOP TAXA + OTHER
# ------------------------------------------------------------
top_taxa <- long_data %>%
  group_by(Group) %>%
  slice_max(order_by = Mean, n = 10, with_ties = FALSE) %>%
  ungroup()

others <- long_data %>%
  filter(!Taxon %in% top_taxa$Taxon) %>%
  group_by(Group) %>%
  summarise(Mean = sum(Mean), .groups = "drop") %>%
  mutate(Taxon = "Other")

plot_data <- bind_rows(top_taxa, others)

# Clean names
plot_data <- plot_data %>%
  mutate(
    Taxon_clean = ifelse(Taxon == "Other", "Other", gsub("^.*__", "", Taxon))
  )

# Add percentages
plot_data <- plot_data %>%
  left_join(group_totals, by = "Group") %>%
  mutate(Percentage = Mean / TotalAbundance * 100)

# ------------------------------------------------------------
### COLOR PALETTE
# ------------------------------------------------------------
taxa_order <- plot_data %>%
  group_by(Taxon_clean) %>%
  summarise(TotalMean = sum(Mean), .groups = "drop") %>%
  arrange(desc(TotalMean)) %>%
  pull(Taxon_clean)

taxa_order <- unique(c(taxa_order, "Other"))

palette_colors <- colorRampPalette(
  brewer.pal(min(12, length(taxa_order)), "Set3")
)(length(taxa_order))

palette_colors <- setNames(palette_colors, taxa_order)

# ------------------------------------------------------------
### PLOT FUNCTION
# ------------------------------------------------------------
plot_pie <- function(df) {

  df <- df %>%
    arrange(desc(Mean)) %>%
    mutate(
      label = paste0(Taxon_clean, " (", round(Percentage, 1), "%)"),
      label = factor(label, levels = unique(label))
    )

  colors_subset <- palette_colors[unique(df$Taxon_clean)]
  names(colors_subset) <- levels(df$label)

  ggplot(df, aes(x = "", y = Mean, fill = label)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors_subset) +
    labs(
      title = unique(df$Group),
      fill = "Taxon"
    ) +
    theme_void() +
    theme(
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# ------------------------------------------------------------
### GENERATE AND OUTPUT PLOTS
# ------------------------------------------------------------
plot_list <- plot_data %>%
  group_split(Group) %>%
  setNames(unique(plot_data$Group))

pie_plots <- lapply(plot_list, plot_pie)

# Save plots
# pdf("composition_pies.pdf", width = 12, height = 8)
for (i in seq_along(pie_plots)) {
  print(pie_plots[[i]])
}
# dev.off()
