# PAOLA FORABOSCO 28/04/2026

# ============================================================================ #
#   DESCRIPTION
# ============================================================================ #
# This script performs compositional analysis, PCA, and k-means clustering
# on macronutrient data (Carbohydrates, Fats, Proteins) from a CSV file.
#
# Workflow:
# 1. Calculate calories and proportions per macronutrient, apply CLR transform
# 2. Perform PCA and plot biplots of individuals and variables
# 3. Determine clusters using k-means and visualize with PCA biplot + clusters
# 4. Compute mean macronutrient values per cluster and save results
#
# Inputs files columns: ID, sample, Carbohydrates, Fats, Proteins
# Outputs:
# - PCA coordinates + cluster labels
# - Cluster mean macronutrients + sample counts
# - JPEG plots (PCA and Kmeans clusters)

# ============================================================================ #
#  INITIALIZATION
# ============================================================================ #

# Load required libraries
library(dplyr)          # Data manipulation
library(tidyr)          # Data tidying
library(compositions)   # CLR transformation for compositional data
library(ggplot2)        # Plotting
library(class)          # KNN function
library(caret)          # Cross-validation and KNN tuning
library(factoextra)     # PCA visualization
library(FactoMineR)     # PCA and multivariate analysis

# ------------------------------- #
# Function: prepare_dataset_macro
# ------------------------------- #

prepare_dataset_macro <- function(df) {
  # Calculate calories contributed by each macronutrient
  df <- df %>%
    mutate(
      kcal_carb = Carbohydrates * 4,
      kcal_fat  = Fats * 9,
      kcal_prot = Proteins * 4,
      kcal_tot_macro = kcal_carb + kcal_fat + kcal_prot,
      # Proportion of each macronutrient in total calories
      prop_carb = kcal_carb / kcal_tot_macro,
      prop_fat  = kcal_fat / kcal_tot_macro,
      prop_prot = kcal_prot / kcal_tot_macro
    )
  
  # Select compositional columns and perform CLR transformation
  compo <- df[, c("prop_carb", "prop_fat", "prop_prot")]
  compo_clr <- clr(as.matrix(compo))
  compo_clr <- data.frame(compo_clr)
  colnames(compo_clr) <- c("Carbohydrates", "Fats", "Proteins")
  rownames(compo_clr) <- df$ID
  
  return(compo_clr)
}

# ------------- #
# LOAD DATA
# ------------- #

label <- "Diaries_TM" # FFQ
DATA <- read.csv(paste0(label, ".csv"), header = TRUE)
rownames(DATA) <- DATA$ID
dat <- prepare_dataset_macro(DATA)  # CLR-transformed dataset

# ============================================================================ #
#  PCA ANALYSIS
# ============================================================================ #

# Perform PCA
pca <- PCA(dat, graph = FALSE, scale = TRUE)

# Extract variance explained
eig <- get_eigenvalue(pca)

# Prepare individual coordinates for plotting and clustering
ind.coord <- as.data.frame(pca$ind$coord)[, 1:3]
colnames(ind.coord) <- c("PCA1", "PCA2", "PCA3")
ind.coord$sample <- as.factor(DATA$sample)

# ------------------- #
#  INITIAL PCA BIPLOT  #
# ------------------- #

X <- 1  # Select PCA 1
Y <- 2  # Select PCA 2
x_var <- round(eig[X, 2], 0)  # Variance explained by PC1
y_var <- round(eig[Y, 2], 0)  # Variance explained by PC2

pca_biplot <- fviz_pca_biplot(
  pca,
  axes = c(X, Y),
  col.var = "black",       # Variable arrows in black
  col.ind = ind.coord$sample,
  palette = c("royalblue3", "magenta", "cyan3"),
  label = "var",
  repel = TRUE,
  pointshape = 16,
  pointsize = 3,
  labelsize = 4,
  alpha.var = 1
) +
  labs(
    x = paste0("PCA", X, " (", x_var, "%)"),
    y = paste0("PCA", Y, " (", y_var, "%)"),
    color = NULL,
    shape = NULL
  ) +
  guides(shape = "none") +
  ggtitle(paste0(label, " - PCA Biplot with Samples"))


# --------------------- #
# Save the plot as JPEG #
# --------------------- #
jpeg(
  file = paste0(label, "_PCA_sample.jpeg"),  # File name
  width = 650, height = 650, units = "px", res = 100      # Image size and resolution
)
print(pca_biplot)  # Print the plot to the JPEG device
dev.off()   # Close the device


# ============================================================================ #
#  K-MEANS CLUSTERING
# ============================================================================ #

# Determine optimal clusters using the elbow method
fviz_nbclust(ind.coord[, 1:3], kmeans, method = "wss")

# Perform k-means clustering with 3 clusters
set.seed(123)
kmeans_clust <- kmeans(ind.coord[, 1:3], centers = 3)
kmeans_labels <- as.factor(kmeans_clust$cluster)

# Save cluster assignments with PCA coordinates
ind.coord$cluster <- kmeans_labels

# Add ID column from rownames
ind.coord$ID <- rownames(ind.coord)

# Reorder columns so ID is first
ind.coord <- ind.coord[, c("ID", setdiff(names(ind.coord), "ID"))]

# Save to CSV
write.csv(ind.coord, paste0(label, "_clusters.csv"), row.names = FALSE)

# -------------------------- #
#  PCA BIPLOT WITH CLUSTERS
# -------------------------- #

levels(kmeans_labels) <- c("1", "2", "3")
# levels(kmeans_labels) <- c("High-Fat", "Balanced", "High-Carb")

cluster_biplot <- fviz_pca_biplot(
  pca,
  axes = c(X, Y),
  col.var = "black",
  habillage = kmeans_labels,
  col.ind = kmeans_labels,
  palette = c("royalblue3", "magenta", "cyan3"),
  pointshape = 16,
  pointsize = 3,
  label = "var",
  addEllipses = TRUE,
  ellipse.type = "convex",
  labelsize = 3,
  alpha.var = 1
) +
  labs(
    x = paste0("PCA", X, " (", x_var, "%)"),
    y = paste0("PCA", Y, " (", y_var, "%)")
  ) +
  guides(color = guide_legend(override.aes = list(label = NULL)), shape = "none") +
  ggtitle(paste0(label, " - PCA Biplot & K-means Clusters"))


# --------------------- #
# Save the plot as JPEG #
# --------------------- #
jpeg(
  file = paste0(label, "_KMEANS_clusters.jpeg"),  # File name
  width = 650, height = 650, units = "px", res = 100      # Image size and resolution
)
print(cluster_biplot)  # Print the plot to the JPEG device
dev.off()   # Close the device


# ============================================================================ #
#  CALCULATE MEANS PER CLUSTER
# ============================================================================ #

dat.ORIG <- DATA %>%
  mutate(
    kcal_carb = Carbohydrates * 4,
    kcal_fat = Fats * 9,
    kcal_prot = Proteins * 4,
    kcal_tot_macro = kcal_carb + kcal_fat + kcal_prot,
    prop_carb = kcal_carb / kcal_tot_macro,
    prop_fat = kcal_fat / kcal_tot_macro,
    prop_prot = kcal_prot / kcal_tot_macro,
    kmeans_clust = as.factor(kmeans_clust$cluster)
  )

# Compute mean values per cluster with sample count
dat.means <- dat.ORIG %>%
  group_by(kmeans_clust) %>%
  summarise(
    across(where(is.numeric), ~ round(mean(.), 2)),  # Mean of numeric columns
    N = n()                                          # Number of samples per cluster
  ) %>%
  relocate(N, .after = kmeans_clust)               # Move N next to cluster

write.csv(dat.means, paste0(label, "_clusters_means.csv"), row.names = FALSE)





