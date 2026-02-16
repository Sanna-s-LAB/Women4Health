library(factoextra)
set.seed(123)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
out_basedir <- "results12/intensity_shared_prots_261125/"
#
# Cluster non-linear trajectories
#
# Read the pre-saved trajectories (rows names - proteins, column names - 100 time points in the interval from 1 to 4, values - predicted GAM values)
all_prots_traj <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_93_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)

n_points = 100
fitted_matrix = all_prots_traj [,1:n_points]

# scale the trajectories
scaled_matrix <- t(apply(fitted_matrix, 1, function(x) (x - min(x)) / (max(x) - min(x))))
colnames(scaled_matrix) <- seq(1,4, length.out=n_points)

# Calculate the correlation distance: (1 - Pearson Correlation)
dist_mat <- as.dist(1 - cor(t(scaled_matrix)))

# Perform clustering 

# check the optimal number of clusters
hclust_cor <- function(x, k) {
  dist_mat <- as.dist(1 - cor(t(x))) 
  hc <- hclust(dist_mat, method = "ward.D2")
  list(data = x, cluster = cutree(hc, k = k))
}
fviz_nbclust(scaled_matrix, FUNcluster = hclust_cor, method = "silhouette")
fviz_nbclust(scaled_matrix, FUNcluster = hclust_cor, method = "wss")
fviz_cluster(hclust_cor(scaled_matrix, 6), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal())

# Cluster using k = 6
hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 6)
write.table(as.data.frame(clusters), file = paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.k6.txt"),quote = F, sep = "\t")


#
# Cluster linear trajectories
#
all_prots_traj_lin <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_70_linear_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)

n_points = 100
fitted_matrix_lin = all_prots_traj_lin [,1:n_points]

# scale trajectories
scaled_matrix_lin <- t(apply(fitted_matrix_lin, 1, function(x) (x - min(x)) / (max(x) - min(x))))
colnames(scaled_matrix_lin) <- seq(1,4, length.out=n_points)

# Calculate the correlation distance: (1 - Pearson Correlation)
dist_mat_lin <- as.dist(1 - cor(t(scaled_matrix_lin)))

# Cluster into 2 (linearly increasing and decreasing) clusters
hc_lin <- hclust(dist_mat_lin, method = "ward.D2")
clusters_lin <- cutree(hc_lin, k = 2)

write.table(as.data.frame(clusters_lin), file = paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.k2.txt"),quote = F, sep = "\t")

#
# Combine linear and non-linear clusters for plotting
#
clusters <- clusters + 2
clusters_combined <- c(clusters_lin, clusters)

fitted_matrix_combined <- rbind(fitted_matrix, fitted_matrix_lin)

write.table(as.data.frame(clusters_combined), file = paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.combined.txt"),quote = F, sep = "\t")
write.table(as.data.frame(fitted_matrix_combined), file = paste0(out_basedir, "trajectories_gam/gam_linear_trajectories_163signif.txt"),quote = F, sep = "\t", row.names = TRUE, col.names = NA)

clusters_combined <- as.data.frame(clusters_combined) %>%
  rownames_to_column("protein")
fitted_matrix_combined <- as.data.frame(fitted_matrix_combined) %>%
  rownames_to_column("protein")

plot_data <-  left_join(fitted_matrix_combined, clusters_combined, by = "protein") %>%
  tidyr::pivot_longer(cols = -c(protein, clusters_combined), 
                      names_to = "time_point", 
                      values_to = "value") %>%
  mutate(time_point = as.numeric(gsub("X", "", time_point)))

#
# Plot trajectories by cluster
#
pdf(paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.combined.pdf"))

ggplot(plot_data, aes(x = time_point, y = value, group = prot,)) +
  geom_line(alpha = 0.3, color = 'dodgerblue3') +
  stat_summary(aes(group = clusters_combined), fun = mean, geom = "line", size = 1.5, color = 'dodgerblue4') +
  facet_wrap(~ clusters_combined) +
  labs(x = "Phase", y = "Scaled GAM fitted protein trajectories") +
  theme_minimal()

dev.off()