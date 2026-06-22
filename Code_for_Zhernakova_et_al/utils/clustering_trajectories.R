library(factoextra)
all_prots_traj <- read.delim(paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_235_prots.txt"), as.is = T, check.names = F, sep = "\t", row.names = 1)

non_linear_prots <- signif[signif$gam_edf_round > 1, "prot"]
linear_prots <- signif[signif$gam_edf_round == 1, "prot"]

n_points = 100
points <- as.character(seq(1,4,length.out = n_points))
fitted_matrix = all_prots_traj [non_linear_prots, points]
scaled_matrix <- t(scale(t(fitted_matrix)))

colnames(scaled_matrix) <- points
colnames(fitted_matrix) <- points
# Correlation distance: (1 - Pearson Correlation)
dist_mat <- as.dist(1 - cor(t(scaled_matrix)))

hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 4)

write.table(as.data.frame(clusters), file = paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.k5.txt"),quote = F, sep = "\t")


# Linear trajectories

fitted_matrix_lin = all_prots_traj [linear_prots,points]
scaled_matrix_lin <- t(scale(t(fitted_matrix_lin)))

colnames(scaled_matrix_lin) <- points
colnames(fitted_matrix_lin) <- points
# Correlation distance: (1 - Pearson Correlation)
dist_mat_lin <- as.dist(1 - cor(t(scaled_matrix_lin)))

hc_lin <- hclust(dist_mat_lin, method = "ward.D2")
clusters_lin <- cutree(hc_lin, k = 2)


write.table(as.data.frame(clusters_lin), file = paste0(out_basedir, "trajectories_gam/gam_linear_clustering_inv_correl_dist.k2.txt"),quote = F, sep = "\t")

clusters <- clusters + 2
clusters_combined <- c(clusters_lin, clusters)

fitted_matrix_combined <- rbind(fitted_matrix, fitted_matrix_lin)
scaled_matrix_combined <- rbind(scaled_matrix, scaled_matrix_lin)

write.table(as.data.frame(fitted_matrix_combined), file = paste0(out_basedir, "trajectories_gam/fitted_matrix_combined.txt"),quote = F, sep = "\t", row.names = TRUE, col.names = NA)
write.table(as.data.frame(scaled_matrix_combined), file = paste0(out_basedir, "trajectories_gam/scaled_matrix_combined.txt"),quote = F, sep = "\t", row.names = TRUE, col.names = NA)

clusters_combined <- as.data.frame(clusters_combined) %>%
  rownames_to_column("protein")
fitted_matrix_combined <- as.data.frame(fitted_matrix_combined) %>%
  rownames_to_column("protein")
scaled_matrix_combined <- as.data.frame(scaled_matrix_combined) %>%
  rownames_to_column("protein")

plot_data <-  left_join(scaled_matrix_combined, clusters_combined, by = "protein") %>%
  tidyr::pivot_longer(cols = -c(protein, clusters_combined), 
                      names_to = "time_point", 
                      values_to = "value") %>%
  mutate(time_point = as.numeric(gsub("X", "", time_point)))

plot_data_reordered <- plot_data %>%
  mutate(clusters_combined = recode(clusters_combined,
                                    `4` = 5,
                                    `5` = 6,
                                    `6` = 4))
write.table(plot_data_reordered, file = paste0(out_basedir, "trajectories_gam/plot_data_reordered.txt"),quote = F, sep = "\t")

pdf(paste0(out_basedir, "trajectories_gam/gam_clustering_inv_correl_dist.hcclust.k4.combined.pdf"), height = 6, width = 4)
# Plot trajectories by cluster
ggplot(plot_data_reordered, aes(x = time_point, y = value, group = protein,)) +
  geom_line(alpha = 0.2, color = 'dodgerblue3') +
  stat_summary(aes(group = clusters_combined), fun = mean, geom = "line", linewidth = 1.5, color = 'dodgerblue4') +
  facet_wrap(~ clusters_combined,ncol = 3 ) +
  labs(x = "Phase", y = "Scaled GAM fitted protein trajectories") +
  theme_minimal()


hclust_cor <- function(x, k) {
  dist_mat <- as.dist(1 - cor(t(x))) 
  hc <- hclust(dist_mat, method = "ward.D2")
  list(data = x, cluster = cutree(hc, k = k))
}

pam_cor <- function(x, k) {
  dist_mat <- as.dist(1 - cor(t(x))) 
  pam_res <- cluster::pam(dist_mat, diss = T, k = k)
  list(cluster = pam_res$clustering)
}



# 2. Run Silhouette method using this custom function
set.seed(123)
p1 <- fviz_nbclust(scaled_matrix, 
                   FUNcluster = hclust_cor, 
                   method = "silhouette")

p2 <- fviz_nbclust(scaled_matrix, 
                   FUNcluster = hclust_cor, 
                   method = "wss")



print(p1+p2)
dev.off()





#####
# clustering of individual trajectories
#####
library(factoextra)
set.seed(123) 



prot_name = 'PROK1'

p = cluster_individual_trajectories_per_prot(d_wide_adj_covar, prot_name, n_clusters = 3, all_prots_traj[prot_name,], center =T)
p$p_traj
p$p_clust

pdf(paste0(out_basedir, "trajectories_gam/PROK1_individual_traj.pdf"), height = 5, width = 4)
prot_name = 'PROK1'
p = cluster_individual_trajectories_per_prot(d_wide_adj_covar, prot_name, n_clusters = 3, all_prots_traj[prot_name,], center =T)
p$p_traj 
dev.off()

cluster_individual_trajectories_per_prot <- function(d_wide_adj_covar, prot, n_clusters = 3, gam_pred = NULL, center = T){
  tmp_complete <- d_wide_adj_covar %>%
    filter(phase %in% all_phases) %>%
    group_by(ID) %>%
    filter(n_distinct(phase) == 4 & !any(is.na(.data[[prot_name]]))) %>%
    ungroup()
  tmp_complete$phase <- factor(tmp_complete$phase, levels = all_phases)
  
  wide_data <- tmp_complete %>%
    dplyr::select(ID, phase, all_of(prot)) %>%
    pivot_wider(names_from = phase, values_from = all_of(prot)) %>%
    dplyr::select(ID, all_of(all_phases))
  
  # Extract the numeric matrix (exclude ID column)
  mat <- as.matrix(wide_data[, -1])
  row.names(mat) <- wide_data$ID
  # center per ID
  if (center) {
    mat_centered <- t(apply(mat, 1, function(row) row - mean(row)))
  } else {
    mat_centered = mat
  }
  
  # Elbow + silhouette + gap statistic in one plot
  p1 <- fviz_nbclust(mat_centered, kmeans, method = "wss")   # elbow
  p2 <- fviz_nbclust(mat_centered, kmeans, method = "silhouette", print.summary = T)
  #fviz_nbclust(mat_centered, kmeans, method = "gap_stat", nboot = 50)
  
  nclusters = n_clusters
  km_centered <- kmeans(mat_centered, centers = nclusters, nstart = 25)
  mat_centered <- as.data.frame(mat_centered)
  mat_centered$cluster <- as.numeric(km_centered$cluster)
  
  plot_data <- mat_centered %>%
    rownames_to_column("ID") %>%
    dplyr::select(ID, F, O, EL, LL, cluster) %>%
    pivot_longer(cols = c(F, O, EL,LL),
                 names_to = "phase",
                 values_to = "value") %>%
    mutate(phase = factor(phase, levels = all_phases), phase_num = as.numeric(phase)) 
  
  # Mean trajectory per cluster
  cluster_means <- plot_data %>%
    group_by(cluster, phase_num) %>%
    summarise(mean = mean(value, na.rm = TRUE), .groups = "drop")
  
  if (!is.null(gam_pred)) {
    gam_pred = as.data.frame(t(gam_pred)) %>%
      rownames_to_column("phase_num") %>%
      mutate(phase_num = as.numeric(phase_num))
    colnames(gam_pred)[2] <- "GAM_curve"
    if (center) gam_pred$GAM_curve <- gam_pred$GAM_curve - mean(gam_pred$GAM_curve)
    p_traj <- ggplot(plot_data, aes(x = phase_num, y = value, group = ID, color = factor(cluster))) +
      geom_line(alpha = 0.4) +
      geom_line(data = cluster_means, aes(x = phase_num, y = mean, group = cluster, color = factor(cluster)),
                linewidth = 1, inherit.aes = FALSE) +
      geom_line(data = gam_pred, aes(x = phase_num, y = GAM_curve),
                linewidth = 1, inherit.aes = FALSE)
    
  } else {
    p_traj <- ggplot(plot_data, aes(x = phase, y = value, group = ID, color = factor(cluster))) +
      geom_line(alpha = 0.4) +
      geom_line(data = cluster_means, aes(x = phase, y = mean, group = cluster, color = factor(cluster)),
                linewidth = 1, inherit.aes = FALSE) 
  }
  p_traj <- p_traj +
    labs(title = paste0(prot, " trajectories by cluster"),
         color = "Cluster") +
    ylab("adjusted centered abundance") +
    theme_minimal() +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    ) +
    theme(legend.position = "none") +
    scale_color_manual(values = c( "#E69F00", "#56B4E9", "#009E73", "dodgerblue4", "#FF6B6B"))
  return(list(p_traj = p_traj, p_clust = (p1 + p2)))
}





