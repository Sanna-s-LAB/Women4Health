library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(pheatmap)
library(patchwork)

set.seed(123)
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
script_dir <- "/Users/Dasha/work/Sardinia/W4H/olink/scripts/"
out_basedir <- "results12/intensity_shared_prots_261125/"

source(paste0(script_dir, "utils/utility_functions.R"))


################################################################################
# Read in the data and format it
################################################################################

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.txt", as.is = T, check.names = F, sep = "\t")

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    dplyr::select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates <- read.delim("results12/covariates_olink_batch12.phase_avg.txt", sep = "\t", check.names = F, as.is = T)

if (! "ID" %in% colnames(covariates)){
  covariates$ID <- gsub("_.*", "", covariates$SampleID)
  covariates$phase <- gsub(".*_", "", covariates$SampleID)
  covariates <- covariates %>%
    dplyr::select(SampleID, ID, phase, everything())
}
covariates$phase <- relevel(factor(covariates$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

covariates[] <- lapply(covariates, function(col) {
  if (length(unique(col)) < 3) {
    return(factor(col))
  } else {
    return(col)
  }
})


# Make a dataframe with proteins adjusted for all covariates per visit
covariate_names <- c("Age","BMI","batch", "from", "storage_months")
covariates$from <- relevel(as.factor(covariates$from), ref = "X")

shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 50]
d_wide_shared <- d_wide[ ,shared_prots]
all_prots <- colnames(d_wide_shared)[! colnames(d_wide_shared) %in% c("ID", "TP","phase","SampleID")]

d_wide_b2 <- d_wide[d_wide$SampleID %in% covariates[covariates$batch == 'batch2', "SampleID"],]

# remove overlapping prots from the b2 dataset
d_wide_b2 <- d_wide_b2[, !colnames(d_wide_b2) %in% all_prots]

dim(d_wide)
dim(d_wide_b2)
dim(d_wide_shared)

d_wide_full <- d_wide
d_wide <- d_wide_shared

all_phases <- c("F", "O", "EL", "LL")

joined_data <- full_join(covariates, d_wide, by = c("SampleID", "ID", "phase"))
d_wide_adj_covar <- regress_covariates_lmm_phase(d_wide, covariates, covars_longitudinal = T)

joined_data_adj_covar <- full_join(covariates, d_wide_adj_covar, by = c("SampleID") )

write.table(d_wide_adj_covar, file = paste0(out_basedir, "olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.adj_all_covariates.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# ICC for each protein
################################################################################

icc <- data.frame(matrix(nrow = (ncol(d_wide) -4), ncol = 3))
colnames(icc) <- c("prot", "ICC", "R2_TP")
cnt <- 1
for (prot in all_prots){
  res <- get_ICC(d_wide, prot)
  icc[cnt,] <- c(prot, unlist(res))
  cnt <- cnt + 1
}

icc <- na.omit(icc) %>%
  mutate(across(-c( prot), as.numeric)) 
icc <- icc[order(icc$ICC, decreasing = F),]

cat("ICC ranges from", min(icc$ICC), "to", max(icc$ICC), "with a median of", median(icc$ICC), "\n")
ggplot(icc, aes(x=ICC, y =R2_TP)) + geom_point() + theme_minimal() + xlab("ICC (variance explained by ID)") + ylab("Marginal R2 (var explained by phase)")

write.table(icc, file = paste0(out_basedir, "ICC_per_protein.txt"), quote = F, sep = "\t", row.names = FALSE)

pdf(paste0(out_basedir, "plots/ICC_per_protein.pdf"))
p <- ggplot(icc, aes(x = ICC, y = R2_TP)) + 
  geom_point() + 
  theme_minimal() + 
  xlab("ICC (variance explained by Sample ID)") + 
  ylab("Marginal R2 (variance explained by phase)")

# Add marginal boxplots
ggMarginal(p, type = "boxplot", margins = "both", 
           size = 5, fill = 'white', alpha = 0.7)

dev.off()


################################################################################
# Differentially expressed proteins between visits
################################################################################

limma_res_all <- data.frame()
phase_comb <- t(combn(all_phases, 2))
for (i in 1:nrow(phase_comb)){
  tp1 = phase_comb[i,1]
  tp2 = phase_comb[i,2]

  # limma
  limma_res <- run_limma(joined_data, tp1, tp2) %>%
    rownames_to_column(var = 'prot')
  if (nrow(limma_res[limma_res$adj.P.Val < 0.05,]) > 0) {
    limma_res_all <- rbind(limma_res_all, cbind(paste0(tp1, "_", tp2), limma_res))
  }
}
colnames(limma_res_all)[1] <- "phase1_phase2"
limma_res_all$sign <- ifelse(limma_res_all$adj.P.Val < 0.05, T, F)

write.table(limma_res_all, file = paste0(out_basedir, "limma_DEPs.txt"), quote = F, sep = "\t", row.names = FALSE)

# Plot significant DEPs

# stacked barplot
results <- limma_res_all %>%
  mutate(
    Regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

summary_data <- results %>%
  filter(Regulation != "Not significant") %>%  # Exclude non-significant proteins
  group_by(phase1_phase2, Regulation) %>%
  summarize(Count = n(), .groups = "drop") %>%
  mutate(phase1_phase2 = factor(phase1_phase2, 
                                levels = c("F_O", "F_EL", "F_LL", "O_EL", "O_LL", "EL_LL")))

pdf(paste0(out_basedir, "plots/limma_DEPs_barplot.pdf"))
ggplot(summary_data, aes(x = phase1_phase2, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Phase Comparison",
    y = "Number of DEP",
    fill = "Effect direction"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = my_colors[c(3,2)])
dev.off()

limma_res_all$signif_direction <- ifelse(limma_res_all$sign, ifelse(limma_res_all$logFC < 0, -1, 1), 0)

prot_subs <- unique(limma_res_all[limma_res_all$sign == T, "prot"])
limma_heatmap2 <- my_pivot_wider(limma_res_all[limma_res_all$prot %in% prot_subs,], "prot", "phase1_phase2", "logFC")
signif_labels <- my_pivot_wider(limma_res_all[limma_res_all$prot %in% prot_subs,], "prot", "phase1_phase2", "adj.P.Val")
signif_labels <- ifelse(signif_labels < 0.05, "*", "")

max_val <- max(abs(min(limma_heatmap2[row.names(limma_heatmap2) != 'PROK1',])), max(limma_heatmap2[row.names(limma_heatmap2) != 'PROK1',]))
limma_heatmap2[limma_heatmap2 > max_val] <- max_val
breaksList = seq(-max_val, max_val, by = 0.05)
colorList <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
h2 <- pheatmap(limma_heatmap2, display_numbers = signif_labels, fontsize_number = 14, 
              color = colorList, breaks = breaksList, cutree_rows = 4)

clusters <- as.data.frame(cutree(h2$tree_row,4)) %>%
  rownames_to_column("prot")
colnames(clusters)[2] <- "cluster"

pdf(paste0(out_basedir, "plots/limma_DEPs_heatmap.pdf"), width = 6, height = 10)
grid::grid.newpage()
grid::grid.draw(h2$gtable)
dev.off()


# CSF3

pdf(paste0(out_basedir, "plots/limma_DEPs_CSF3.pdf"), width = 5, height = 5)
ggplot(d_wide_adj_covar, aes(x = phase, y = CSF3, group = phase)) + 
  geom_boxplot(width = 0.3, color = my_colors[3], outliers = F) + 
  geom_jitter(, alpha = 0.2, width = 0.1, color = my_colors[3]) + 
  theme_minimal() +
  ylab("CSF3 adjusted levels")
dev.off()


################################################################################
# Protein vs phase GAM and LMM
################################################################################

gam_res_prot_tp <- data.frame(matrix(nrow = length(all_prots), ncol = 7))
colnames(gam_res_prot_tp) <- c("prot", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")
pb <- txtProgressBar(min = 0, max = length(all_prots), style = 3)

cnt <- 1
for (prot in all_prots){
  res_gam <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = F)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(d_wide, prot, covariates)
  gam_res_prot_tp[cnt,] <- c(prot, unlist(res_gam), res_lmm)
  cnt <- cnt + 1
  setTxtProgressBar(pb, cnt)
}
close(pb)

gam_res_prot_tp <- na.omit(gam_res_prot_tp) %>%
  mutate(across(-c( prot), as.numeric)) 

gam_res_prot_tp$gam_BH_pval <- p.adjust(gam_res_prot_tp$gam_pval, method = 'BH')
gam_res_prot_tp$lmm_BH_pval <- p.adjust(gam_res_prot_tp$lmm_pval, method = 'BH')

gam_res_prot_tp <- gam_res_prot_tp[order(gam_res_prot_tp$gam_pval),]
gam_res_prot_tp$gam_edf_round <- round(gam_res_prot_tp$gam_edf)

write.table(gam_res_prot_tp, file = paste0(out_basedir, "prot_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = FALSE)
#gam_res_prot_tp <- read.delim(paste0(out_basedir, "prot_vs_tp_gam.txt"), sep = "\t", as.is = T, check.names = F)
signif <- gam_res_prot_tp[gam_res_prot_tp$gam_BH_pval < 0.05,]

cat ("Number of proteins that change significantly with time:", nrow(signif), "\n")
cat("Of them, the number of proteins with a non-linear change: ", nrow(signif[signif$gam_edf_round > 1,]), "\n")
cat("Of them, the number of proteins showing a significant association with time also in LMMs:", nrow(signif[signif$lmm_BH_pval < 0.05,]))


####################NULL################################################################################
# Cluster longitudinal trajectories: GAM
################################################################################

n_points = 100

# Take only significant non-linear trajectories:
all_sign_prots <- signif[signif$gam_edf_round > 1, "prot"]

prot_trajs <- data.frame(matrix(nrow = length(all_sign_prots) , ncol = n_points))
row.names(prot_trajs) <- all_sign_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)

# Pre-compute all GAM fits once
pb <- txtProgressBar(min = 0, max = length(all_sign_prots), style = 3)
cnt <- 1
for (prot in all_sign_prots) {
  gam_fit <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = T, n_points = n_points)
  prot_trajs[prot,] <- gam_fit$predicted
  cnt <- cnt + 1
  setTxtProgressBar(pb, cnt)
}
close(pb)

write.table(prot_trajs, file = paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_", length(all_sign_prots), "_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)

# Take only significant linear trajectories:
all_sign_prots <- signif[signif$gam_edf_round == 1, "prot"]

prot_trajs <- data.frame(matrix(nrow = length(all_sign_prots) , ncol = n_points))
row.names(prot_trajs) <- all_sign_prots
colnames(prot_trajs) <- seq(1,4, length.out = n_points)

# Pre-compute all GAM fits once
pb <- txtProgressBar(min = 0, max = length(all_sign_prots), style = 3)
cnt <- 1
for (prot in all_sign_prots) {
  gam_fit <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = T, n_points = n_points)
  prot_trajs[prot,] <- gam_fit$predicted
  cnt <- cnt + 1
  setTxtProgressBar(pb, cnt)
}
close(pb)
write.table(prot_trajs, file = paste0(out_basedir, "trajectories_gam/protein_trajectories_gam_", length(all_sign_prots), "_linear_prots.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)

# Run clustering separately:
source(paste0(script_dir, "utils/clustering_trajectories.R"))

################################################################################
# GAM  phenotype vs TP
################################################################################
pheno <- read.delim("../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt", as.is = T, check.names = F, sep = "\t")            

all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
all_phenos <- c("INS", "HOMA_B", "HOMA_IR", "GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")
covariate_names_pheno <- c("from", "Age", "BMI")

if (! "ID" %in% colnames(pheno)){
  pheno$ID <- gsub("_.*", "", pheno$SampleID)
  pheno$phase <- gsub(".*_", "", pheno$SampleID)
  pheno <- pheno %>%
    dplyr::select(SampleID, ID, phase, everything())
}

pheno$phase <- relevel(factor(pheno$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

pheno_adjusted <- regress_covariates_lmm_phase(pheno, subset(covariates, select = -c(storage_months, batch)), covars_longitudinal = T)
pheno_adjusted$TP <- as.numeric(pheno_adjusted$phase)
## Hormone and lipid trajectories
all_phenos_combined <- c(all_hormones, all_phenos)

gam_res_pheno_tp <- data.frame(matrix(nrow = length(all_phenos_combined), ncol = 7))
colnames(gam_res_pheno_tp) <- c("pheno", "gam_pval", "gam_edf", "gam_fval", "n", "n_samples", "lmm_pval")

pheno_trajs <- data.frame(matrix(nrow = length(all_phenos_combined) , ncol = n_points))
row.names(pheno_trajs) <- all_phenos_combined
colnames(pheno_trajs) <- seq(1,4, length.out = n_points)

plot_list <- list()
cnt <- 1
for (ph in c(all_hormones, all_phenos)){
  print(ph)
  res_gam <- gam_prot_tp_adj_covar(pheno, ph, subset(covariates, select = -c(storage_months,batch)), scale = T, predict = T)
  res_lmm <- lmm_prot_tp_poly3_adj_covar(pheno, ph, subset(covariates, select = -c(storage_months,batch)))
  gam_res_pheno_tp[cnt,] <- c(ph, unlist(res_gam[1:5]), res_lmm)
  cnt <- cnt + 1
  
  pheno_trajs[ph,] <- res_gam$predicted
  traj <- data.frame(TP = seq(1,4, length.out = length(res_gam$predicted)), pheno = res_gam$predicted, lower = res_gam$lower, upper = res_gam$upper)
  traj <- full_join(pheno_adjusted[,c("TP", "ID", ph)], traj, by = "TP")
  colnames(traj)[3] <- "values"
  traj$values <- scale(traj$values)
  plot_list[[ph]] <- ggplot(traj) + 
    geom_boxplot(aes(x = TP, y = values, group = TP), width = 0.1, color = my_colors[3], outliers = F) +
    geom_line(aes(x = TP, y = pheno), color = my_colors[4]) + 
    geom_ribbon(aes(x = TP, ymin = lower, ymax = upper), alpha = 0.2, fill = my_colors[4]) +
    geom_jitter(aes(x = TP, y = values), alpha = 0.2, width = 0.1, color = my_colors[3]) +
    labs(x = "Phase ", y = ph) +
    ggtitle(paste0(ph, "; GAM P = ", formatC(res_gam$pval, digits = 3))) + 
    theme_minimal() + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
}
pdf(paste0(out_basedir,"/plots/hormones_gam_withpoints2.pdf"), height = 7, width = 10)
grid.arrange(grobs = plot_list[names(plot_list) %in% all_hormones], ncol = 3, nrow = 2)  
dev.off()

pdf(paste0(out_basedir,"/plots/phenotypes_gam_withpoints2.pdf"), height = 10, width = 7)
grid.arrange(grobs = plot_list[names(plot_list) %in% all_phenos], ncol = 2, nrow = 5)  
dev.off()

gam_res_pheno_tp <- na.omit(gam_res_pheno_tp) %>%
  mutate(across(-c( pheno), as.numeric)) 

gam_res_pheno_tp$gam_BH_pval <- p.adjust(gam_res_pheno_tp$gam_pval, method = 'BH')
gam_res_pheno_tp$lmm_BH_pval <- p.adjust(gam_res_pheno_tp$lmm_pval, method = 'BH')

gam_res_pheno_tp <- gam_res_pheno_tp[order(gam_res_pheno_tp$gam_pval),]
gam_res_pheno_tp$gam_edf_round <- round(gam_res_pheno_tp$gam_edf)
gam_res_pheno_tp$gam_BH_sign <- ifelse(gam_res_pheno_tp$gam_BH_pval < 0.05,T,F)

write.table(gam_res_pheno_tp, file = paste0(out_basedir, "pheno_vs_tp_gam.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(pheno_trajs, file = paste0(out_basedir, "pheno_traj_gam.txt"), quote = F, sep = "\t", col.names = NA, row.names = T)


################################################################################
# GAM protein vs hormone/phenotype 
################################################################################

gam_res <- data.frame(matrix(nrow = length(all_prots) * length(all_hormones), ncol = 8))
colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

cnt <- 1
for (ph in all_hormones) {
  cat(ph, "\n")
  pb <- txtProgressBar(min = 1, max = length(all_prots), style = 3)
  i = 1
  for (prot in all_prots){
    res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'spline', anova_pval = F)
    res_lmm <- lmm_pheno_prot_adj_covar(d_wide, pheno, prot, ph, covariates, scale = T, adjust_timepoint = 'cubic')
    gam_res[cnt,] <- c(prot, ph, unlist(res_gam), res_lmm$pval)
    cnt <- cnt + 1
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval) %>%

cat("Number of BH significant associations:\n")
cat(nrow(gam_res[gam_res$BH_pval < 0.05,]), "\n")

write.table(gam_res, file = paste0(out_basedir, "prot_vs_allpheno_spline_gam.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)

gam_res <- read.delim(paste0(out_basedir,"prot_vs_allpheno_spline_gam.shared_prots.txt"), as.is = T, check.names = F, sep = "\t")
colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

gam_res <- gam_res[!gam_res$pheno %in% c("TSH", "FT4", "TST"),]

gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval) 
  
gam_res_hormones <- gam_res[gam_res$pheno %in% all_hormones,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)

gam_res_phenos <- gam_res[gam_res$pheno %in% all_phenos,] %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval)

nrow(gam_res_hormones[gam_res_hormones$BH_pval < 0.05,])
nrow(gam_res_phenos[gam_res_phenos$BH_pval < 0.05,])

write.table(gam_res_hormones, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(gam_res_phenos, file = paste0(out_basedir, "prot_vs_phenotypes.gam.spline.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)


### Plot Heatmaps

# signif in at least 1 hormone
prot_subs <- gam_res_hormones[gam_res_hormones$BH_pval < 0.05, "prot"]
h <- plot_association_heatmap(gam_res_hormones, prot_subs, transpose = T, cluster_cols = F)
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.BH0.05.transposed.pdf"), width = 4, height = 10, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

# signif in at least 2 hormones
prot_subs2 <- with(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, ], 
                   names(which(table(prot) >= 2)))
h <- plot_association_heatmap(gam_res_hormones, prot_subs2, transpose = T, cluster_cols = F)
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.BH0.05.transposed.min2horm.pdf"), width = 4, height = 5, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

# volcano
pdf(paste0(out_basedir, "plots/prot_vs_hormones.gam.spline.volcano.pdf"), width = 8, height = 6, useDingbats = F)
plot_association_volcano(gam_res_hormones)
dev.off()

### Phenotypes

# signif for at least 1 pheno
prot_subs <- gam_res_phenos[gam_res_phenos$BH_pval < 0.05, "prot"]
col_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
gam_res_phenos$pheno <- factor(gam_res_phenos$pheno, levels = col_order)
h <- plot_association_heatmap(gam_res_phenos, prot_subs, transpose = T, 
                              cluster_cols = F, col_order = col_order)
pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.BH0.05.transposed.pdf"), width = 6, height = 10, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

# signif for at least 2 pheno
prot_subs2 <- with(gam_res_phenos[gam_res_phenos$BH_pval < 0.05, ], 
                   names(which(table(prot) >= 2)))
col_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
gam_res_phenos$pheno <- factor(gam_res_phenos$pheno, levels = col_order)
h <- plot_association_heatmap(gam_res_phenos, prot_subs, transpose = T, 
                              cluster_cols = F, col_order = col_order)
pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.BH0.05.transposed.min2pheno.pdf"), width = 6, height = 10, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(h$gtable)
dev.off()

pdf(paste0(out_basedir, "plots/prot_vs_phenos.gam.spline.volcano.pdf"), width = 8, height = 8, useDingbats = F)
plot_association_volcano(gam_res_phenos)
dev.off()



hormone_assoc_count <- as.data.frame(table(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, "pheno"]))
pdf(paste0(out_basedir, "plots/hormone_assoc_barplot.pdf"), width = 4, height = 6, useDingbats = F)

ggplot(hormone_assoc_count, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Hormone", y = "Number of associated proteins") +
  theme_minimal()
dev.off()


# Upset
upset_matrix <- gam_res_hormones %>%
  filter(BH_pval < 0.05) %>%
  dplyr::select(prot, pheno) %>%
  distinct() %>%
  mutate(value = 1) %>%
  tidyr::pivot_wider(
    id_cols = prot,
    names_from = pheno,
    values_from = value,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("prot")

pdf(paste0(out_basedir, "plots/hormone_assoc_upset.pdf"), width = 6, height = 4, useDingbats = F)
# Create plot
upset(upset_matrix,
      nsets = ncol(upset_matrix),
      nintersects = 20,
      order.by = "freq",
      mainbar.y.label = "Number of Proteins",
      sets.x.label = "Proteins per Hormone")

dev.off()

################################################################################
# GAM hormones vs phenotypes all with all
################################################################################

gam_res_hormones_pheno <- data.frame(matrix(nrow = length(all_phenos_combined) * length(all_phenos_combined), ncol = 8))
colnames(gam_res_hormones_pheno) <- c("hormone", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")

tested_comb <- c()
cnt <- 1
for (hormone in all_phenos_combined) {
  cat(hormone, "\n")
  pb <- txtProgressBar(min = 1, max = length(all_phenos), style = 3)
  i = 1
  for (ph in all_phenos_combined){
    if ( (paste(hormone, ph) %in% tested_comb) ||  (paste(ph, hormone) %in% tested_comb)) next
    if (hormone == ph) next
    res_gam <- gam_prot_pheno_adj_covar(pheno, pheno, hormone, ph, subset(covariates, select = -c(storage_months,batch)), scale = T, adjust_timepoint = 'spline')
    res_lmm <- lmm_pheno_prot_adj_covar(pheno, pheno, hormone, ph, subset(covariates, select = -c(storage_months,batch)), scale = T, adjust_timepoint = 'cubic')
    
    gam_res_hormones_pheno[cnt,] <- c(hormone, ph, unlist(res_gam), res_lmm$pval)
    cnt <- cnt + 1
    i <- i + 1
    
    tested_comb <- c(tested_comb, paste(hormone, ph))
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

gam_res_hormones_pheno <- gam_res_hormones_pheno %>%
  na.omit(gam_res_hormones_pheno) %>%
  mutate(across(-c(pheno, hormone), as.numeric)) %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH')) %>%
  mutate(BH_lmm_pval = p.adjust(lmm_pval, method = 'BH')) %>%
  relocate(BH_pval, .after = pval) %>%
  relocate(BH_lmm_pval, .after = lmm_pval) %>%
  arrange(pval)

cat("Number of BH significant associations:\n")
cat(nrow(gam_res_hormones_pheno[gam_res_hormones_pheno$BH_pval < 0.05,]), "\n")

write.table(gam_res_hormones_pheno, file = paste0(out_basedir, "all_pheno_vs_all_pheno_spline_gam.shared_prots.txt"), quote = F, sep = "\t", row.names = FALSE)


################################################################################
# GAM hormones vs proteins adjust for PRS
################################################################################
prs <- read.delim("data/merged_protein_PRS.tsv", sep = "\t", as.is = T, check.names = F)
prs[grepl("^[0-9]",prs$IID), "IID"] <- paste0("X", prs[grepl("^[0-9]",prs$IID), "IID"])

prs <- prs[, colSums(is.na(prs)) <= 100]

gam_res_prs <- data.frame(matrix(nrow = length(all_prots) * length(all_phenos_combined), ncol = 7))
colnames(gam_res_prs) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples")

cnt <- 1
for (ph in all_phenos_combined) {
  cat(ph, "\n")
  
  pb <- txtProgressBar(min = 1, max = length(all_prots), style = 3)
  i = 1
  for (prot in all_prots){
    if (prot %in% colnames(prs)){
      cov_prs <- left_join(covariates, prs[,c("IID", prot)], by = c("ID" = "IID"))
      colnames(cov_prs)[ncol(cov_prs)] <- "PRS"
      res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, cov_prs, scale = T, adjust_timepoint = 'spline')
      
      gam_res_prs[cnt,] <- c(prot, ph, unlist(res_gam))
    } else {
      gam_res_prs[cnt,] <- c(prot, ph, rep(NA,5))
    }
    cnt <- cnt + 1
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

gam_res_prs <- gam_res_prs %>%
  na.omit(gam_res_prs) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval)

write.table(gam_res_prs, file = paste0(out_basedir, "prot_vs_allpheno_spline_gam.shared_prots.PRS.txt"), quote = F, sep = "\t", row.names = FALSE)

gam_res_prs = read.delim(paste0(out_basedir, "prot_vs_allpheno_spline_gam.shared_prots.PRS.txt"), sep = "\t", as.is = T, check.names = F)
gam_res = read.delim(paste0(out_basedir, "prot_vs_allpheno_spline_gam.shared_prots.txt"), sep = "\t", as.is = T, check.names = F)
colnames(gam_res) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples", "lmm_pval")
gam_res <- gam_res %>%
  na.omit(gam_res) %>%
  mutate(across(-c(pheno, prot), as.numeric)) %>% 
  arrange(pval) 

cmp <- left_join(gam_res[,c("prot", "pheno", "pval", "estimate", "n_samples")], gam_res_prs[,c("prot", "pheno", "pval", "estimate", "n_samples")], by = c("pheno", "prot"))

cmp_hormones <- cmp[cmp$pheno %in% all_hormones,]
cmp_pheno <- cmp[cmp$pheno %in% all_phenos,]

formatC(cor(cmp_hormones$estimate.x, cmp_hormones$estimate.y, use = 'complete.obs'), digits = 2)
formatC(cor(cmp_pheno$estimate.x, cmp_pheno$estimate.y, use = 'complete.obs'), digits = 2)

p1 <- ggplot(cmp[cmp$pheno %in% all_hormones,], aes(estimate.x, estimate.y)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Hormone - protein associations\ncomparison of estimates") + 
  theme_minimal()

p2 <- ggplot(cmp[cmp$pheno %in% all_hormones,], aes(-log10(pval.x), -log10(pval.y) )) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Hormone - protein associations\ncomparison of -log10(P)") +
  theme_minimal()

p3 <- ggplot(cmp[cmp$pheno %in% all_phenos,], aes(estimate.x, estimate.y)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Phenotype - protein associations\ncomparison of estimates") + 
  theme_minimal()

p4 <- ggplot(cmp[cmp$pheno %in% all_phenos,], aes(-log10(pval.x), -log10(pval.y) )) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = 'lm', color = 'dodgerblue4', size = 0.5) +
  xlab("no correction for PRS") + ylab ("with correction for PRS") + 
  ggtitle("Phenotype - protein associations\nComparison of -log10(P)") +
  theme_minimal()

pdf(paste0(out_basedir, "plots/PRS_comparison_associations.pdf"), height = 8, width = 8, useDingbats = F)
(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, face="bold"))
dev.off()



################################################################################
# Network
################################################################################



