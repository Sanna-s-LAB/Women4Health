my_colors <- c("#b71f57", "#96d1aa", "#099197", "#112f2c")
setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/")

library(ggplot2)
library(dplyr)
library(patchwork)
library(lubridate)
library(purrr)
library(tibble)

set.seed(123)
out_basedir <- "results12/"

#
# NB! Functions are at the end of the code, need to be run first
#


################################################################################
# Read and format the data
################################################################################

d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide.txt", as.is = T, check.names = F, sep = "\t", colClasses = c(ID = "character"))
batch_info <- read.delim("data/batch_info.txt", as.is = T, check.names = F, sep = "\t")

d_wide$TP <- gsub(".*_","", d_wide$SampleID)
d_wide$ID <- gsub("_.*","", d_wide$SampleID)
d_wide <- d_wide %>% select(SampleID, ID, TP, everything())

batch2_shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 50]
d_wide_shared <- d_wide[ ,batch2_shared_prots]
d_wide_b2 <- d_wide[d_wide$SampleID %in% batch_info[batch_info$Batch == 'batch2', "SampleID"],]

dim(d_wide)
dim(d_wide_b2)
dim(d_wide_shared)


################################################################################
# PCA on proteins with missing data (nipals) on shared proteins
################################################################################

res_pca <- run_pca_nipals_per_tp(d_wide_shared, nPCs = 70)
num_pcs_80_b12 <- res_pca$num_pcs_80
pca_per_tp_b12 <- res_pca$pca_per_tp

max(as.numeric(num_pcs_80_b12))
# [1] 60
write.table(pca_per_tp_b12, file = paste0(out_basedir, "olink_batch12_shared_prot_rm_outliers_4sd.PCA.txt"), quote = F, sep = "\t", row.names = FALSE)

################################################################################
# PCA on proteins with missing data (nipals) on batch2 data
################################################################################
res_pca_b2 <- run_pca_nipals_per_tp(d_wide_b2, nPCs = 60)

num_pcs_80_b2 <- res_pca_b2$num_pcs_80
pca_per_tp_b2 <- res_pca_b2$pca_per_tp

max(as.numeric(num_pcs_80_b2))
# [1] 49

write.table(pca_per_tp_b2, file = paste0(out_basedir, "olink_batch12_only_batch2_rm_outliers_4sd.PCA.txt"), quote = F, sep = "\t", row.names = FALSE)

################################################################################
# Technical covariates
################################################################################

#
# Season and storage time.
#

collect_date <- read.delim("../../phenotypes/batch12/date_of_collection_clean.txt", as.is = T, check.names = F, sep = "\t")
collect_date$date_collection <- as.Date(collect_date$date_collection, format = "%Y-%m-%d")
colnames(collect_date) <- gsub("SampleID", "ID", colnames(collect_date))
collect_date$SampleID <- paste0(collect_date$ID, "_", collect_date$visit_number)

batch_info <- read.delim("data/batch_info.txt", as.is = T, check.names = F, sep = "\t")
collect_date <- left_join(collect_date, batch_info, by = "SampleID")
colnames(collect_date) <- gsub("Batch", "batch", colnames(collect_date))

collect_date$shipment_date <- as.Date(ifelse(collect_date$batch == 'batch1', "24/09/2024", "01/09/2025"), format = "%d/%m/%Y")

collect_date$storage_months<- interval(collect_date$date_collection, collect_date$shipment_date) %/% months(1)
collect_date$storage_quarters <- round(collect_date$storage_months / 4)

colnames(collect_date) <- gsub("visit_number","TP",colnames(collect_date))
getSeason <- function(DATES) {
  WS <- as.Date("2025-12-21", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2025-3-20",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2025-6-21",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2025-9-23",  format = "%Y-%m-%d") # Fall Equinox
  
  d <- as.Date(strftime(DATES, format="2025-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Autumn")))
}
collect_date$season <- getSeason(collect_date$date_collection)

collect_date$date_collection <- NULL
collect_date$shipment_date <- NULL

write.table(collect_date, file = paste0(out_basedir, "correlations_with_covariates/season_storage_time.txt"), quote = F, sep = "\t", row.names = FALSE)
collect_date$season <- factor(collect_date$season, levels = c("Winter", "Spring", "Summer", "Autumn"))
collect_date$season_num <- as.numeric(collect_date$season)

season_kw <-  run_kruskal_test_each_TP(pca_per_tp_b12, collect_date[,c("ID", "TP","season")])
storage_lm <-  run_lm_each_TP(pca_per_tp_b12, collect_date[,c("ID", "TP","storage_months")])

write.table(season_kw, file = paste0(out_basedir, "correlations_with_covariates/shared_prots_season_KW.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(storage_lm, file = paste0(out_basedir, "correlations_with_covariates/shared_prots_storage_lm_per_tp.txt"), quote = F, sep = "\t", row.names = FALSE)

tmp <- left_join(pca_per_tp_b12, collect_date[,c("ID", "TP","season", "storage_months")], by = c("ID", "TP"))
ggplot(tmp, aes(y = PC5, x = PC6, color = season)) + 
  geom_point() + 
  scale_color_manual(values = my_colors) +
  theme_minimal()

pls <- plot_PCA_boxplot(pca_per_tp_b12, collect_date[,c("ID", "storage_months")], 'PC6', 'storage_months')
pdf(paste0(out_basedir, "correlations_with_covariates/plots/shared_prots_PC_vs_storage_months.pdf"), width = 10, height = 5)
pls[[1]] + pls[[2]]
dev.off()

pls <- plot_PCA_boxplot(pca_per_tp_b12, collect_date[,c("ID", "season", "storage_quarters")], 'PC1', 'season')
pdf(paste0(out_basedir, "correlations_with_covariates/plots/shared_prots_PC_vs_season.pdf"), width = 10, height = 5)
pls[[1]] + pls[[2]]
dev.off()

################################################################################
# Main covariates
################################################################################


pheno_0 <- read.delim("../../phenotypes/batch12/questionnaire_201025_selected_visit_0.txt", sep = "\t", check.names = F, as.is = T)
pheno_0[pheno_0 == 'NA'] <- NA
pheno_0$Visit_number <- NULL
pheno_0$Code <- NULL
pheno_0$ID <- gsub("X", "", pheno_0$ID )
pheno_0 <- pheno_0[pheno_0$ID %in% d_wide_shared$ID,]

factor_counts <- lapply(pheno_0, function(column) {
  if (length(unique(column)) < 5) {
    return(table(column, useNA = "ifany"))
  } else {
    return(NULL)
  }
})

for (col_name in names(factor_counts)) {
  cat("Factor levels for", col_name, ":\n")
  print(factor_counts[[col_name]])
  cat("\n")
}

continuous_pheno <- c("Age", "BMI", "Waist_circumference", "Hip_circumference","WHR", "numero_gravidanze","altezza_cm", "peso_kg", "anni_stop_pillola")


# Kruskal Wallis test on categorical covariates
covar_kw_res <- run_kruskal_test_each_TP(pca_per_tp_b12, pheno_0[,! colnames(pheno_0) %in% continuous_pheno])
colnames(covar_kw_res)[4] <- "pval"

# Linear regression on Age and BMI
pheno_cont_0 <- pheno_0[,c("ID", continuous_pheno)]
covar_lm_res <- run_lm_each_TP(pca_per_tp, pheno_cont_0)
colnames(season_kw)[4] <- "pval"

combined_covar_res <- rbind(covar_kw_res, subset(covar_lm_res, select = -c(beta)),
                            season_kw, subset(storage_lm, select = -c(beta)) )

num_tests_per_tp <- nrow(combined_covar_res[combined_covar_res$TP == '1',])
combined_covar_res$bonf_sign <- ifelse(combined_covar_res$pval < 0.05/num_tests_per_tp, T, F)

combined_covar_res$logp <- -log10(combined_covar_res$pval)
combined_covar_res$pval <- as.numeric(combined_covar_res$pval)

combined_covar_res <- combined_covar_res[order(combined_covar_res$pval),]

write.table(combined_covar_res, file = paste0(out_basedir, 'correlations_with_covariates/PC_vs_covariates_visit_0.txt'), quote = F, sep = "\t", row.names = FALSE)

col_order <- unique(combined_covar_res$covariate)
row_order <- paste("PC", seq(1,10), sep = "")

color_palette <- colorRampPalette(c("white", "#E6F2FF", "#99CCFF", "#3399FF", "#0066CC", "#004C99"))(100)
breaks <- seq(0, max(combined_covar_res$logp, na.rm = TRUE), length.out = 100)
color_fn <- colorRamp2(breaks, color_palette)

plot_list <- list()
for (tp in 1:4){
  logp_mat <- my_pivot_wider(combined_covar_res[combined_covar_res$TP == tp,c("PC", "covariate", "logp")], row_names = 'PC', names_from = 'covariate', values_from = 'logp')
  pval_mat <- my_pivot_wider(combined_covar_res[combined_covar_res$TP == tp,c("PC", "covariate", "pval")], row_names = 'PC', names_from = 'covariate', values_from = 'pval')
  logp_mat <- logp_mat[row_order, col_order]
  pval_mat <- pval_mat[row_order, col_order]
  labels_matrix <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  labels_matrix[pval_mat < 0.05/num_tests_per_tp] <- "*"
  
  plot_list[[tp]] <- Heatmap(
    logp_mat,
    name = "-log10(p)",
    border = TRUE,
    rect_gp = gpar(col = "grey", lwd = 1),
    col = color_fn,
    cluster_rows = F,
    cluster_columns = F,
    row_names_gp = gpar(fontsize = 14),  # Row font size
    column_names_gp = gpar(fontsize = 14),  # Column font size
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(labels_matrix[i, j] != "") {
        grid.text(labels_matrix[i, j], x, y, gp = gpar(fontsize = 12))
      }
    },
    column_title = paste("Visit", tp),
    show_heatmap_legend = if(tp == 1) TRUE else FALSE  
  )
}
pdf(paste0(out_basedir, 'correlations_with_covariates/plots/corrplots_PC_vs_covariates.pdf'), width = 20, height  = 5)
draw(plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]], gap = unit(1, "cm"))
dev.off()


################################################################################
# Functions
################################################################################


run_kruskal_test_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  kw_res <- data.frame(matrix(ncol = 4))
  colnames(kw_res) <- c("TP", "PC", "covariate", "KW_test_pval")
  cnt <- 1
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]

    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      merged <- left_join(pca10, tmp_covariates, by = c("ID", "TP"))
    } else {
      tmp_covariates <- covariates
      merged <- left_join(pca10, tmp_covariates, by = "ID")
    }
    
    all_covs <- colnames(tmp_covariates)[! colnames(tmp_covariates) %in% c("SampleID", "ID", "TP")]
    all_pcs <- colnames(pca10)[! colnames(pca10) %in% c("SampleID", "ID", "TP")]
    for (pc in all_pcs){
      for (cov in all_covs){
        subs <- merged[,c(pc,cov) ]
        colnames(subs) <- c("PC", "cov")
        kw <- kruskal.test(PC ~ cov, data=subs)$p.value
        
        kw_res[cnt,] <- c(tp, pc, cov, kw)
        cnt <- cnt + 1
      }
    }
  }
  
  kw_res <- na.omit(kw_res) %>%
    mutate(across(-c(PC, covariate), as.numeric))   
  
  kw_res
}

run_lm_each_TP <- function(pca_per_tp, covariates, num_pcs = 10){
  lm_res <- data.frame(matrix(ncol = 5))
  colnames(lm_res) <- c("TP", "PC", "covariate", "beta", "pval")
  cnt <- 1
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      merged <- left_join(pca10, tmp_covariates, pca10, by = c("ID", "TP"))
    } else {
      tmp_covariates <- covariates
      merged <- left_join(pca10, tmp_covariates, pca10, by = "ID")
    }
    
    all_covs <- colnames(tmp_covariates)[! colnames(tmp_covariates) %in% c("SampleID", "ID", "TP")]
    all_pcs <- colnames(pca10)[! colnames(pca10) %in% c("SampleID", "ID", "TP")]
    
    for (pc in all_pcs){
      for (cov in all_covs){
        subs <- merged[,c(pc,cov) ]
        colnames(subs) <- c("PC", "cov")
        lm_fit <- lm(PC ~ cov, data = subs)
        lm_res[cnt,] <- c(tp, pc, cov, summary(lm_fit)$coefficients['cov', 1], summary(lm_fit)$coefficients['cov', 4])
        cnt <- cnt + 1
      }
    }
  }
  
  lm_res <- na.omit(lm_res) %>%
    mutate(across(-c(PC, covariate), as.numeric))   
  
  lm_res
}

plot_PCA_boxplot <- function(pca_per_tp, covariates, pc, cov){
  pca_all_tps <- data.frame()
  for (tp in 1:4){
    pca10 <- pca_per_tp[pca_per_tp$TP == tp,]
    pca10$TP <- NULL
    
    if('TP' %in% colnames(covariates)){
      tmp_covariates <- covariates[covariates$TP == tp,]
      tmp_covariates$TP = NULL
    } else {
      tmp_covariates <- covariates
    }
    
    merged <- inner_join(tmp_covariates, pca10, by = c("ID"))
    pca_all_tps <- rbind(pca_all_tps, cbind(tp, merged))
  }
  
  subs <- pca_all_tps[c("PC1", "PC2", pc, cov, "tp")]
  colnames(subs) <- c("PC1", "PC2", "PC", "cov", "TP")
  if (length(unique(subs$cov)) < 5){
    p1 <- ggplot(subs, aes(x = cov, y = PC, group = cov)) + geom_boxplot() + theme_bw() + xlab(cov) + ylab(pc) + facet_wrap(~TP)
    p2 <- ggplot(pca_all_tps, aes(x = PC1, y = PC2, color = season)) + geom_point() + scale_color_manual(values = my_colors) + theme_bw()
  } else {
    p1 <- ggplot(subs, aes(x = cov, y = PC)) + geom_point() + geom_smooth(method = 'lm')+ theme_bw() + xlab(cov) + ylab(pc) + facet_wrap(~TP)
    p2 <- ggplot(subs, aes(x = PC1, y = PC2, color = cov)) + geom_point() + theme_bw() + facet_wrap(~TP)
  }
  return(list(p1,p2))
}

run_pca_nipals_per_tp <- function(d_wide, nPCs = 70){
  num_pcs_80 <- c()
  pca_per_tp <- data.frame()
  for (tp in 1:4){
    print(tp)
    tmp_wide <- as.data.frame(subset(d_wide[d_wide$TP == tp,], select = -c(TP, SampleID, ID)))
    row.names(tmp_wide) <- d_wide[d_wide$TP == tp, "ID"]
    pca <- pcaMethods::pca(tmp_wide, method = 'nipals', nPcs = nPCs, center = T, scale = 'vector')
    cumulative_variance <- cumsum(pca@R2)
    cat("10 PCs explain", cumulative_variance[10], " of variance\n")
    num_pcs_80 <- c(num_pcs_80, which(cumulative_variance >= 0.80)[1])
    
    pca10 <- as.data.frame(pca@scores)[,1:10] %>%
      rownames_to_column(var = 'ID') 
    pca_per_tp <- rbind(pca_per_tp, data.frame(TP = tp, pca10))
  }
  pca_per_tp$SampleID <- paste0(pca_per_tp$ID, "_", pca_per_tp$TP)
  pca_per_tp <- pca_per_tp %>% select(SampleID, ID, TP, everything())
  return(list(pca_per_tp = pca_per_tp, num_pcs_80 = num_pcs_80))
}

