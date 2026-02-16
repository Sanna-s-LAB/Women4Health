library(reshape2)
library(rlang)
library(parallel)
library(mgcv)

################################################################################
# Read the data
################################################################################

all_phases <- c("F", "O", "EL", "LL")

# Read protein data
d_wide <- read.delim("data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.phase_avg.txt", as.is = T, check.names = F, sep = "\t")

if (! "ID" %in% colnames(d_wide)){
  d_wide$ID <- gsub("_.*", "", d_wide$SampleID)
  d_wide$phase <- gsub(".*_", "", d_wide$SampleID)
  d_wide <- d_wide %>%
    dplyr::select(SampleID, ID, phase, everything())
}

d_wide$phase <- relevel(factor(d_wide$phase, levels = c("F", "O", "EL", "LL")), ref = "F")
# Make a df with proteins available for both batch 1 and 2
shared_prots <- colnames(d_wide)[colSums(is.na(d_wide)) < 50]
d_wide <- d_wide[ ,shared_prots]

# Read covariate data
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
covariate_names <- c("Age","BMI","batch", "from", "storage_months")
covariates$from <- relevel(as.factor(covariates$from), ref = "X")

# Read hormone + phenotype data:
pheno <- read.delim("../../phenotypes/batch12/cleaned_phenotypes_251125_uniformed_adjusted.withHOMA.log_some.phase_avg.txt", as.is = T, check.names = F, sep = "\t")            

all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
all_phenos <- c("INS", "HOMA_B", "HOMA_IR", "GL", "AST", "ALT", "TRI", "COL", "HDL", "LDL")

if (! "ID" %in% colnames(pheno)){
  pheno$ID <- gsub("_.*", "", pheno$SampleID)
  pheno$phase <- gsub(".*_", "", pheno$SampleID)
  pheno <- pheno %>%
    dplyr::select(SampleID, ID, phase, everything())
}

pheno$phase <- relevel(factor(pheno$phase, levels = c("F", "O", "EL", "LL")), ref = "F")

################################################################################
# Function to run the cross-lagged temporal causality model
################################################################################

run_wide_models <- function(d_wide, pheno, covariates, pair){
  protein <- pair$prot
  ph <- pair$pheno
  
  subs <- left_join(
    left_join(d_wide[,c("SampleID", "ID", "phase", protein)], pheno[,c("SampleID", ph)], by = "SampleID"), 
    covariates, by = c("SampleID", "ID", "phase"))
  
  colnames(subs)[c(4,5)] <- c("prot", "pheno")
  subs$SampleID <- NULL
  subs$pheno <- scale(subs$pheno)
  subs$prot <- scale(subs$prot)
  
  wide_data <- reshape(
    subs,
    timevar = "phase",
    idvar = "ID",
    direction = "wide"
  )
  wide_data$ID <- as.factor(wide_data$ID)
  
  # 2 -> 1
  formula_direct_1 <- as.formula("prot.O ~ pheno.F + prot.F + 
                                 Age.F + BMI.F +batch.F + from.F + storage_months.F + s(ID, bs='re')")
  formula_inverse_1 <- as.formula("pheno.O ~ prot.F + pheno.F + 
                                  Age.F + BMI.F +batch.F + from.F + storage_months.F + s(ID, bs='re')")
  
  fit_direct_1 <- gam(formula_direct_1, method = 'REML', data = wide_data)
  fit_inverse_1 <- gam(formula_inverse_1, method = 'REML', data = wide_data)
  
  p_direct_1 <- summary(fit_direct_1)$p.table["pheno.F", 4]
  p_inverse_1 <- summary(fit_inverse_1)$p.table["prot.F", 4]
  
  est_direct_1 <- formatC(summary(fit_direct_1)$p.table["pheno.F", 1], digits = 2)
  est_inverse_1 <- formatC(summary(fit_inverse_1)$p.table["prot.F", 1], digits = 2)
    
  direction1 <- ifelse(p_direct_1 < 0.05 & p_inverse_1 >= 0.05, paste(ph, "→", protein),
                         ifelse(p_inverse_1 < 0.05 & p_direct_1 >= 0.05, paste(protein, "→", ph),
                                ifelse(p_direct_1 < 0.05 & p_inverse_1 < 0.05, "both directions significant", "nothing significant")))
  
  estimate1 <- ifelse(p_direct_1 < 0.05 & p_inverse_1 >= 0.05, est_direct_1,
                       ifelse(p_inverse_1 < 0.05 & p_direct_1 >= 0.05, est_inverse_1,
                              ifelse(p_direct_1 < 0.05 & p_inverse_1 < 0.05, paste(est_direct_1, est_inverse_1), "NA")))
  
  
  # 3 -> 2
  formula_direct_2 <- as.formula("prot.EL ~ pheno.O + prot.O + 
                                 Age.O + BMI.O +batch.O + from.O + storage_months.O + s(ID, bs='re')")
  formula_inverse_2 <- as.formula("pheno.EL ~ prot.O + pheno.O +
                                  Age.O + BMI.O +batch.O + from.O + storage_months.O + s(ID, bs='re')")
  
  fit_direct_2 <- gam(formula_direct_2, method = 'REML', data = wide_data)
  fit_inverse_2 <- gam(formula_inverse_2, method = 'REML', data = wide_data)
  
  p_direct_2 <- summary(fit_direct_2)$p.table["pheno.O", 4]
  p_inverse_2 <- summary(fit_inverse_2)$p.table["prot.O", 4]
  
  est_direct_2 <- formatC(summary(fit_direct_2)$p.table["pheno.O", 1], digits = 2)
  est_inverse_2 <- formatC(summary(fit_inverse_2)$p.table["prot.O", 1], digits = 2)
  
  direction2 <- ifelse(p_direct_2 < 0.05 & p_inverse_2 >= 0.05, paste(ph, "→", protein),
                       ifelse(p_inverse_2 < 0.05 & p_direct_2 >= 0.05, paste(protein, "→", ph),
                              ifelse(p_direct_2 < 0.05 & p_inverse_2 < 0.05, "both directions significant", "nothing significant")))
  estimate2 <- ifelse(p_direct_2 < 0.05 & p_inverse_2 >= 0.05, est_direct_2,
                      ifelse(p_inverse_2 < 0.05 & p_direct_2 >= 0.05, est_inverse_2,
                             ifelse(p_direct_2 < 0.05 & p_inverse_2 < 0.05, paste(est_direct_2, est_inverse_2), "NA")))
  
    
  # 4 -> 3
  formula_direct_3 <- as.formula("prot.LL ~ pheno.EL + prot.EL  +
                                 Age.EL + BMI.EL +batch.EL + from.EL + storage_months.EL + s(ID, bs='re')")
  formula_inverse_3 <- as.formula("pheno.LL ~ prot.EL + pheno.EL  +
                                  Age.EL + BMI.EL +batch.EL + from.EL + storage_months.EL + s(ID, bs='re')")
  
  fit_direct_3 <- gam(formula_direct_3, method = 'REML', data = wide_data)
  fit_inverse_3 <- gam(formula_inverse_3, method = 'REML', data = wide_data)
  
  p_direct_3 <- summary(fit_direct_3)$p.table["pheno.EL", 4]
  p_inverse_3 <- summary(fit_inverse_3)$p.table["prot.EL", 4]
  
  est_direct_3 <- formatC(summary(fit_direct_3)$p.table["pheno.EL", 1], digits = 2)
  est_inverse_3 <- formatC(summary(fit_inverse_3)$p.table["prot.EL", 1], digits = 2)
  
  direction3 <- ifelse(p_direct_3 < 0.05 & p_inverse_3 >= 0.05, paste(ph, "→", protein),
                       ifelse(p_inverse_3 < 0.05 & p_direct_3 >= 0.05, paste(protein, "→", ph),
                              ifelse(p_direct_3 < 0.05 & p_inverse_3 < 0.05, "both directions significant", "nothing significant")))
  estimate3 <- ifelse(p_direct_3 < 0.05 & p_inverse_3 >= 0.05, est_direct_3,
                      ifelse(p_inverse_3 < 0.05 & p_direct_3 >= 0.05, est_inverse_3,
                             ifelse(p_direct_3 < 0.05 & p_inverse_3 < 0.05, paste(est_direct_3, est_inverse_3), "NA")))
  
  
  l<- list (
    data.frame(
      phenotype = ph,
      protein = protein,
      summary = paste0("1:", direction1,";2:", direction2, ";3:", direction3),
      summary_short = gsub(";?[1-3]:nothing significant","",paste0("1:", direction1,";2:", direction2, ";3:", direction3) ),
      estimate = gsub(";?t[1-3]:NA;?","", paste0("t1:", estimate1,";t2:", estimate2, ";t3:", estimate3) )
    ),
    data.frame(
      direction_1_to_2 = direction1,
      p_sem_X1_to_Y2 = p_direct_1,
      p_sem_Y1_to_X2 = p_inverse_1,
      est_sem_X1_to_Y2 = est_direct_1,
      est_sem_Y1_to_X2 = est_inverse_1
    ),
    data.frame(
      direction_2_to_3 = direction2,
      p_sem_X2_to_Y3 = p_direct_2,
      p_sem_Y2_to_X3 = p_inverse_2,
      est_sem_X2_to_Y3 = est_direct_2,
      est_sem_Y2_to_X3 = est_inverse_2
    ),
    data.frame(
      direction_3_to_4 = direction3,
      p_sem_X3_to_Y4 = p_direct_3,
      p_sem_Y3_to_X4 = p_inverse_3,
      est_sem_X3_to_Y4 = est_direct_3,
      est_sem_Y3_to_X4 = est_inverse_3
    )
  )
  return(do.call(cbind,l))
}

################################################################################
# Run the cross-lagged temporal causality model on significant associations
################################################################################

# proteins vs phenotypes/hormones

gam_res <- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.shared_prots.txt"), as.is = T, check.names = F, sep = "\t")
pairs <- gam_res[gam_res$BH_pval < 0.05, c("prot", "pheno")]

ncores <- detectCores()-1
results_list <- mclapply(seq_len(nrow(pairs)), function(i) {
  run_wide_models(d_wide, pheno, covariates, pairs[i, ])
}, mc.cores = ncores)
results_df <- do.call(rbind, results_list)

write.table(results_df, file = paste0(out_basedir, "prot_vs_phenotypes.causality_wide.linear.txt"), quote = F, sep = "\t", row.names = F)

# phenotypes vs hormones
pairs <- as.data.frame(t(combn(all_phenos_combined, 2)))
colnames(pairs) <- c("prot", "pheno")

ncores <- detectCores()-1
results_list <- mclapply(seq_len(nrow(pairs)), function(i) {
  run_wide_models(pheno, pheno, covariates, pairs[i, ])
}, mc.cores = ncores)

results_df <- do.call(rbind, results_list)

write.table(results_df, file = paste0(out_basedir, "hormone_vs_phenotypes.causality_wide.linear.txt"), quote = F, sep = "\t", row.names = F)


