
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")
#my_colors <- c("#eda048", "#4d9bb0", "#075c62", "#86e6ca")
#my_colors <- c("#f5b026", "#e26413", "#9a4020", "#63552d")
#my_colors <- c("#097054", "#FFDE00", "#6599FF", "#FF9900")

#library(tidyverse)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyverse)
library(mgcv)
library(ggtext)
library(RColorBrewer)
library(limma)



#' Performs association analysis between protein levels and phase or visit using GAMs
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows). 
#' @param prot protein name to run the GAM for 
#' @param covariates data frame with all covariates to add to the model
#' @param scale Logical. Whether to scale the data. Default is FALSE.
#' @param rm_outliers Logical. Whether to remove outliers. Default is FALSE.
#' @param predict Logical. Whether to generate predicted fitted values. Default is TRUE.
#' @param anova_pval Logical. Whether to compute ANOVA p-values instead of normal GAM reported. Default is FALSE.
#' @param n_points Numeric. Number of time points to use for predictions. Default is 20.
#' 
gam_prot_tp_adj_covar <- function(d_wide, prot, covariates, scale = F, rm_outliers = F, predict = T, anova_pval = F, n_points = 20){
  # if d_wide has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase = NULL
  }
  
  # combine protein and covariate datasets
  d_subs <- inner_join(d_wide[,c(prot, "SampleID", "ID", "TP")], covariates, by = c("SampleID", "ID", "TP"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$ID <- as.factor(d_subs$ID)
  d_subs <- na.omit(d_subs)
  
  if (rm_outliers) d_subs <- remove_outliers_zscore(d_subs, "prot")
  
  if (scale) d_subs$prot <- scale(d_subs$prot)
  
  covariate_names = colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  if ("batch" %in% colnames(d_subs)){
    if (table(d_subs$batch)["batch1"] == 0) covariate_names = covariate_names[covariate_names != "batch"]
  }
  
  # make GAM formula
  fo_gam <- as.formula(paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+")))
  fo_gam_null <- as.formula(paste("prot ~ s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+")))
  
  # Run the model
  model <- gam(fo_gam, data = d_subs,  method = 'REML')
  
  if (anova_pval){
    model0 <- gam(fo_gam_null, data = d_subs, method = 'REML')
    an <- anova.gam(model, model0)
    pval <- an$`Pr(>F)`[2]
  } else {
    pval <- summary(model)$s.table["s(TP)","p-value"]
  }
  edf <- summary(model)$s.table["s(TP)","edf"]
  fval <- summary(model)$s.table["s(TP)","F"]
  
  if (predict){
    covar_means <- as.data.frame(lapply(covariates[,covariate_names], function(x) {
      if(is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else {
        levels(x)[1]  # Use first factor level
      }
    }))
    
    new_data <- expand.grid(
      TP = seq(1, 4, length.out = n_points),
      ID = unique(d_subs$ID),
      predicted = NA
    ) %>%
      bind_cols(
        covar_means[1,] 
      )
    
    predictions <- predict.gam(model, newdata = new_data,  exclude = "s(ID)", se.fit = T)
    new_data$predicted <- predictions$fit
    new_data$SE <- predictions$se.fit
    new_data$lower <- new_data$predicted - 1.96 * new_data$SE
    new_data$upper <- new_data$predicted + 1.96 * new_data$SE
    
    new_data2 <- unique(new_data[,c("TP", "predicted", "lower", "upper")])
    
    return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), predicted = new_data2$predicted, lower = new_data2$lower, upper = new_data2$upper))
  } 
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

#' Performs association analysis between protein and hormone/phenotype levels
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows). 
#' @param pheno data frame with phenotypes (in columns) for all samples (in rows). 
#' @param prot protein name to use in the association
#' @param ph phenotype name to use in the association
#' @param covariates data frame with all covariates to add to the model
#' @param scale Logical. Whether to scale the data. Default is FALSE.
#' @param rm_outliers Logical. Whether to remove outliers. Default is FALSE.
#' @param adjust_timepoint how to adjust for the phase/visit. Can be one of "none" (do not adjust for timepoint), "linear" (add timepoint as a parameteric term), "spline" (add timepoint as a spline term)
#' @param adjust_pheno how to add the phenotype to the model: "linear" - as a linear parameteric term (default) or "spline" - as a spline term
#' @param longitudinal Logical. Whether to run a GAM with random intercept (default) or a simple lm with no random effect
#' @param add_age_interaction Logical. Whether to add interaction with age to the model
#' 
gam_prot_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, rm_outliers = F, adjust_timepoint = 'spline', adjust_pheno = 'linear', anova_pval = F, predict = F, add_age_interaction = F, longitudinal = T){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
    
    #covariates$TP <- as.numeric(covariates$phase)
    #covariates$phase <- NULL
  }
  
  if(! ("TP" %in% colnames(covariates) || "phase" %in% colnames(covariates)) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID", "TP"))
    d_subs$SampleID.y <- NULL
    covariates$TP = NULL
  }
  
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  if ("batch" %in% colnames(d_subs)){
    if (table(na.omit(d_subs)$batch)["batch1"] == 0) covariate_names = covariate_names[covariate_names != "batch"]
  }
  # generate the GAM formula
  if (longitudinal){
    if (adjust_timepoint == 'spline'){
      fo_gam <- paste("prot ~ s(pheno) + s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
      fo_gam_null <- paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
    } else if (adjust_timepoint == 'linear') {
      fo_gam <- paste("prot ~ s(pheno) + TP + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
      fo_gam_null <- paste("prot ~ TP + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
    } else if (adjust_timepoint == 'none') {
      fo_gam <- paste("prot ~ s(pheno) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
      fo_gam_null <- paste("prot ~  s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+"))
    } else {
      stop ("Wrong adjust_timepoint argument. Should be one of spline, linear or none.")
    }
  } else {
    fo_gam <- paste("prot ~ s(pheno) + ", paste(covariate_names, collapse = "+"))
    fo_gam_null <- paste("prot ~ ", paste(covariate_names, collapse = "+"))
  }
  
  
  if (add_age_interaction) {
    if (! "Age" %in% colnames(d_subs)) {cat ("No Age covariate provided for the interaction!\n")}
    fo_gam <- paste0(fo_gam, " + pheno * Age")
    d_subs$Age <- scale(d_subs$Age)
  }

  # Linear relation between protein and phenotype
  if (adjust_pheno != 'spline'){
    fo_gam <- gsub("s\\(pheno\\)", "pheno", fo_gam)
    
    model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
    
    est <- summary(model)$p.table["pheno","Estimate"]
    se <- summary(model)$p.table["pheno","Std. Error"]
    pval <- summary(model)$p.table["pheno","Pr(>|t|)"]
    
    if (add_age_interaction) {
      interaction_pval <- summary(model)$p.table["pheno:Age","Pr(>|t|)"]
      return(list(pval = pval,  est = est, se = se, n = nrow(d_subs), n_samples = length(unique(d_subs$ID)), age_inter_pval = interaction_pval))
    } 
    return(list(pval = pval,  est = est, se = se, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
  }
  
  # NON-linear relation between protein and phenotype
  model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
  
  edf <- round(summary(model)$s.table["s(pheno)","edf"])
  fval <- summary(model)$s.table["s(pheno)","F"]
  
  if (anova_pval){
    model0 <- gam(as.formula(fo_gam_null), data = d_subs, method = 'REML')
    an <- anova.gam(model, model0)
    pval <- an$`Pr(>F)`[2]
  } else {
    pval <- summary(model)$s.table["s(pheno)","p-value"]
  }
  
  return(list(pval = pval,  edf = edf, fval = fval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

#' test for association of a protein vs all hormones together
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows)
#' @param pheno data frame with hormones in columns and samples in rows
#' @param covariates data frame of covariates to add to the model
#' @param adjust_timepoint how to adjust for the phase/visit. Can be one of "none" (do not adjust for timepoint), "linear" (add timepoint as a parameteric term), "spline" (add timepoint as a spline term)
#' 
gam_prot_all_pheno_together_adj_covar <- function(d_wide, pheno, prot, covariates, scale = F, adjust_timepoint = 'spline'){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
  }
  
  colnames(pheno) <- gsub("17BES", "X17BES",colnames(pheno))
  if(! "TP" %in% colnames(covariates) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno, by = c("SampleID", "TP", "ID")),
                         covariates, by = c("SampleID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno, by = c("SampleID", "TP", "ID")),
                         covariates, by = c("SampleID","ID", "TP"))
    covariates$TP = NULL
  }
  colnames(d_subs)[1:4] <- c("SampleID", "ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  d_subs$ID <- as.factor(d_subs$ID)
  
  #if (scale) {
  #  d_subs$prot <- scale(d_subs$prot)
  #  d_subs$pheno <- scale(d_subs$pheno)
  #}

  pheno_names <- colnames(pheno)[!colnames(pheno ) %in% c("ID", "SampleID", "TP")]
  if (adjust_timepoint == 'spline'){
    fo_gam <- paste("prot ~  s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
    fo_gam_null <- paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
  } else if (adjust_timepoint == 'linear') {
    fo_gam <- paste("prot ~  TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
    fo_gam_null <- paste("prot ~ TP + s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
  } else if (adjust_timepoint == 'none') {
    fo_gam <- paste("prot ~  s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
    fo_gam_null <- paste("prot ~  s(ID,  bs = 're') + ", paste(colnames(covariates)[-1], collapse = "+"), "+", paste(pheno_names, collapse = "+"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of spline, linear or none.")
  }
  
  # Run the GAM
  model <- gam(as.formula(fo_gam), data = d_subs, method = 'REML')
  
  ests <- summary(model)$p.table[pheno_names,"Estimate"]
  ses <- summary(model)$p.table[pheno_names,"Std. Error"]
  pvals <- summary(model)$p.table[pheno_names,"Pr(>|t|)"]
  
  return(list(pvals = pvals,  ests = ests, ses = ses, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

#' Performs association analysis between protein levels and phase or visit using LMMs
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows). 
#' @param prot protein name to run the GAM for 
#' @param covariates data frame with all covariates to add to the model
#' @param scale Logical. Whether to scale the data. Default is FALSE.
#' 
lmm_prot_tp_poly3_adj_covar <- function(d_wide, prot, covariates, scale = F){
  # if phases not visits, convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase = NULL
  }
  d_subs <- inner_join(d_wide[,c(prot, "SampleID","ID", "TP")], covariates, by = c("SampleID", "ID", "TP"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  covariate_names = colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
  fo_lmm <- as.formula(paste("prot ~ poly(TP,3) +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  model <- lmer(fo_lmm, data = d_subs)
  fo_lmm_base <- as.formula(paste("prot ~ ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  model0 <- lmer(fo_lmm_base, data = d_subs)
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return(pval)
}

#' Performs association analysis between protein levels and phase or visit using LMMs treating the phase/visit as a factor
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows). 
#' @param prot protein name to run the LMM for 
#' @param covariates data frame with all covariates to add to the model
#' @param scale Logical. Whether to scale the data. Default is TRUE
#' 
lmm_prot_tp_factor_adj_covar <- function(d_wide, prot, covariates, scale = T){
  # if phases not visits, convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase = NULL
  }
  d_subs <- inner_join(d_wide[,c(prot, "SampleID","ID", "TP")], covariates, by = c("SampleID", "ID", "TP"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) d_subs$prot <- scale(d_subs$prot)
  
  covariate_names = colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
  fo_lmm <- as.formula(paste("prot ~ TP +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  model <- lmer(fo_lmm, data = d_subs)
  coefs <- summary(model)$coefficients
  
  return(list(betas = coefs[c("TP2", "TP3", "TP4"), "Estimate"], pvals = coefs[c("TP2", "TP3", "TP4"), "Pr(>|t|)"]))
}

#' Performs association analysis between protein and hormone/phenotype levels using LMMs
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows). 
#' @param pheno data frame with phenotypes (in columns) for all samples (in rows). 
#' @param prot protein name to use in the association
#' @param ph phenotype name to use in the association
#' @param covariates data frame with all covariates to add to the model
#' @param scale Logical. Whether to scale the data. Default is FALSE.
#' @param adjust_timepoint how to adjust for the phase/visit. Can be one of "none" (do not adjust for timepoint), "linear" (add timepoint as a linear term), "cubic" (add timepoint as 3rd degree polynomial terms)
#' @param longitudinal Logical. Whether to run a LMM with random intercept (default) or a simple lm with no random effect
#' 
lmm_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic", longitudinal = T){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase <- NULL
  }
  
  if(! ("TP" %in% colnames(covariates) || "phase" %in% colnames(covariates)) ){
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID"))
  } else {
    d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                         covariates, by = c("ID", "TP"))
    d_subs$SampleID.y <- NULL
    covariates$TP = NULL
  }
  
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  if ("batch" %in% colnames(subs)){
    if (table(subs$batch)["batch1"] == 0) covariate_names = covariate_names[covariate_names != "batch"]
  }
  if (longitudinal){
    if (adjust_timepoint == 'cubic'){
      fo_lmm <- as.formula(paste("prot ~ poly(TP, 3) + pheno +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
      fo_lmm_base <- as.formula(paste("prot ~ poly(TP, 3) + ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
    } else if (adjust_timepoint == 'linear') {
      fo_lmm <- as.formula(paste("prot ~ TP + pheno +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
      fo_lmm_base <- as.formula(paste("prot ~ TP + ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
    } else if (adjust_timepoint == 'none') {
      fo_lmm <- as.formula(paste("prot ~ pheno +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
      fo_lmm_base <- as.formula(paste("prot ~ ", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
    } else {
      stop ("Wrong adjust_timepoint argument. Should be one of cubic, linear or none.")
    }
    model <- lmer(fo_lmm, data = d_subs)

  } else {
    fo_lmm <- as.formula(paste("prot ~ pheno +", paste(covariate_names, collapse = "+")))
    fo_lmm_base <- as.formula(paste("prot ~ ", paste(covariate_names, collapse = "+")))
    model <- lm(fo_lmm, data = d_subs)
  }
  est <- summary(model)$coefficients["pheno", "Estimate"]
  se <- summary(model)$coefficients["pheno","Std. Error"]
  tval <- summary(model)$coefficients["pheno","t value"]
  pval <- summary(model)$coefficients["pheno","Pr(>|t|)"]
  
  return(list(estimate = est, pval = pval, se = se, tval = tval, n = nrow(d_subs), n_samples = length(unique(d_subs$ID))))
}

lmm_prot_tp_interaction_pheno_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F, adjust_timepoint = "cubic"){
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("ID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  if (adjust_timepoint == 'cubic'){
    fo_lmm <- as.formula(paste("prot ~ poly(TP, 3) + pheno +  poly(TP, 3) * pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ poly(TP, 3) + pheno +", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else if (adjust_timepoint == 'linear') {
    fo_lmm <- as.formula(paste("prot ~ TP + pheno + TP * pheno + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
    fo_lmm_base <- as.formula(paste("prot ~ TP + pheno + ", paste(colnames(covariates)[-1], collapse = "+"), "+ (1|ID)"))
  } else {
    stop ("Wrong adjust_timepoint argument. Should be one of cubic or linear.")
  }
  
  model <- lmer(fo_lmm, data = d_subs)
  model0 <- lmer(fo_lmm_base, data = d_subs)
  
  #est <- summary(model)$coefficients["pheno", "Estimate"]
  #se <- summary(model)$coefficients["pheno",2]
  #tval <- summary(model)$coefficients["pheno",3]
  an <- suppressMessages(anova(model, model0))
  pval <- an$`Pr(>Chisq)`[2]
  
  return( pval)
}

# compare GAM and LMM
compare_gam_lmm <- function(d_wide, prot, covariates){
  # if d_wide has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase = NULL
  }
  
  # combine protein and covariate datasets
  d_subs <- inner_join(d_wide[,c(prot, "SampleID", "ID", "TP")], covariates, by = c("SampleID", "ID", "TP"))
  colnames(d_subs)[1] <- "prot"
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$ID <- as.factor(d_subs$ID)
  d_subs <- na.omit(d_subs)
  
  d_subs$prot <- scale(d_subs$prot)
  
  covariate_names = colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  if ("batch" %in% colnames(d_subs)){
    if (table(d_subs$batch)["batch1"] == 0) covariate_names = covariate_names[covariate_names != "batch"]
  }
  
  # make GAM formula
  fo_gam <- as.formula(paste("prot ~ s(TP, k = 4) + s(ID,  bs = 're') + ", paste(covariate_names, collapse = "+")))
  fo_lmm <- as.formula(paste("prot ~ poly(TP,3) +", paste(covariate_names, collapse = "+"), "+ (1|ID)"))
  # run the models
  model_gam <- gam(fo_gam, data = d_subs,  method = 'ML')
  model_lmm <- lmer(fo_lmm, data = d_subs, REML = FALSE)
  
  aic <-AIC(model_lmm, model_gam)
  bic <- BIC(model_lmm, model_gam)
   
  return(list(aic,bic))
}



# Run association between protein and phenotype using lm per phase/visit 
lm_per_tp_pheno_prot_adj_covar <- function(d_wide, pheno, prot, ph, covariates, scale = F){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
    pheno$TP <- as.numeric(pheno$phase)
    pheno$phase <- NULL
    
    covariates$TP <- as.numeric(covariates$phase)
    covariates$phase <- NULL
  }
  
  d_subs <- inner_join(inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID")),
                       covariates, by = c("SampleID"))
  colnames(d_subs)[1:5] <- c("SampleID", "ID", "TP", "prot", "pheno")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }

  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP","phase")]
  res_table <- data.frame()
  for (tp in unique(d_subs$TP)){

    fo_lm <- as.formula(paste("prot ~ pheno +", paste(covariate_names, collapse = "+")))
    model <- lm(fo_lm, data = d_subs[d_subs$TP == tp,])
    coefs <- summary(model)$coefficients
    res_table <- rbind(res_table, c(ph, prot, tp, coefs['pheno', 1], coefs['pheno', 4]))
  }
  colnames(res_table) <- c("pheno", "prot", "TP","estimate", "pval")
  if (phases){
    res_table <- res_table %>%
      mutate(
        phase = case_when(
            TP == "1" ~ "F",
            TP == "2" ~ "O", 
            TP == "3" ~ "EL",
            TP == "4" ~ "LL",
            .default = "other"
        ),
        .before = TP
      ) 
    res_table$TP <- NULL
  }
  return(res_table)
}

#' Calculate ICC using LMM
#'
#' @param d_wide data frame with proteins (in columns) for all samples (in rows). 
#' @param prot protein name
#' 
get_ICC <- function(d_wide, prot){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    #cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    d_wide$phase <- NULL
  }
  
   d_subs <- d_wide[,c(prot, "ID", "TP")]
   colnames(d_subs)[1] <- "prot"
   
   d_subs$TP <- as.numeric(d_subs$TP)
   d_subs <- na.omit(d_subs)
   d_subs$prot <- scale(d_subs$prot)
   
   m <- lmer(prot ~ 1 + TP + (1|ID), data = d_subs)
   
   vc <- as.data.frame(VarCorr(m))
   var_ID <- vc$vcov[vc$grp == "ID"]  # Variance due to random effect (ID)
   var_residual <- vc$vcov[vc$grp == "Residual"]  # Residual variance
   total_var <- var_ID + var_residual  # Total variance (excluding fixed effects)
   prop_ID <- var_ID / total_var  # Proportion of variance explained by ID
   
   R2m <- performance::r2(m)$R2_marginal
   if (var_ID == 0){
     prop_ID = NA
     R2m = NA
   }
   return (list(ICC = prop_ID, var_tp = R2m))
}



#' Run differential abundance analysis using limma
#'
#' @param joined_data data frame with proteins and covariates together
#' @param tp1 first phase to compare
#' @param tp2 second phase to compare
#' 
run_limma<-function(joined_data, tp1, tp2, all_prots, covariate_names) {
  df <-joined_data[joined_data$phase %in% c(tp1, tp2),]
  df$ID <- as.factor(df$ID)
  df$SampleID <- NULL
  df$phase <- factor(df$phase, levels = c(tp1,tp2))
  
  # design a model 
  formula <- reformulate(termlabels = c("0 + as.factor(phase)", covariate_names), 
                         response = NULL)
  design<-model.matrix(formula, data = df)
  colnames(design)[c(1,2)] <- c("phase1", "phase2")
  
  # specify the pairing
  corfit <- duplicateCorrelation(t(df[,all_prots]), design, block = df$ID)
  
  # make contrast - what to compare
  contrast<- makeContrasts(Diff = phase2 - phase1, levels=design)
  
  # apply linear model to each protein
  # Robust regression provides an alternative to least squares regression that works with less restrictive assumptions. Specifically, it provides much better regression coefficient estimates when outliers are present in the data
  fit<-lmFit(t(df[,all_prots]), design=design,  method="robust", correlation =
               corfit$consensus )
  
  # Extract group means directly from the fit coefficients
  #group_means <- as.data.frame(fit$coefficients)[,c("phase1", "phase2")]
  #colnames(group_means) <- paste0("Adjusted_mean_", colnames(group_means))
  
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=length(all_prots), adjust.method="BH", confint=TRUE)
  #DE_results <- cbind(DE_results, group_means[rownames(DE_results), ])
  #DE_results$Bonferroni_signif <- ifelse(DE_results$P.Value < 0.05 / nrow(DE_results), T, F)
  return(DE_results)
}

#' Run paired wilcoxon test to compare protein levels between 2 phases
#'
#' @param joined_data_adj_covar data frame with proteins adjusted for covariates
#' @param tp1 first phase to compare
#' @param tp2 second phase to compare
#' 
run_wilcox <- function(joined_data_adj_covar, tp1, tp2) {
  joined_data_adj_covar$SampleID <- NULL
  wilcox_pvals <- data.frame(matrix(ncol = 3))
  colnames(wilcox_pvals) <- c("TP1_TP2", "prot", "wilcox_pval")
  cnt <- 1
  for (prot in all_prots){
    df <-joined_data_adj_covar[joined_data_adj_covar$phase %in% c(tp1, tp2), c("ID", "phase", prot)]
    df_wide <- na.omit(my_pivot_wider(df, row_names = "ID", names_from = "phase", values_from = prot))
    pval <- wilcox.test(df_wide[,1], df_wide[,2], paired = T)$p.value
    wilcox_pvals[cnt,] <- c(paste0(tp1, "_", tp2), prot, pval)
    cnt <- cnt + 1
  }
  wilcox_pvals$wilcox_pval <- as.numeric(wilcox_pvals$wilcox_pval)
  wilcox_pvals$BH_qval <- p.adjust(wilcox_pvals$wilcox_pval, method = 'BH')
  #wilcox_pvals$Bonferroni_signif <- ifelse(wilcox_pvals$wilcox_pval < 0.05 / nrow(wilcox_pvals), T, F)
  return(wilcox_pvals)
}

# Run Levene test to comapre protein variances between phases
compare_variances_levene <- function(data) {
  data_long <- data %>%
    pivot_longer(-c(SampleID, ID, phase), 
                 names_to = "protein", values_to = "abundance")
  
  results <- data_long %>%
    group_by(protein) %>%
    summarise(
      levene_p = leveneTest(abundance ~ phase)$`Pr(>F)`[1],
      # Variance statistics
      var_by_phase = list({
        group_by(cur_data(), phase) %>%
          summarise(variance = var(abundance), .groups = 'drop')
      }),
      .groups = 'drop'
    ) %>%
    filter(!is.na(levene_p)) %>%
    mutate(
      adj_p = p.adjust(levene_p, "BH"),
      significant = adj_p < 0.05,
      max_var = map_dbl(var_by_phase, ~max(.x$variance)),
      min_var = map_dbl(var_by_phase, ~min(.x$variance)),
      fold_change = max_var / min_var
    ) %>%
    select(-var_by_phase) %>%
    arrange(levene_p)
  
  return(results)
}

# Get mean protein abundance per phase
get_mean_per_tp <- function(d_wide, prot){
  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  mean_prot_by_TP <- d_subs %>%
    group_by(TP) %>%
    summarize(mean_prot = mean(prot, na.rm = TRUE))
  return(mean_prot_by_TP)
}

#' A custom pivot wider function setting row names
#'
#' @param d long data frame
#' @param row_names column to take the row names from
#' @param names_from column to take the col names from
#' @param values_from column to take the cell values from
#' 
my_pivot_wider <- function(d, row_names, names_from, values_from){
  d2 <- d[,c(row_names, names_from, values_from)] %>%
    pivot_wider(names_from = {{names_from}}, values_from = {{values_from}})
  d2 <- as.data.frame(d2)
  row.names(d2) <- d2[,row_names]
  d2[,row_names ] <- NULL
  return(d2)
}

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

#' Regress covariates using a LMM
#'
#' @param data data frame to regress covariates from
#' @param covar_data data frame with covariates to regress
#' @param covars_longitudinal Logical. True if there are repeated measures
#' @param keep_scale Logical. True if we want to keep the original scale after adjustment

regress_covariates_lmm_phase <- function(data, covar_data, covars_longitudinal = T, keep_scale = F){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(data) & "phase" %in% colnames(data)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    data$TP <- as.numeric(data$phase)
    data$phase <- NULL
  }
  
  if (!"SampleID" %in% colnames(covar_data) & covars_longitudinal) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,"TP"] <- NULL
  
  covar_names = colnames(covar_data)[! colnames(covar_data) %in% c("SampleID", "ID", "TP", "phase")]
  all_pheno = colnames(data)[! colnames(data) %in% c("SampleID", "ID", "TP", "phase")]
  cnt <- 1
  for (ph in all_pheno){
    #print(ph)
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    
    cur_covar_names = covar_names
    if ("batch" %in% colnames(subs)){
      if (table(subs$batch)["batch1"] == 0) cur_covar_names = cur_covar_names[cur_covar_names != "batch"]
    }
    if (length(unique(subs$ID)) == length(subs$ID)) { # if no repeated measurements
      fo_lm <- as.formula(paste("pheno ~ ", paste(cur_covar_names, collapse = "+")))
      lm_fit <- lm(fo_lm, data = subs)
      if (!keep_scale){
        subs[,ph] <- residuals(lm_fit)
      } else { # keep the original scale and global mean
        intercept <- coef(lm_fit)[1]
        subs[,ph] <- subs$pheno - (predict(lm_fit) - intercept)
      }

    } else {
      fo_lmm <- as.formula(paste("pheno ~ ", paste(cur_covar_names, collapse = "+"), "+ (1|ID)"))
      lmm_fit <- lmer(fo_lmm, data = subs)
      if (!keep_scale){
        subs[,ph] <- subs$pheno - lme4:::predict.merMod(lmm_fit, re.form = NA)
      } else { # keep the original scale and global mean
        intercept <- fixef(lmm_fit)[1]
        predicted <- lme4:::predict.merMod(lmm_fit, re.form = NA)
        subs[,ph] <- subs$pheno - (predicted - intercept)
      }
    }
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  
  if(phases){
    cat ("renaming TP to phase\n")
    d_adj <- rename_TP_to_phase(d_adj)
  }
  
  return(d_adj)
}

#' Regress covariates using a LM
#'
#' @param data data frame to regress covariates from
#' @param covar_data data frame with covariates to regress
#' @param covars_longitudinal Logical. True if there are repeated measures
#' @param keep_scale Logical. True if we want to keep the original scale after adjustment

regress_covariates_lm <- function(data, covar_data, covars_longitudinal = T, keep_scale = F){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(data) & "phase" %in% colnames(data)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    data$TP <- as.numeric(data$phase)
    data$phase <- NULL
  }
  
  if (!"SampleID" %in% colnames(covar_data) & covars_longitudinal) {
    covar_data <- cbind(paste0(covar_data$ID, "_",covar_data$TP), covar_data)
    colnames(covar_data)[1] <- "SampleID"
  }
  
  d_adj <- data[,c("SampleID", "ID", "TP")]
  
  data[,"TP"] <- NULL
  covar_data[,"TP"] <- NULL
  
  covar_names = colnames(covar_data)[! colnames(covar_data) %in% c("SampleID", "ID", "TP", "phase")]
  all_pheno = colnames(data)[! colnames(data) %in% c("SampleID", "ID", "TP", "phase")]
  cnt <- 1
  for (ph in all_pheno){
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    
    
    fo_lm <- as.formula(paste("pheno ~ ", paste(covar_names, collapse = "+")))
    lm_fit <- lm(fo_lm, data = subs)
    if (!keep_scale){
      subs[,ph] <- residuals(lm_fit)
    } else { # keep the original scale and global mean
      intercept <- coef(lm_fit)[1]
      subs[,ph] <- subs$pheno - (predict(lm_fit) - intercept)
    }
      
    d_adj <- left_join(d_adj, subs[, c("SampleID", ph)], by = "SampleID")
  }
  
  if(phases){
    cat ("renaming TP to phase\n")
    d_adj <- rename_TP_to_phase(d_adj)
  }
  
  return(d_adj)
}

#' Regress covariates using a LM per phase
#'
#' @param data data frame to regress covariates from
#' @param covar_data data frame with covariates to regress
#' @param covars_longitudinal Logical. True if there are repeated measures
#' @param keep_scale Logical. True if we want to keep the original scale after adjustment

regress_covariates_lm_per_phase <- function(data, covar_data, covars_longitudinal = T, keep_scale = F){
  # if data has phases instead of visits convert phase letter into phase number
  phases = F
  if(! "TP" %in% colnames(data) & "phase" %in% colnames(data)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    data$TP <- as.numeric(data$phase)
    data$phase <- NULL
  }
  
  # Prepare covariate data based on longitudinal structure
  if (covars_longitudinal) {
    # Covariates measured at each timepoint (need SampleID)
    if (!"SampleID" %in% colnames(covar_data)) {
      covar_data <- covar_data %>%
        mutate(SampleID = paste0(ID, "_", TP))
    }
    join_by <- "SampleID"
    covar_data <- covar_data %>% select(-any_of(c("ID", "TP", "phase")))
  } else {
    # Time-invariant covariates (join by ID)
    join_by <- "ID"
    covar_data <- covar_data %>% select(-any_of(c("TP", "phase", "SampleID")))
  }
  
  covar_names <- covar_data %>%
    select(-any_of(c("SampleID", "ID", "TP", "phase"))) %>%
    colnames()
  
  all_pheno <- data %>%
    select(-any_of(c("SampleID", "ID", "TP", "phase"))) %>%
    colnames()
  
  d_adj <- data[, c("SampleID", "ID", "TP")]
  
  for (ph in all_pheno) {
    
    # Merge data with covariates
    merged_data <- data %>%
        select(SampleID, ID, TP, all_of(ph)) %>%
        inner_join(covar_data, by = join_by) %>%
        na.omit()
    colnames(merged_data)[colnames(merged_data) == ph] <- "pheno"
    
    tp_results <- merged_data %>%
      group_by(TP) %>%
      group_modify(~ {
        # Fit linear model for this timepoint
        formula <- as.formula(paste("pheno ~", paste(covar_names, collapse = " + ")))
        lm_fit <- lm(formula, data = .x)
        
        # Calculate adjusted values
        if (!keep_scale) {
          # Return residuals (covariate-adjusted, centered at 0)
          .x$adjusted <- residuals(lm_fit)
        } else {
          # Keep original scale and global mean
          intercept <- coef(lm_fit)[1]
          .x$adjusted <- residuals(lm_fit)+ intercept
        }
        
        # Return SampleID and adjusted values
        .x %>% select(SampleID, adjusted)
      }) %>%
      ungroup() %>%
      # Rename the adjusted column to the phenotype name
      rename_with(~ ph, adjusted) %>%
      select(-TP)
    
    d_adj <- d_adj %>%
      left_join(tp_results, by = "SampleID")
  }
  if (phases) {
    d_adj$ID <- gsub("_.*", "", d_adj$SampleID)
    d_adj$phase <- gsub(".*_", "", d_adj$SampleID)
    d_adj <- d_adj %>%
      dplyr::select(SampleID, ID, phase, everything())
  }
  
  return(d_adj)
}
    
    

# Rename phase or visit number  number into phase letter
rename_TP_to_phase <- function(d) {
  if (! "TP" %in% colnames(d)) {
    cat("error during converting visit to phase: no TP column!\n")
    return (d)
  }
  d %>%
    mutate(
      phase = factor(
        case_when(
          TP == 1 ~ "F",
          TP == 2 ~ "O", 
          TP == 3 ~ "EL",
          TP == 4 ~ "LL",
          .default = "other"
        ),
        levels = c("F", "O", "EL", "LL", "other")  # Specify factor levels
      ),
      .before = TP
    ) %>%
    dplyr::select(-TP)
}

#' Plot 2 trajectories together
plot_together <- function(d_wide = NULL, pheno = NULL, prot, ph, annot = "", method = "gam", scale = T, trajectories = NULL){
  if(! is.null(trajectories)){
    prot_name <- sym(prot)
    ph_name <- sym(ph)
    traj_t <- as.data.frame(t(trajectories)) %>%
      rownames_to_column(var = "TP")
    traj_t$TP <- as.numeric(traj_t$TP)
    g <- ggplot(traj_t) +
      geom_line(aes(x = TP, y = !!prot_name), color = my_colors[2]) +
      geom_line(aes(x = TP, y = !!ph_name), color = my_colors[3]) +
      labs(
        x = "Timepoint",
        y = "",
        title = paste0(
          "<span style='color:", my_colors[2], "'>", prot_name, "</span> vs ",
          "<span style='color:", my_colors[3], "'>", ph_name, "</span>"
        )
      ) +
      theme_minimal() + theme(plot.title = element_markdown())
    return(g)
  }
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    pheno$TP <- as.numeric(pheno$phase)
    d_wide$phase <- NULL
    pheno$phase <- NULL
  }
  
  if (is.character(d_wide$TP)) d_wide$TP <- as.numeric(d_wide$TP)
  d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[,c("ID", "TP", ph)], by = c("ID", "TP"))
  colnames(d_subs) <- c("ID", "TP", "prot", "pheno")
  
  d_subs$TP <- as.numeric(d_subs$TP)
  if (scale){
    d_subs$prot <- scale( d_subs$prot)
    d_subs$pheno <- scale( d_subs$pheno)
  }
  plot_title <- paste0(ph, " - ", prot)
  if (annot != "") plot_title <- paste0(plot_title, ", ", annot)
  if (method == "smooth"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  }  else if (method == "poly3"){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE), color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  } else if(method == 'boxplot'){
    d_subs_long <- d_subs[,-1] %>% pivot_longer(names_to = 'type', cols = c('prot', 'pheno'), values_to = 'value')
    colnames(d_subs_long)[3] <- "value"
    d_subs_long$TP <- as.factor(d_subs_long$TP)
    g <- ggplot(d_subs_long, aes(x = TP,  fill = type, y = value)) +
      geom_boxplot(position = position_dodge(width = 0.85), width = 0.8) + 
      labs(x = "Timepoint ", y = "", title = plot_title) +
      theme_minimal() +
      stat_summary(
        fun = median,
        geom = 'line',
        aes(group = type, color = type),
        position = position_dodge(width = 0.85),
        linewidth = 1
      ) +
      scale_fill_manual(values = c(my_colors[c(2,3)])) + 
      scale_color_manual(values = c('darkgoldenrod4', 'deepskyblue4'))
  } else if (method == 'gam'){
    g <- ggplot(d_subs, aes(x = TP, y = prot)) +
      geom_smooth(method = 'gam', formula=y ~ s(x, k = 4) , color = my_colors[2], aes(x = TP, y = pheno)) +  
      geom_smooth(method = 'gam', formula=y ~ s(x, k = 4), color = my_colors[3], aes(x = TP, y = prot)) +  
      labs(x = "Timepoint ", y = "", 
           title = plot_title) +
      theme_minimal()
  }
  
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
  
}

plot_traj_many_prots2 <- function(prot_trajs = NULL, prots, colored = T, signif = NULL, title = "", phases = F){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  tmp <- prot_trajs[prots,]
  colnames(tmp) <- colnames(prot_trajs)
  d_subs <- tmp %>%
    rownames_to_column(var = 'prot') %>%
    pivot_longer(cols = -prot, names_to = 'TP')
  
  
  d_subs$TP <- as.numeric(d_subs$TP)
  
  if (colored){
    g <- ggplot(d_subs, aes(x = TP, color = prot, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F) +
      theme_minimal()
  } else {
    g <- ggplot(d_subs, aes(x = TP, y = value, group = prot)) +
      geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, alpha = 0.5) +
      theme_minimal() 
    if (!is.null(signif)){
      g <- g + geom_line(data = d_subs[d_subs$prot %in% signif,],aes(x = TP, y = value, group = prot), stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F,  color = 'red')
    }
    
    if (title != ""){
      g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
    }
  }
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
}

plot_traj_prots_and_pheno <- function(d_wide, pheno, prots, ph, title = "", method = 'gam', prot_trajs = NULL, ph_trajs = NULL, phases = F){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  if (!is.null(prot_trajs)) {
    res_trajs <- rbind(ph_trajs, prot_trajs) %>%
      rownames_to_column(var = 'feature')
    
    res_trajs[1, "feature"] <- "pheno"
  } else {
    res_trajs <- data.frame(matrix(nrow = length(prots) + 1, ncol = 101))
    ph_fit <- fit_lmm_poly3_adj_covar(pheno, ph, n = 100, covariates, scale = T)
    res_trajs[1,] <- c("pheno", ph_fit$predicted)
    
    cnt <- 2
    for (prot in prots){
      if (method == 'lmm') {
        fit <- fit_lmm_poly3_adj_covar(d_wide, prot, n = 100, covariates, scale = T)
      } else if (method == 'gam') {
        fit <- gam_prot_tp_adj_covar(d_wide, prot, covariates, scale = T, predict = T)
      } else {
        stop("Error! Wrong method, should be gam or lmm.")
      }
      res_trajs[cnt,] <- c(prot, fit$predicted)
      cnt <- cnt + 1
    }
    colnames(res_trajs) <- c("feature", seq(1,4, length.out = 100))
    
  }
  d_subs <- res_trajs %>%
    pivot_longer(cols = -feature,names_to = 'TP')
  
  
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs$value <- as.numeric(d_subs$value)
  d_subs$feature_type <- ifelse(d_subs$feature == 'pheno', "phenotype" ,"proteins")  
  g <- ggplot(d_subs, aes(x = TP, color = feature_type, y = value, group = feature)) +
    geom_line(stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, linewidth = 0.5, alpha = 0.4) +
    geom_line(data = d_subs[d_subs$feature_type == 'phenotype',], color = my_colors[2], stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F, linewidth = 1) +
    theme_minimal() +
    scale_color_manual(values = my_colors[c(2,3)])
  
  if (title != ""){
    g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
  }
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
}

plot_medians_prots_and_pheno <- function(d_wide, pheno, prots, ph, title = "", scale = T,phases = F){
  #tmp <- as.data.frame(t(apply(prot_trajs[prots,], 1, scale)))
  
  d_subs <- inner_join(pheno[,c("ID", "TP", ph)], d_wide[,c("ID", "TP", prots)], by = c("ID", "TP"))
  colnames(d_subs)[1:3] <- c("ID", "TP", "pheno")
  
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs$ID <- NULL
  
  d_subs <- na.omit(d_subs)
  if (scale){
    d_subs[,-1] <- scale(d_subs[,-1])
  }
  
  medians <- aggregate(. ~ TP, data=d_subs, FUN=median) %>%
    pivot_longer(-TP, names_to = 'feature')
  
  
  medians$TP <- as.numeric(medians$TP)
  medians$value <- as.numeric(medians$value)
  medians$feature_type <- ifelse(medians$feature == 'pheno', "phenotype" ,"proteins")  
  g <- ggplot(medians, aes(x = TP, color = feature_type, y = value, group = feature)) + 
    geom_point() + 
    geom_line(linewidth = 0.5, alpha = 0.4) + 
    geom_line(data = medians[medians$feature_type == 'phenotype',], color = my_colors[2], linewidth = 1) +
    theme_minimal() +
    scale_color_manual(values = my_colors[c(2,3)])
  
  if (title != ""){
    g <- g + ggtitle(title) + theme(plot.title = element_text(size=10)) + theme_minimal()
  }
  if(phases) {
    g <-g + xlab("phase") + scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  }
  g
}

make_radian_plot <- function(d_wide, prots, value = 'mean'){
  library(ggradar)
  if (value == 'mean'){
    d_long <- d_wide %>%
      pivot_longer(names_to = 'Assay', cols = -c('ID', 'TP', 'SampleID'), values_to = 'NPX')
    d_subs <- d_long[d_long$Assay %in% prots,]
    mean_prot_by_TP <- d_subs %>%
      group_by(TP, Assay) %>%
      summarize(value = 1 + mean(NPX, na.rm = TRUE))
  }
  
  data_wide <- mean_prot_by_TP %>%
    pivot_wider(names_from = Assay, values_from = value)
  
  data_wide$TP <- as.factor(data_wide$TP)
  
  g <- ggradar(
    data_wide,
    group.colours = my_colors,
    legend.title = "Timepoints",
    axis.label.size = 2,
    grid.label.size = 4,
    group.line.width = 1,
    group.point.size = 0,
    background.circle.colour = "white",
    gridline.mid.colour = "gray"
  )
  g
}

plot_clusters <- function(cl, method = "", num_k = "", colored = F, signif = NULL, save_pdf = T, add_cluster_name = F, out_path = NA, prot_trajs = NULL){
  plot_list = list()
  for (cluster in unique(cl)){
    title = ifelse(add_cluster_name, paste0(cluster, ", N = ", length(cl[cl == cluster])), "")
    plot_list[[cluster]] <- plot_traj_many_prots2(prot_trajs, names(cl[cl == cluster]), colored = colored, signif, title)
  }
  
  pdf_path = ifelse(is.na(out_path), 
                    paste0("../plots/clustering_signif_v2/", method, "_k", num_k, ".pdf"),
                    out_path)
  
  
  ncols = 4
  nrows = ceiling(length(unique(cl))/4)
  
  if(save_pdf) pdf(pdf_path, width = 4*ncols, height = 4*nrows)
  grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)
  
  
  #if (length(unique(cl)) < 10){
  #  grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)  
  #} else if(length(unique(cl)) < 10) {
  #  grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)  
  #} else {
  #  grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)  
  #}
  if(save_pdf) dev.off()
  
}

plot_association_heatmap <- function(assoc_df, prot_subs, rows = 'pheno', cols = 'prot', vals = 'estimate', signif_vals = 'BH_pval', transpose = F, cutrows = NA, cutcols = NA, cluster_cols = T, col_order = NULL, fontsize = 10){
  assoc_df_wide <- my_pivot_wider(assoc_df[assoc_df$prot %in% prot_subs,], rows, cols, vals)
  signif_labels <- my_pivot_wider(assoc_df[assoc_df$prot %in% prot_subs,], rows, cols, signif_vals)
  signif_labels <- ifelse(signif_labels < 0.05, "*", "")
  
  if (transpose){
    assoc_df_wide <- as.data.frame(t(assoc_df_wide))
    signif_labels <- as.data.frame(t(signif_labels))
    fontsize_row = fontsize
    fontsize_col = 10
  } else {
    fontsize_col = fontsize
    fontsize_row = 10
  }
  if (!is.null(col_order)) {
    assoc_df_wide <- assoc_df_wide[, col_order, drop = FALSE]
    signif_labels <- signif_labels[, col_order, drop = FALSE]
  }
  
  max_val <- max(abs(min(assoc_df_wide)), max(assoc_df_wide))
  breaksList = seq(-max_val, max_val, by = 0.01)
  if(!0 %in% breaksList) breaksList <- sort(c(breaksList, 0))
  
  full_palette <- rev(brewer.pal(n = 11, name = "RdYlBu"))
  full_palette[ceiling(length(full_palette)/2)] <- "#FFFFFF"
  colorList <- colorRampPalette(full_palette)(length(breaksList))
  
  h <- pheatmap::pheatmap(assoc_df_wide, display_numbers = signif_labels,  
                fontsize_col = fontsize_col, fontsize_row = fontsize_row,
                color = colorList, breaks = breaksList, cutree_rows = cutrows, 
                cutree_cols = cutcols, cluster_cols = cluster_cols)
  
  h
}

plot_association_volcano <- function(assoc_df){
  if ("P4" %in% assoc_df$pheno) desired_order <- c("P4", "E2", "LH", "FSH","PRL")
  if (!"P4" %in% assoc_df$pheno) desired_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
  ggplot(assoc_df, aes(x = estimate, y = -log10(BH_pval))) +
    geom_point(aes(color = BH_pval < 0.05), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("grey", "red"), 
                       labels = c("FALSE" = "Not significant", "TRUE" = "FDR < 0.05")) +
    ggrepel::geom_text_repel(
      data = subset(assoc_df, BH_pval < 0.05),
      aes(label = prot), 
      size = 2,
      max.overlaps = 10
    ) +
    facet_wrap(~ factor(pheno, levels = desired_order), ncol = 2) +
    theme_bw() +
    labs(color = "Significance") + 
    theme(legend.position = "none")
}

# scatter colored by visit to see the relationship at each visit
scatter_col_tp <- function(d_wide, pheno, prot, ph, scale = F, add_points = F){
  phases = F
  if(! "TP" %in% colnames(d_wide) & "phase" %in% colnames(d_wide)){
    phases = T
    cat("Working with phases not visit numbers!\n")
    d_wide$TP <- as.numeric(d_wide$phase)
    pheno$TP <- as.numeric(pheno$phase)
    d_wide$phase <- NULL
    pheno$phase <- NULL
  }
  
  if ("SampleID" %in% colnames(d_wide)){
    d_subs <- inner_join(d_wide[,c("SampleID", "ID", "TP", prot)], pheno[,c("SampleID" ,ph)], by = c("SampleID"))
    colnames(d_subs) <- c("SampleID", "ID", "TP", "prot", "pheno")
  } else {
    d_subs <- inner_join(d_wide[,c("ID", "TP", prot)], pheno[,c("ID", "TP" ,ph)], by = c("ID", "TP"))
    colnames(d_subs) <- c("ID", "TP", "prot", "pheno")
  }
  d_subs$TP <- as.factor(d_subs$TP)
  d_subs <- na.omit(d_subs)
  
  if(scale){
    d_subs$pheno <- scale(d_subs$pheno)
    d_subs$prot <- scale(d_subs$prot)
  }
  g <- ggplot(d_subs, aes(x = prot, y = pheno, colour = TP)) + 
    geom_smooth(method = 'lm', alpha = 0.2) + 
    stat_smooth(method = 'lm', se = F) +
    theme_minimal() +
    labs(x = prot, y = ph, 
         title = paste0(ph, " - ", prot))  +
    scale_color_manual(values = my_colors)
  
  if (add_points) g <- g + geom_point(alpha = 0.5)
  if(phases) {
    g <-g + scale_color_manual(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL"),
      values = my_colors
    ) + labs(colour="Phase")
  }
  g
}


plot_boxplot_with_traj <- function(d_adj, prot) {
  d_adj_gam <- na.omit(d_adj[,c("SampleID", "ID", "phase", prot)])
  colnames(d_adj_gam)[4] <- "prot"
  if("phase" %in% colnames(d_adj_gam)) {
    d_adj_gam$TP <- as.numeric(d_adj_gam$phase)
  }
  d_adj_gam$ID <- as.factor(d_adj_gam$ID)

  model <- gam(prot ~ s(TP, k = 4) + s(ID, bs = 're'), 
               data = d_adj_gam, method = 'REML')
  
  new_data <- data.frame(
    TP = seq(1, 4, length.out = 20),
    ID = unique(d_adj_gam$ID)[1] 
  )
  
  predictions <- predict(model, newdata = new_data, 
                         exclude = "s(ID)", se.fit = TRUE)
  
  traj <- data.frame(
    TP = new_data$TP,
    pheno = predictions$fit,
    lower = predictions$fit - 1.96 * predictions$se.fit,
    upper = predictions$fit + 1.96 * predictions$se.fit
  )
  
  p <- ggplot(d_adj_gam) +
    geom_boxplot(aes(x = as.numeric(TP), y = prot, group = TP), 
                 width = 0.1, color = my_colors[3], outliers = FALSE) +
    geom_line(data = traj, aes(x = TP, y = pheno), color = my_colors[4]) +
    geom_ribbon(data = traj, aes(x = TP, ymin = lower, ymax = upper), 
                alpha = 0.2, fill = my_colors[4]) +
    geom_jitter(aes(x = as.numeric(TP), y = prot), 
                alpha = 0.2, width = 0.1, color = my_colors[3]) +
    labs(x = "Phase", y = paste0(prot, "\n(covariate-adjusted)")) +
    theme_minimal() +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL")
    )
  
  return(p)
}

rename_p4_e2 <- function(vector){
  vector <- gsub("^17BES$","E2", vector)
  vector <- gsub("^PROG$","P4", vector)
  return(vector)
}
