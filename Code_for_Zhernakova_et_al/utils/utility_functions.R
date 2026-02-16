
my_colors <- c("#eddb6d", "#ed9f47", "#4b9aaf", "#3a6887")

library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyverse)
library(mgcv)
library(ggtext)
library(RColorBrewer)
library(limma)

######################################################
#####           ASSOCIATION ANALYSIS             #####
######################################################

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
  d_subs$ID <- as.factor(d_subs$ID)
  
  if (scale) {
    d_subs$prot <- scale(d_subs$prot)
    d_subs$pheno <- scale(d_subs$pheno)
  }
  
  covariate_names <- colnames(covariates)[! colnames(covariates) %in% c("SampleID", "ID", "TP", "phase")]
  
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

######################################################
#####                     ICC                    #####
######################################################

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
   
   return (list(ICC = prop_ID, var_tp = R2m))
}

######################################################
#####                  DAP analysis              #####
######################################################

#' Run differential abundance analysis using limma
#'
#' @param joined_data data frame with proteins and covariates together
#' @param tp1 first phase to compare
#' @param tp2 second phase to compare
#' 
run_limma<-function(joined_data, tp1, tp2) {
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
  # apply contrast
  contrast_fit<-contrasts.fit(fit, contrast)
  # apply empirical Bayes smoothing to the SE
  ebays_fit<-eBayes(contrast_fit)
  # summary
  print(summary(decideTests(ebays_fit)))
  # extract DE results
  DE_results<-topTable(ebays_fit, n=length(all_prots), adjust.method="BH", confint=TRUE)
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

######################################################
#####                   HELPERS                  #####
######################################################

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
  cnt <- 1
  for (ph in colnames(data)[3: (ncol(data))]){
    if (covars_longitudinal){
      covar_data$ID <- NULL
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "SampleID"))
    } else {
      subs <- na.omit(inner_join(data[, c("ID","SampleID", ph)], covar_data, by = "ID"))
    }
    colnames(subs)[3] <- 'pheno'
    
    if (length(unique(subs$ID)) == length(subs$ID)) { # if no repeated measurements
      fo_lm <- as.formula(paste("pheno ~ ", paste(covar_names, collapse = "+")))
      lm_fit <- lm(fo_lm, data = subs)
      if (!keep_scale){
        subs[,ph] <- residuals(lm_fit)
      } else { # keep the original scale and global mean
        intercept <- coef(lm_fit)[1]
        subs[,ph] <- subs$pheno - (predict(lm_fit) - intercept)
      }

    } else {
      fo_lmm <- as.formula(paste("pheno ~ ", paste(covar_names, collapse = "+"), "+ (1|ID)"))
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


######################################################
#####                  PLOTTING                  #####
######################################################

# plot the association result as a heatmap
plot_association_heatmap <- function(assoc_df, prot_subs, rows = 'pheno', cols = 'prot', vals = 'estimate', signif_vals = 'BH_pval', transpose = F, cutrows = NA, cutcols = NA, cluster_cols = T, col_order = NULL){
  assoc_df_wide <- my_pivot_wider(assoc_df[assoc_df$prot %in% prot_subs,], rows, cols, vals)
  signif_labels <- my_pivot_wider(assoc_df[assoc_df$prot %in% prot_subs,], rows, cols, signif_vals)
  signif_labels <- ifelse(signif_labels < 0.05, "*", "")
  
  if (transpose){
    assoc_df_wide <- as.data.frame(t(assoc_df_wide))
    signif_labels <- as.data.frame(t(signif_labels))
    fontsize_row = 8
    fontsize_col = 10
  }
  if (!is.null(col_order)) {
    assoc_df_wide <- assoc_df_wide[, col_order, drop = FALSE]
    signif_labels <- signif_labels[, col_order, drop = FALSE]
  }
  
  max_val <- max(abs(min(assoc_df_wide)), max(assoc_df_wide))
  breaksList = seq(-max_val, max_val, by = 0.01)
  colorList <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
  h <- pheatmap(assoc_df_wide, display_numbers = signif_labels, fontsize_number = 10, 
                fontsize_col = fontsize_col, fontsize_row = fontsize_row,
                color = colorList, breaks = breaksList, cutree_rows = cutrows, 
                cutree_cols = cutcols, cluster_cols = cluster_cols)
  
  h
}

# Plot the association results as a volcano plot
plot_association_volcano <- function(assoc_df){
  if ("PROG" %in% assoc_df$pheno) desired_order <- c("PROG", "17BES", "LH", "FSH","PRL")
  if (!"PROG" %in% assoc_df$pheno) desired_order <- c("ALT", "AST", "TRI", "HDL", "COL", "LDL", "INS", "HOMA_B", "HOMA_IR", "GL")
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
      max.overlaps = 20
    ) +
    facet_wrap(~ factor(pheno, levels = desired_order)) +
    theme_minimal() +
    labs(color = "Significance")
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
  
  if (add_points) g <- g + geom_point()
  if(phases) {
    g <-g + scale_color_manual(
      breaks = c(1, 2, 3, 4),
      labels = c("F", "O", "EL", "LL"),
      values = my_colors
    )
  }
  g
}


# Plots protein/phenotype boxplots per phase together with predicted trajectories and points
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
