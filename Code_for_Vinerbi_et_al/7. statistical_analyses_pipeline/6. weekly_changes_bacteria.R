# Rscript
# This script run the linear models between hormones and visit number
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 31/01/2025
# R version: R v4.4.1



library(ggplot2)
library(lme4)
library(MASS)
library(bestNormalize)
library(lmerTest)
library(gridExtra)
library(foreach)
library(doParallel)
library(glmmTMB)


# GET DATA
setwd('data/corrected_clr/filtered/species/')
list_of_files <- list.files()

get_data <- function(n){
  #take one adjusted dataframe to grab names of the bacteria that are in more than 20% of the samples
  dummy_data <- read.csv(list_of_files[2])[,-c(1,ncol(read.csv(list_of_files[2])))] #remove 1th and last column (Code)
  bacteria <- colnames(dummy_data) #take names of bacteria
  data <- read.csv(list_of_files[n]) #take the right dataframe
  if(n == 1){
    colnames(data)[1] <- 'Code' #in the clr dataframe the first column represents Code
    data <- data[,c('Code', bacteria)] #take just the bacteria determined above
  } else {
    data <- data[,-1] #in adjusted dataframes we have to remove the first column (1,2,...)
  }
  
  data$Visit_number <- sapply(data$Code,function(stringa){substring(stringa,6,6)})
  data$Woman <- sapply(data$Code,function(stringa){substring(stringa,1,4)})
  data <- data[,c('Code','Visit_number','Woman',
                  setdiff(colnames(data),c('Code','Visit_number','Woman')))]
  
  return(data)
} #function to grab the correct dataframe

#titles
title_of_plots_short <- c(
  '1' = 'No Corr',
  '2' = 'Tech Vars Corr',
  '3' = 'Tech-Age',
  '4' = 'Tech-Age-Preg',
  '5' = 'Tech-Age-Preg-Pill',
  '6' = 'Tech-Age-Preg-Pill-SwabM',
  '7' = 'Tech-Age-Preg-Pill-SwabM/A',
  '8' = 'Tech-Age-Preg-Pill-SwabM/A-WHR',
  '9' = 'Tech-Age-Preg-Pill-SwabM/A-WHR-Bris',
  '10' = 'Tech-Age-Preg-Pill-SwabM/A-WHR-Bris-Coitus'
)

check_the_variation <- function(n, num_perm) {
  
  abundance <- get_data(n)  # get the data
  
  model_results <- matrix(ncol = 10)  # initialize the matrix
  colnames(model_results) <- c('Bacteria', 'Estimate', 'Std.Error', 'p.value', 
                               'p.permuted', 'Lower', 'Upper', 'ICC', 'Status', 'Correction')
  
  for (bact in colnames(abundance[, -(1:3)])) {
    df <- abundance[, c(bact, 'Code', 'Woman', 'Visit_number')]  # take just the bacterium under analysis
    df$Visit_number <- as.numeric(df$Visit_number)  # convert visit to numeric
    
    model <- tryCatch({
      print(paste(bact, '~ Visit_number'))
      lmer(as.formula(paste(bact, '~ Visit_number + (1|Woman)')), data = df)
    }, warning = function(x) {
      cat(crayon::red('Warning here\n'))
      return(NULL)
    }, message = function(x) {
      cat(crayon::red('Message here\n'))
      return(NULL)
    }, error = function(x) {
      cat(crayon::red('Error here\n'))
      return(NULL)
    })
    
    if (!is.null(model)) {
      cat(crayon::green('SUCCESS\n'))
      s <- summary(model)
      coeffs <- s$coefficients['Visit_number', c('Estimate', 'Std. Error', 'Pr(>|t|)')]
      conf_interval <- confint(model, parm = "Visit_number", level = 0.95)
      icc_value <- performance::icc(model)$ICC_adjusted
      
      if (coeffs['Pr(>|t|)'] < 0.05) {
        p_less_perms <- sapply(1:num_perm, function(perm) {
          cat(crayon::bgYellow$bold$red("<<>> randomization is happening here <<>>\n"))
          df[[bact]] <- sample(df[[bact]])
          model_perm <- tryCatch({
            lmer(as.formula(paste(bact, '~ Visit_number + (1|Woman)')), data = df)
          }, warning = function(x) {
            cat(crayon::red('Warning here\n'))
            return(NULL)
          }, message = function(x) {
            cat(crayon::red('Message here\n'))
            return(NULL)
          }, error = function(x) {
            cat(crayon::red('Error here\n'))
            return(NULL)
          })
          
          if (!is.null(model_perm)) {
            cat(crayon::green('SUCCESS\n'))
            s_perm <- summary(model_perm)
            coeffs_perm <- s_perm$coefficients['Visit_number', 'Pr(>|t|)']
            return(coeffs_perm <= coeffs['Pr(>|t|)'])
          } else {
            return(NA)
          }
        })
        
        # Filter out NA values and compute pseudo p-value
        valid_perms <- p_less_perms[!is.na(p_less_perms)]
        pseudo_pvalue <- ifelse(length(valid_perms) > 0, sum(valid_perms) / length(valid_perms), NA)
      } else {
        pseudo_pvalue <- NA
      }
      
      row <- c(bact, coeffs['Estimate'], coeffs['Std. Error'], coeffs['Pr(>|t|)'], 
               pseudo_pvalue, conf_interval[1], conf_interval[2], icc_value, 
               "Success", title_of_plots_short[n])
    } else {
      cat(crayon::red('FAILURE\n'))
      row <- c(bact, rep(NA, 6), "Failure", title_of_plots_short[n])
    }
    
    model_results <- rbind(model_results, row)
  }
  
  model_results <- as.data.frame(model_results[-1, ])
  rownames(model_results) <- NULL
  
  # Mark significant results
  model_results$p.adj <- p.adjust(as.numeric(model_results$p.value), method = 'fdr')
  model_results$significance <- ifelse(as.numeric(model_results$p.value) < 0.05, 
                                       "Significant", "Not Significant")
  
  return(model_results)
}
complete_function <- function(num_perm){
  res_temp <- list()
  for(corr in 1:length(title_of_plots_short)){
    res <- check_the_variation(n = corr,num_perm = num_perm)
    res_temp[[corr]] <- res
  }
  complete_res <- do.call(rbind,res_temp)
  return(complete_res)
}

complete_results <- complete_function(1000)

for(col in c('Estimate','Std.Error','p.value','p.permuted','Lower','Upper','ICC','p.adj')){
  complete_results[[col]] <- round(as.numeric(complete_results[[col]]),5)
}


complete_results$Correction <- factor(complete_results$Correction, levels = as.character(title_of_plots_short))

complete_results$Correction_Number <- as.factor(as.numeric(factor(complete_results$Correction)))

# Replace dots with spaces for consistency
complete_results$Bacteria <- gsub("\\.", " ", complete_results$Bacteria)

# Handle the case with three "Unclassified" levels (family, genus, species)
complete_results$Bacteria <- sub("^Unclassified Unclassified Unclassified$", "f_uncl g_uncl sp_uncl", complete_results$Bacteria)

# Handle the case with a single "Unclassified" (species)
complete_results$Bacteria <- sub(" Unclassified$", " sp_uncl", complete_results$Bacteria)

g_combined <- ggplot(complete_results, aes(x = Bacteria, y = Estimate, color = significance)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), linewidth = 1, width = 0.2) +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed', alpha = 0.2) +
  coord_flip() +
  scale_color_manual(values = c("Significant" = "#ffa32b", "Not Significant" = "#666460")) +
  
  # Adjust text size in geom_text() to match 5-7 pt
  geom_text(aes(label = sprintf("%.3f", p.value)),
            hjust = 0.5, vjust = -1.5, size = 2.5, color = 'black') +  # Reduced to ~5 pt
  geom_text(aes(label = sprintf("%.3f", p.permuted)),
            hjust = 0.5, vjust = 2, size = 2.5, color = 'blue') +  # Reduced to ~5 pt
  
  labs(x = "", y = "Coefficient", color = "Significance") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 7, family = 'sans'), 
    axis.text.y = element_text(size = 7, face = 'italic', family = 'sans'),
    axis.title.x = element_text(size = 7, family = 'sans'),
    legend.position = 'bottom',
    legend.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank()
  ) +
  facet_wrap(~ Correction_Number, nrow = 1)


complete_results$ICC <- as.numeric(complete_results$ICC)

write.csv(x = complete_results[,-8],file = 'results/linear_models/results_across_week.csv')
write.csv(x = subset(complete_results, Correction == 'Tech Vars Corr')[,c(1,8)],'results/linear_models/icc_table.csv')


png("results/linear_models/bacteria_weekly_changes.jpg",height=3500,width=4000,res = 300)
plot(g_combined)
dev.off()

pdf("results/linear_models/bacteria_weekly_changes.pdf",width = 12,height = 10.5)
plot(g_combined)
dev.off()

