# Rscript
# This script run the linear models between hormones and bacteria
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 31/01/2025
# R version: R v4.4.1


library(ggplot2)
library(lme4)
library(lmerTest)
library(bestNormalize)

#get the metadata
df_feno <- read.csv('data/clean_phenotypes.csv')
df_feno <- df_feno[,-1]
df_hormones <- df_feno[,c('PROG','LH','FSH','BES17','PRL','Code')]

quest <- read.csv('data/clean_questionnaire.csv')
quest <- quest[,c('Code','Visit_number','Age','Smoking_habits','BMI_ranges','Age_category','BMI','Pregnancy_category')]

df_lab <- read.csv('data/Laboratory_data.csv')

#function to remove the outliers relative to a specific column
remove_outliers <- function(data, column_name) {
  
  Q1 <- quantile(data[[column_name]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column_name]], 0.75, na.rm = TRUE)
  
  IQR_val <- Q3 - Q1
  
  lower_bound <- Q1 - 4 * IQR_val
  upper_bound <- Q3 + 4 * IQR_val
  
  filtered_data <- data[data[[column_name]] >= lower_bound & data[[column_name]] <= upper_bound, ]
  
  return(filtered_data)
}

for(col in colnames(df_hormones)[1:5]){
  df_hormones <- remove_outliers(df_hormones, col)
}

complete_data <- merge(df_lab, df_hormones, by = 'Code')
complete_data <- merge(complete_data,quest, by = 'Code')

#create a function to fit the linear mixed model
fit_function_gaussian <- function(abundance_data,covariates,use_zeros = TRUE){
  
  #prepare data
  df <- merge(abundance_data,complete_data, by = 'Code')
  bact_to_check <- colnames(abundance_data[,-which(colnames(abundance_data) == 'Code')])
  
  model_results <- list()
  #iteration over hormones
  for(horm in colnames(df_hormones[,-6])){
    
    res <- matrix(ncol = 9, nrow = 0)  # Initialize matrix for the results of the current hormone
    colnames(res) <- c('Bacteria','Hormone','Estimate','Intercept','Std_error',
                       'p.value','lower','upper','Status')
    #iteration over bacteria
    for(bact in bact_to_check) {
      
      df_temp <- df[,c(horm,bact,'Code','Woman','Pregnancy_category','Visit_number','Qubit_DNA','Qubit_Library','Total_counts')]
      b <- bestNormalize::orderNorm(df_temp[[horm]]) #normalize the hormones
      df_temp[[horm]] <- b$x.t
      
      if(!use_zeros){ #the parameter tell us whether use the zero samples or not
        intersect_genus <- intersect(colnames(abundance_data),colnames(rel_ab_genus))
        intersect_species <- intersect(colnames(abundance_data),colnames(rel_ab_species))
        if(length(intersect_genus) !=0){
          rel_ab <- rel_ab_genus[,intersect_genus]
        }else{
          rel_ab <- rel_ab_species[,intersect_species]
        }
        zeros_samples <- rownames(subset(rel_ab, rel_ab[[bact]] != 0)) #grab zeros samples
        df_temp <- subset(df_temp, !(Code %in% zeros_samples)) #take just the non-zeros samples
      }
      
      df_temp <- na.omit(df_temp) #remove the NA
      
      # Define formula
      formula <- paste0(horm,' ~ ', bact, ' + (1|Woman)')
      formula <- paste0(formula,' + ',paste0(covariates,collapse = ' + '))
      cat(crayon::strip_style(formula)) #display the formula
      cat('\n')
      
      if(nrow(df_temp) == 0){
        cat(crayon::bgBlue$bold$yellow("<<>> no samples available <<>>\n"))
        cat(crayon::bgBlue$bold$yellow("Skipping...")) 
        cat('\n')
        next
      }
      
      df_temp$Visit_number <- as.numeric(df_temp$Visit_number) #make the visit column numerical
      
      # Scale the numeric columns
      df_temp[,bact] <- scale(df_temp[,bact]) #scale the bacteria column
      
      
      if(nrow(df_temp) == 0){ #if there's not sample go to the next bacteria
        cat(crayon::bgBlue$bold$yellow("<<>> no samples available <<>>\n"))
        cat(crayon::bgBlue$bold$yellow("Skipping...")) 
        next
      }
      
      # Fit the model and handle warnings/errors
      model <- tryCatch({
        lmer(formula = as.formula(formula), data = df_temp)
      }, warning = function(w) {
        cat(crayon::red('Warning here'))
        cat('\n')
        return(NULL)
      },message = function(e) {
        cat(crayon::red('Warning here'))
        cat('\n')
        return(NULL)
      },
      error = function(e) {
        cat(crayon::red('Error here'))
        cat('\n')
        return(NULL)
      })
      if (is.null(model)) {
        cat(crayon::red('FAILURE'))  # Mark as failure
        cat('\n')
        row <- c(bact, horm, rep(NA, 6), 'Failure', nrow(df_temp))
      } else {
        cat(crayon::green('SUCCESS'))
        cat('\n')
        s <- summary(model)
        d <- confint(model, method = 'Wald')
        d <- d[bact, ]
        
        # Get the model results
        est <- s$coefficients[bact, 'Estimate']
        int <- s$coefficients['(Intercept)', 'Estimate']
        std <- s$coefficients[bact, 'Std. Error']
        p <- s$coefficients[bact, 'Pr(>|t|)']
        low <- d[1]
        upp <- d[2]
        
        # Store the results in a row and add it to the matrix
        row <- c(bact, horm, est, int, std, p, low, upp, 'Success')
      }
      
      res <- rbind(res, row)  # Append the row to the results matrix
    }
    
    res <- as.data.frame(res)  # Convert the matrix to data frame
    model_results[[horm]] <- res  # Store the results in the list
  }
  to_ret <- do.call(rbind,model_results) #create the dataframe to return
  
  return(to_ret)
}

#to add the permuted p value column
permut_analysis <- function(results, n_perm = 1000, covariates) {
  
  cat(crayon::bgYellow$bold$blue("<<<< Randomization here >>>>")) #display it's going to fitting randomized models
  cat('\n')
  
  results$p.value <- as.numeric(results$p.value) #make the p-value column as numerical
  
  # Filter significant results
  sig_results <- subset(results, !is.na(p.value) & (p.value < 0.05))
  
  #create an association column
  sig_results$association <- paste0(sig_results$Hormone, '-', sig_results$Bacteria)
  
  # Create a list where names are hormones and elements are the associated bacteria
  hormone_bacteria_list <- split(sig_results$Bacteria, sig_results$Hormone)
  
  # Remove duplicates from each list element
  hormone_bacteria_list <- lapply(hormone_bacteria_list, unique)
  
  # GENUS
  genus_to_check <- unique(intersect(unlist(unname(hormone_bacteria_list)), colnames(clr_genus)))
  df_genus <- merge(clr_genus, complete_data, by = 'Code')
  
  # SPECIES
  species_to_check <- unique(intersect(unlist(unname(hormone_bacteria_list)), colnames(clr_species)))
  df_species <- merge(clr_species, complete_data, by = 'Code')
  
  permutation_results <- list()
  
  # Function to perform analysis for a given level (genus or species)
  analyze_level <- function(level, df, taxa_to_check) {
    for (bact in taxa_to_check) {

      # Get hormones to test for the current bacterium
      hormones_to_test <- names(hormone_bacteria_list)[sapply(hormone_bacteria_list, function(x) any(x == bact))]
      
      for (horm in hormones_to_test) {
        
        # Get the original p-value from the input results
        original_p <- sig_results$p.value[sig_results$Hormone == horm & sig_results$Bacteria == bact]
        
        # Build the formula for the model
        formula <- paste0(horm, ' ~ ', bact, ' + (1|Woman)')
        formula <- paste0(formula, ' + ', paste0(covariates, collapse = ' + '))
        cat(crayon::bgYellow$bold$blue(formula))
        cat('\n')
        
        # Perform permutations
        permuted_p_values <- numeric(n_perm)
        
        for (n in 1:n_perm) {
          # Permute the bacterium column
          df[[bact]] <- scale(df[[bact]])
          df[[bact]] <- sample(df[[bact]])
          
          # Fit the permuted model
          permuted_model <- tryCatch({
            glmer(formula = as.formula(formula), data = df, family = Gamma(link = "log"),
                  glmerControl(optimizer = 'bobyqa'))
          }, warning = function(w) {
            cat(crayon::red('Warning here'))
            cat('\n')
            return(NULL)
          }, error = function(e) {
            cat(crayon::red('Warning here'))
            cat('\n')
            return(NULL)
          })
          
          if (!is.null(permuted_model)) {
            permuted_p_values[n] <- summary(permuted_model)$coefficients[bact, 'Pr(>|z|)']
          } else {
            permuted_p_values[n] <- NA  # Mark failures as NA
          }
        }
        
        # Count how many permuted p-values are less than the original p-value
        p_less_count <- sum(permuted_p_values < original_p, na.rm = TRUE)
        
        # Calculate the permutation-based significance
        permutation_significance <- p_less_count / length(permuted_p_values[!is.na(permuted_p_values)])
        
        # Append results
        permutation_results <<- rbind(
          permutation_results,
          data.frame(
            Hormone = horm,
            Bacteria = bact,
            p.permuted = permutation_significance
          )
        )
      }
    }
  }
  
  # Perform analysis for genus
  analyze_level("Genus", df_genus, genus_to_check)
  
  # Perform analysis for species
  analyze_level("Species", df_species, species_to_check)
  
  # Merge permutation results back into the original results dataset
  results <- merge(results, permutation_results, by = c("Hormone", "Bacteria"), all.x = TRUE)
  
  return(results)
} 

#get the corrected (by tech vars) abundances
clr_genus <- read.csv('data/corrected_clr/filtered/genus/1.Tech.csv')[,-1]
clr_species <- read.csv('data/corrected_clr/filtered/species/1.Tech.csv')[,-1]

#fit model for genus
genus_res_gauss <- fit_function_gaussian(clr_genus,covariates = c('Pregnancy_category','Visit_number'))

#fit model for species
species_res_gauss <- fit_function_gaussian(clr_species,covariates = c('Pregnancy_category','Visit_number'))

#wrap the results of genus and species
linear_results <- rbind(genus_res_gauss,species_res_gauss)

results_divided_by_hormones <- list()
for(horm in c('PROG','LH','FSH','BES17','PRL','Code')){
  dummy_df <- subset(linear_results, Hormone == horm)
  dummy_df$p.adj <- p.adjust(dummy_df$p.value, method = 'BH')
  results_divided_by_hormones[[horm]] <- dummy_df
}

linear_results <- do.call(rbind, results_divided_by_hormones)

#add the permuted pvalue column
linear_results <- permut_analysis(linear_results, covariates = c('Pregnancy_category','Visit_number'))


#### PLOT RESULTS ####

# Convert columns to numeric
for(col in c('p.value','p.adj','Estimate','upper','lower','Std_error','Intercept', 'p.permuted')){
  linear_results[[col]] <- round(as.numeric(linear_results[[col]]),5)
}
 
write.csv(x = linear_results,file = 'results/linear_models/results_hormones.csv')
 
#Filter significant rows
res_df <- subset(linear_results, !is.na(p.value) & p.value < 0.05)
  
#Add asterisks for significant adjusted p-values
res_df$association <- ifelse(res_df$p.adj < 0.2, 
                               paste0(res_df$Hormone, '-', res_df$Bacteria, '*'), 
                               paste0(res_df$Hormone, '-', res_df$Bacteria))
  
res_df$association <- gsub("\\.", " ", res_df$association)
res_df$association <- sub("Unclassified.*$", "spp_uncl", res_df$association)
#Create the plot
g <- ggplot(data = res_df, aes(y = factor(association, levels = association[order(Bacteria)]), x = Estimate, xmin = lower, xmax = upper)) +
  geom_point(aes(color = Hormone), size = 3) +
  geom_errorbarh(height = 0.13, aes(color = Hormone), linewidth = 1.5) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  xlab('Coefficient') + ylab('') +
  
  # Adjust text size in geom_text to match 5-7 pt
  geom_text(aes(label = sprintf("%.3f", p.value)), hjust = 0.5, vjust = -1.5, size = 2.5) +  # Approx. 5 pt
  geom_text(aes(label = sprintf("%.3f", p.permuted)), hjust = 0.5, vjust = 1.6, size = 2.5, color = 'blue') +  # Approx. 5 pt
  
  scale_color_manual(values = c("BES17" = "#FFA07A", "FSH" = "#89CFF0", "LH" = "#B3E5B6", "PRL" = "#FFD700", 'PROG' = "#D1B3E5"), 
                     name = "Hormone") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7,family = 'sans'),  # Keep within 5-7 pt
    axis.text.y = element_text(face = 'italic', size = 7, family = 'sans'),
    axis.title = element_text(size = 7,family = 'sans'),
    legend.title = element_text(size = 7,family = 'sans'),
    legend.text = element_text(size = 7,family = 'sans'),
    panel.border = element_blank()
  )

png("results/linear_models/results.png",
    height=3000,width=2000,
    res = 300)
plot(g)
dev.off()

pdf("results/linear_models/results.pdf")
plot(g)
dev.off()
