# Rscript
# This script run the linear models between hormone changes and bacteria abundance changes
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 31/01/2025
# R version: R v4.4.1



library(lme4)
library(lmerTest)
library(ggplot2)

#get the data
df_feno <- read.csv('data/clean_phenotypes.csv')[,c('PROG','PRL','LH','FSH','BES17','Code')]

df_feno$Visit_number <- sapply(df_feno$Code, function(x){substring(x,6,6)})
df_feno$Woman <- sapply(df_feno$Code, function(x){substring(x,1,4)})

clr_species <- read.csv('data/corrected_clr/filtered/species/1.Tech.csv')[,c('Lactobacillus.crispatus',
                                                                                                 'Lactobacillus.gasseri',
                                                                                                 'Lactobacillus.iners',
                                                                                                 'Lactobacillus.jensenii',
                                                                                                 'Streptococcus.agalactiae','Code')]
#add useful columns to the data
clr_species$Woman <- sapply(clr_species$Code,function(x){substring(x,1,4)})
clr_species$Visit_number <- sapply(clr_species$Code,function(x){substring(x,6,6)})

species_to_test <- c('Lactobacillus.crispatus',
                     'Lactobacillus.gasseri',
                     'Lactobacillus.iners',
                     'Lactobacillus.jensenii',
                     'Streptococcus.agalactiae')

woman_to_analyse <- clr_species$Woman
df_feno <- subset(df_feno,Code %in% clr_species$Code)


#define a formula to calculate variations
calculate_deltas <- function(df, exclude_col) {
  dummy_df <- data.frame(Woman = woman_to_analyse)
  
  #Identify the columns to process (excluding specified columns)
  amount_columns <- colnames(df[, -exclude_col])
  
  #Create a copy of the original dataframe to store the deltas
  df_result <- df
  
  #Loop over the columns (amounts) to calculate deltas
  for (amount in amount_columns) {
    #Loop over visits 1 to 3 (assuming visits are sequential)
    for (vis in 1:3) {
      #Subset data for the current and next visit
      df_start <- subset(df, df[['Visit_number']] == vis)
      df_start <- merge(df_start, dummy_df, by = 'Woman',all.y = T)
      
      df_end <- subset(df, df[['Visit_number']] == (vis + 1))
      df_end <- merge(df_end, dummy_df, by = 'Woman',all.y = T)

      #Ensure data is ordered by 'Woman'
      df_start <- df_start[order(df_start[['Woman']]), ]
      df_end <- df_end[order(df_end[['Woman']]), ]

      #Calculate the delta
      delta_name <- paste0(amount, '_delta_', vis, '_', vis + 1)
      df_result[[delta_name]] <- as.numeric(df_end[[amount]]) - as.numeric(df_start[[amount]])
    }
  }
  return(df_result[,-which(colnames(df_result) %in% amount_columns)])
}

#compute differences
df_feno_delta <- calculate_deltas(df_feno, exclude_col = 6:8)
clr_species_delta <- calculate_deltas(clr_species, exclude_col = 6:8)

for(col in colnames(df_feno_delta)[4:18]){
  b <- bestNormalize::orderNorm(df_feno_delta[[col]],warn = F)
  df_feno_delta[[col]] <- b$x.t
}


#create a complete dataset for species level
species_complete_delta <- merge(clr_species_delta,df_feno_delta, by = 'Code')
vect <- c('Code','Woman.x','Woman.y','Visit_number.x','Visit_number.y')
species_complete_delta <- species_complete_delta[,c(vect, setdiff(colnames(species_complete_delta),vect))]
species_complete_delta <- species_complete_delta[,-c(3,5)]
colnames(species_complete_delta)[2:3] <- c('Woman','Visit_number')


set.seed(123)  # Set seed for reproducibility
num_permutations <- 1000  # Number of permutations

#Initialize an empty list to store model results
results_list_species <- list()

#fit the models and store results
for (vis in 1:3) {
  for (horm in c('PROG', 'PRL', 'LH', 'FSH', 'BES17')) {
    for (bact in species_to_test) {
      # Define the delta variable for the current bacteria and visit
      bact_vis <- paste0(bact, '_delta_', vis, '_', vis + 1)
      horm_vis <- paste0(horm, '_delta_', vis, '_', vis + 1)
      
      # Create the formula for the model
      formula_str <- paste0(horm_vis, ' ~ ', bact_vis, ' + Visit_number + (1|Woman)')
      formula <- as.formula(formula_str)
      cat(as.character(formula))
      cat('\n')
      # Fit the model, handling potential errors or warnings gracefully
      model <- lmer(formula, data = species_complete_delta)
      
      summary_model <- summary(model)
      #Extract key information from the model
      real_p_value <- summary_model$coefficients[2, "Pr(>|t|)"]
      result <- data.frame(
        Hormone = horm,
        Bacteria = bact,
        Visit = paste0(vis, '_', vis + 1),
        Estimate = summary_model$coefficients[2, "Estimate"],
        Std_Error = summary_model$coefficients[2, "Std. Error"],
        p.value = real_p_value
      )
      if(real_p_value < 0.05){
        cat(crayon::bgYellow$bold$blue("<<<< Randomization here >>>>")) #display it's going to fitting randomized models
        cat('\n')        
        # Permutation test
        permuted_p_values <- numeric(num_permutations)
        for (i in 1:num_permutations) {
          #Permute the bacteria column
          species_complete_delta[[bact_vis]] <- sample(species_complete_delta[[bact_vis]])
          
          #Fit the model on the permuted data
          permuted_model <- tryCatch({
            lmer(formula, data = species_complete_delta)
          }, error = function(e) {
            NULL
          })
          
          #If model fitting is successful, store the permuted p-value
          if (!is.null(permuted_model)) {
            permuted_summary <- summary(permuted_model)
            permuted_p_values[i] <- permuted_summary$coefficients[2, "Pr(>|t|)"]
          } else {
            permuted_p_values[i] <- NA  # Handle errors by setting NA
          }
        }
        #Count how many permuted p-values are less than the real p-value
        permuted_p_values <- permuted_p_values[!is.na(permuted_p_values)]  # Remove NA values
        num_smaller_p <- sum(permuted_p_values <= real_p_value)
        perm_p_value <- num_smaller_p / length(permuted_p_values)
        
        #Add permutation p-value to the result
        result$p.permuted <- perm_p_value
      }else{
        result$p.permuted <- NA
      }
      
      
      #Append to the results list
      results_list_species[[length(results_list_species) + 1]] <- result
    }
  }
}

# Combine all the results into a single dataframe for species
results_df <- do.call(rbind, results_list_species)

#compute 95% CI
results_df$CI_Lower <- results_df$Estimate - 1.96 * results_df$Std_Error
results_df$CI_Upper <- results_df$Estimate + 1.96 * results_df$Std_Error

#create a label for tested associations
results_df$Label <- paste0(results_df$Hormone, " - ", results_df$Bacteria, " (Visit ", results_df$Visit, ")")

#order the results
results_df$Label <- factor(results_df$Label, levels = results_df$Label[order(results_df$Estimate)])

results_df$Label <- gsub("\\.", " ", results_df$Label)
#make the forest plot
g <- ggplot(subset(results_df, p.value < 0.05), aes(x = Estimate, y = Label)) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%.3f", p.value)),
            hjust = 0.5, vjust = -1.7, size = 1.5) +
  geom_text(aes(label = sprintf("%.3f", p.permuted)),
            hjust = 0.5, vjust = 1.7, size = 1.5,color = 'blue') +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Coefficient",
    y = ""
  ) +
  theme_minimal() + 
  theme(axis.text.y = element_text(face = 'italic',size = 5,family = 'sans'),
        axis.title.x = element_text(size = 5),
        axis.text.x = element_text(size = 5))

for(col in colnames(results_df)){
  if(is.numeric(results_df[[col]])){
    results_df[[col]] <- round(results_df[[col]],5)
  }
}
#save results
write.csv(x = results_df,file = 'results/linear_models/delta_models_results.csv')
png('results/linear_models/delta_models_results.png',height = 900,width = 900, res = 300)
plot(g)
dev.off()

pdf('results/linear_models/delta_models_results.pdf',width = 5,height = 3)
plot(g)
dev.off()


