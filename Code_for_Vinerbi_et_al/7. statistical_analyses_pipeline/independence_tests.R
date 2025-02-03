# Rscript
# This script run the indipendence tests between PCoAs and categorical variables
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 06/12/2024
# R version: R v4.4.1




library(ggplot2)
library(vegan)
library(plotly)
library(reshape2)
library(gridExtra)

#get the abundances
rel_abb_genus <- t(read.csv('data/Relative_abundances_genus.csv', sep = '\t'))
rel_abb_genus <- as.data.frame(rel_abb_genus)
rel_abb_genus$Code <- rownames(rel_abb_genus)
rel_abb_genus$Visit <- sapply(rel_abb_genus$Code,function(stringa){substring(stringa,6,6)})

rel_abb_species <- t(read.table('data/Relative_abundances_species.csv', sep = '\t'))
rel_abb_species <- as.data.frame(rel_abb_species)
rel_abb_species$Code <- rownames(rel_abb_species)
rel_abb_species$Visit <- sapply(rel_abb_species$Code,function(stringa){substring(stringa,6,6)})

#get the laboratory data
df_lab <- read.csv('data/Laboratory_data.csv')

#take from the command line whether use the linear correction of pcoa or not (T = TRUE, F = FALSE)
correction <- commandArgs(trailingOnly = TRUE) 
correction <- ifelse(correction == 'T', TRUE, ifelse(correction == 'F', FALSE, stop('Enter a valid parameter')))

#get the questionnaire
quest <- read.csv('data/clean_questionnaire.csv')

#create a basic function to compute the PCoA of relative abundances
pcoa_simple <- function(choice,visit = 'all'){
  
  if((choice == 'genera') || (choice == 'genus')){
    beta <- vegdist(rel_abb_genus[,-c(ncol(rel_abb_genus)-1,ncol(rel_abb_genus))])
    df <- as.data.frame(cmdscale(beta,k = 5))
  }else if((choice == 'species')){
    beta <- vegdist(rel_abb_species[,-c(ncol(rel_abb_species)-1,ncol(rel_abb_species))])
    df <- as.data.frame(cmdscale(beta,k = 5))
  }
  if(visit != 'all'){ #to filter per visit
    df$Visit <- sapply(rownames(df),function(stringa){substring(stringa,6,6)})
    df <- subset(df, Visit == visit)
    df <- df[,-6]
  }
  colnames(df) <- c('PCoA1','PCoA2','PCoA3','PCoA4','PCoA5')

  return(df)
  
} #function the return the pcoa corresponding to choice

#analize the correlations between PCoAs and techical variables
results_df <- data.frame( #to store the results of correlations
  PCoA = character(),
  Comparison = character(),
  Estimate = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

#calculate the simple pcoa (without correction)
df_lab_pcoa_gen <- pcoa_simple(choice = 'genus',visit = 'all')
df_lab_pcoa_gen$Code <- rownames(df_lab_pcoa_gen)
df_lab_pcoa_gen <- merge(df_lab_pcoa_gen,df_lab,by = 'Code')

df_lab_pcoa_spe <- pcoa_simple(choice = 'species',visit = 'all')
df_lab_pcoa_spe$Code <- rownames(df_lab_pcoa_spe)
df_lab_pcoa_spe <- merge(df_lab_pcoa_spe,df_lab,by = 'Code')

#correlations analyses
for(i in 1:5){
  name_pc <- paste('PCoA', i, sep = '')
  
  test1 <- cor.test(df_lab_pcoa_gen[[name_pc]], df_lab_pcoa_gen$Qubit_DNA, method = 'spearman',exact = F)
  results_df <- rbind(results_df, data.frame(
    PCoA = name_pc,
    Comparison = 'Qubit_DNA vs genus PCoA',
    Estimate = test1$estimate,
    P_Value = test1$p.value
  ))
  test2 <- cor.test(df_lab_pcoa_gen[[name_pc]], df_lab_pcoa_gen$Qubit_Library, method = 'spearman',exact = F)
  results_df <- rbind(results_df, data.frame(
    PCoA = name_pc,
    Comparison = 'Qubit_Library vs genus PCoA',
    Estimate = test2$estimate,
    P_Value = test2$p.value
  ))
  test3 <- cor.test(df_lab_pcoa_gen[[name_pc]], as.numeric(df_lab_pcoa_gen$Total_counts), method = 'spearman',exact = F)
  results_df <- rbind(results_df, data.frame(
    PCoA = name_pc,
    Comparison = 'Total_counts vs genus PCoA',
    Estimate = test3$estimate,
    P_Value = test3$p.value
  ))
  
  test4 <- cor.test(df_lab_pcoa_spe[[name_pc]], as.numeric(df_lab_pcoa_spe$Qubit_DNA), method = 'spearman',exact = F)
  results_df <- rbind(results_df, data.frame(
    PCoA = name_pc,
    Comparison = 'Qubit_DNA vs species PCoA',
    Estimate = test4$estimate,
    P_Value = test4$p.value
  ))
  test5 <- cor.test(df_lab_pcoa_spe[[name_pc]], as.numeric(df_lab_pcoa_spe$Qubit_Library), method = 'spearman',exact = F)
  results_df <- rbind(results_df, data.frame(
    PCoA = name_pc,
    Comparison = 'Qubit_Library vs species PCoA',
    Estimate = test5$estimate,
    P_Value = test5$p.value
  ))
  test6 <- cor.test(df_lab_pcoa_spe[[name_pc]], as.numeric(df_lab_pcoa_spe$Total_counts), method = 'spearman',exact = F)
  results_df <- rbind(results_df, data.frame(
    PCoA = name_pc,
    Comparison = 'Total_counts vs species PCoA',
    Estimate = test6$estimate,
    P_Value = test6$p.value
  ))
}

rownames(results_df) <- NULL

results_df$P_Value <- sapply(results_df$P_Value, function(p){ifelse(p*20 <=1, p*20, 1)})

#create a function that correct linearly the PCoAs from pcoa_simple using results_df calculated above
pcoa <- function(correct = correction, choice, visit = 'all') {
  df <- pcoa_simple(choice = choice, visit = visit)
  
  if (correct) {
    correction_to_make <- subset(results_df, P_Value < 0.05)
    df$Code <- rownames(df)
    dummy_df <- merge(df, df_lab, by = 'Code')
    for (i in 1:5) {
      check <- subset(correction_to_make, PCoA == paste0('PCoA', i))
      tec_lab_sig <- unique(sapply(check$Comparison, function(string) {
        ifelse(unlist(strsplit(string, " "))[3] == choice, unlist(strsplit(string, " "))[1], NA)
      }))
      tec_lab_sig <- tec_lab_sig[!is.na(tec_lab_sig)]
      if (length(tec_lab_sig) != 0) {
        formulaa <- as.formula(paste(paste('PCoA', i, sep = ''), ' ~ ', paste(tec_lab_sig, collapse = ' + ')))
        print(formulaa)
        model <- lm(data = dummy_df, formula = formulaa)
        df[[paste('PCoA', i, sep = '')]] <- residuals(model)
      } else {
        print(paste("There's nothing to correct for ", i, '-th component', sep = ''))
      }
    }
    
    df <- df[, -6]  # Remove the 6th column if needed
  }
  
  return(df)
}

##################
#Independence test
##################

#independence test function basic 
independence_test <- function(col,choice,visit = 1,plot = T,save_plot = F,print_plot = T){

  pcoa_specific <- pcoa(choice = choice,visit = visit)
  values <- setdiff(unique(quest[,col]),c(NA)) #take the unique values of col (categorical)
  values <- sapply(values, as.character) #make them character
  comparison_vector <- list() #initialize a list
  
  for(val in values){ #build a list of dataframes: each of them will be the pcoa data with fixed val
    code <- subset(quest,quest[,col] == val)$Code #grab the codes of this women
    rows <- intersect(code,rownames(pcoa_specific)) #grab the rows of pcoa cointaining those codes
    comparison_vector[[val]] <- pcoa_specific[rows,] #grab the pcoa subset cointaing those rows
  }
  
  pvalues <- NULL #here we will store the pvalues of the tests between pcoa and col
  for(i in 1:5){ #I got i varying across 1,2,3,4,5
    xs <- list() #empty list (it will contain the pcoa relative to val with pcoa fixed)
    for(val in values){
      xs[[val]] <- comparison_vector[[val]][,i] #Take the i-th column of the pcoa relative to val
      xs[[val]] <- t(xs[[val]]) #translate it
    } #I select the columns corresponding to the i-th component and the j-th value of values
    group <- NULL 
    for(val in values){ 
      v_rep <- rep(val, nrow(as.data.frame(comparison_vector[val]))) 
      group <- c(group,v_rep)
    } #create a vector to identify which val the pcoa belong to (group)
    data <- data.frame(value = as.vector(do.call(cbind,xs)),group = group) #create a long dataframe
    #containing all the value of i-th pcoa with specific value of col
    test <- kruskal.test(value ~ group,data = data) #Run the test in order to see if there's significant
    #difference
    pvalues[i] <- test$p.value #save the p value of the test
  }
  names(pvalues) <- colnames(pcoa_specific) #give names to pvalues vector
  df <- pcoa_specific #take pcoa relative to choice (genus or species)
  df$Code <- rownames(pcoa_specific) #define the code column
  df_1 <- quest[,c(col,'Code')] #take the part of quest relative to col and Code
  df_1[,col] <- sapply(df_1[,col],as.factor) #make sure that col is interpreted as factor
  df <- merge(df,df_1,by='Code') #merging the two dataframe
  if(plot){ #check if the user wanna plot the pcoa
    g <- ggplot() + geom_point(aes(df$PCoA1,df$PCoA2,color = df[,col]),alpha = 0.5,size = 3) + 
      xlab('PCoA1') + ylab('PCoA2') + scale_color_discrete(col) + theme_minimal() + 
      theme(legend.position = "bottom") 
    if(print_plot){ #check if the user wanna actually see the plot of pcoa (1,2)
      print(g)
    }
    if(save_plot){ #check if the user wanna save the plot
      if(visit == 'all'){visit_name = 'overall'}else{visit_name = visit}
      name_file <- paste0('independence_tests_',col,'_visit_',visit_name,'.jpg') #name of the file
      path_folder <- paste0('results/independence_test/',choice) #path of the folder
      path_folder <- ifelse(correction,paste0(path_folder,'/corrected/'),
                            paste0(path_folder,'/not_corrected/'))
      if(!dir.exists(path_folder)){ #check if the folder exists
        dir.create(path_folder) #if not exists: create it
      }
      complete_path <- paste(path_folder,name_file, sep = '') #complete path of the file
      ggsave(complete_path,g,dpi = 1000,width = 10) #save the plot
    }
  }
  
  return(pvalues)
} #basic function to make one single independence test

#independence test for a vector of categorical data
independence_test_vect <- function(to_check,choice,visit = 1,save_plots = T,print_plot = F){
  list_pvalues <- list()
  for(col in to_check){
    pvalues <- independence_test(col,choice,visit = visit,plot = T,save_plot = save_plots,print_plot = F)
    list_pvalues[[col]] <- pvalues
  }
  complete_results <- as.data.frame(do.call(rbind,list_pvalues))
  complete_results$Visit_number <- visit
  return(complete_results)
}
  
#create a function that check if there are any correlations with the PCoAs beside the first two and create the relative scatterplot
check_for_extra_graph <- function(results,choice,visit = 1,save_plot = F) { #it wants the table of all cols we've tested
  
  # Intersect the codes and subset the data
  rows <- intersect(quest$Code, rownames(pcoa(choice = choice, visit = visit)))
  quest_inter <- subset(quest, Code %in% rows)
  df_pcoa_inter <- pcoa(choice = choice,visit = visit)[rows, ]
  
  list_plots <- list()
  
  for (i in 1:3) { # Loop over PCoA3, PCoA4, PCoA5 (we've already plot the first two)
    
    pvalues_col <- results[, i + 2] # P-values for the i-th PCoA (adjusted index) (vector)
    names(pvalues_col) <- rownames(results) #names of the cols
    sig_col <- names(pvalues_col[pvalues_col < 0.05]) # Get significant variables (their names)
  
    if (length(sig_col) > 0) { # If there are significant p-values
      
      list_plots_cat <- lapply(sig_col, function(cat) {
        
        name_y <- paste('PCoA', i + 2, sep = '')
        colore <- quest_inter[, cat]
        
        # Create a combined data frame for plotting
        plot_data <- data.frame(PCoA1 = df_pcoa_inter$PCoA1, 
                                PCoA_Y = df_pcoa_inter[, i + 2], 
                                colore = factor(colore))
        
        g <- ggplot(plot_data, aes(x = PCoA1, y = PCoA_Y, color = colore)) + 
          geom_point(alpha = 0.5,size = 3) +
          xlab('PCoA1') + 
          ylab(name_y) + 
          labs(color = cat) + 
          theme_minimal() + 
          theme(legend.position = "bottom")
        visit_name <- visit
        name_file <- paste('independence_tests_',name_y,'_',cat,'_visit',visit_name,'.jpg',sep = '') #name of the file
        path_folder <- paste0('results/independence_test/',
                             choice) #path of the folder
        path_folder <- ifelse(correction,paste0(path_folder,'/corrected/'),
                              paste0(path_folder,'/not_corrected/'))
        complete_path <- paste0(path_folder,name_file) #complete path of the file
        if(save_plot){
          ggsave(complete_path,g,dpi = 1000,width = 10) #save the plot
        }
        return(g)
      })
      names(list_plots_cat) <- sig_col
      list_plots[[i]] <- list_plots_cat
    } else {
      list_plots[[i]] <- 'no sig'
    }
  }
  names(list_plots) <- c('PCoA3', 'PCoA4', 'PCoA5')
  return(list_plots)
}

#create a funcion that create a sintetic plot to visualize the results of independence test on a vector of categorical variables
create_plot_independence <- function(to_check, choice, visit = 1, save_plot = F) {
  results <- independence_test_vect(to_check = to_check, choice = choice, visit = visit, save_plots = F, print_plot = F)
  name_plot <- paste('Independence tests within visit', visit, choice, sep = ' ')
  n_rows <- nrow(results)
  long_results <- data.frame(
    p_values = c(results$PCoA1, results$PCoA2, results$PCoA3, results$PCoA4, results$PCoA5),
    PCoA = c(rep(1, n_rows), rep(2, n_rows), rep(3, n_rows), rep(4, n_rows), rep(5, n_rows)),
    feature = rep(rownames(results), 5)
  )
  long_results$visit <- visit
  long_results$log <- sapply(long_results$p_values, function(x) { return(-log10(x)) })
  long_results$color <- sapply(long_results$p_values, function(p) {
    if (p < 0.01) {
      return('p_value < 0.01')
    } else if ((p < 0.05) & (p >= 0.01)) {
      return('0.01 < p_value < 0.05')
    } else {
      return('p_value > 0.05')
    }
  })
  # Creazione del grafico con geom_raster
  g <- ggplot(long_results, aes(feature, PCoA, fill = color)) +
    geom_tile(color = "black") +  # Set the border color to black
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 10),
          legend.position = 'none') +
    xlab('') + 
    scale_fill_manual(values = c('p_value < 0.01' = 'red', '0.01 < p_value < 0.05' = "pink", 'p_value > 0.05' = "white"),  # Define the colors manually
                      name = "Significance Level")
  # Mostra il grafico
  print(g)
  
  if (save_plot) {
    visit_name <- ifelse(visit == 'all', 'overall', visit)
    name_file <- paste0('complete_independence_tests_', choice, '_visit_', visit_name, '.jpg')  # Name of the file
    path_folder <- paste0('results/independence_test/', choice)  # Path of the folder
    path_folder <- ifelse(correction,paste0(path_folder,'/corrected/'),
                          paste0(path_folder,'/not_corrected/'))
    complete_path <- paste0(path_folder, name_file)  # Complete path of the file
    ggsave(complete_path, g, dpi = 1000, width = 10)
  }
  return(list(g, long_results))
}

#function that wrap the functions above
make_analysis <- function(to_check,choice = c('genus','species'),visit = 1,save_results = T,save_plots = F){
  list_of_grid <- list()
  for(cl in choice){
    list_of_plots <- list()
    results <- independence_test_vect(to_check = to_check,choice = cl,visit = visit,save_plots = save_plots)
    check_for_extra_graph(results = results,choice = cl,save_plot = T)
    create_plot_independence(to_check = to_check,choice = cl,visit = visit,save_plot = save_plots)
    if(save_results){
      visit_name <- visit
      name_file <- paste('results_independence_tests_',cl,'_visit_',visit_name,'.csv',sep = '') #name of the file
      path_folder <- paste0('results/independence_test/',cl) #path of the folder
      path_folder <- ifelse(correction,paste0(path_folder,'/corrected/'),
                            paste0(path_folder,'/not_corrected/'))
      complete_path <- paste0(path_folder,name_file) #complete path of the file
      write.csv(results,complete_path)
    }
  }
}

to_check <- c('Age_category','Smoking_habits','Swab_morning_or_not','Swab_after_feces',
              'sex_less_than_two_days','Pregnancy_category','BMI_ranges',
              'WHR_ranges','Pill_use',
              'Med_or_supp_check','Bristol_stool_scale')

####multi_panel_genus####
list_of_results_genus <- list()
for(v in 1:4){ #analysis per each visit
  results <- independence_test_vect(to_check = to_check, choice = 'genus', visit = v, save_plots = F, print_plot = F)
  results$feature <- rownames(results)
  list_of_results_genus[[v]] <- results
}

all_results_genus <- do.call(rbind, list_of_results_genus) #wrap results per each visit
rownames(all_results_genus) <- NULL
all_results_genus_long <- melt(all_results_genus, id.vars = c("Visit_number", "feature"), 
                               measure.vars = c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5"),
                               variable.name = "PCoA", value.name = "value")
all_results_genus_long$color <- sapply(all_results_genus_long$value, function(p) { #fix the color based on significances
  if (p < 0.01) {
    return('p_value < 0.01')
  } else if ((p < 0.05) & (p >= 0.01)) {
    return('0.01 < p_value < 0.05')
  } else {
    return('p_value > 0.05')
  }
})

all_results_genus_long$PCoA <- sapply(all_results_genus_long$PCoA, function(string){substring(string,5,5)}) 

# Create a named vector for the custom labels
custom_labels <- c(
  "1" = "F",
  "2" = "O",
  "3" = "EL",
  "4" = "LL"
)

genus_grid <- ggplot(all_results_genus_long, aes(PCoA, feature, fill = color)) +
  geom_tile(color = "black") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 14),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.grid.major = element_line(color = NA),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(color = 'black'),
        panel.background = element_rect(color = 'black')) +
  scale_fill_manual(values = c(
    'p_value < 0.01' = '#ffa32b',
    '0.01 < p_value < 0.05' = '#f0dc6c',
    'p_value > 0.05' = "white"
  ), name = "Significance Level")  +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab('PCoA') +  
  ylab('Feature') + 
  facet_wrap(~ Visit_number, nrow = 1, labeller = labeller(Visit_number = custom_labels)) #create the multipanel for genus

####multi_panel_species####

list_of_results_species <- list()
for(v in 1:4){
  results <- independence_test_vect(to_check = to_check, choice = 'species', visit = v, save_plots = F, print_plot = F)
  results$feature <- rownames(results)
  list_of_results_species[[v]] <- results
}

all_results_species <- do.call(rbind, list_of_results_species)
rownames(all_results_species) <- NULL
all_results_species_long <- melt(all_results_species, id.vars = c("Visit_number", "feature"), 
                               measure.vars = c("PCoA1", "PCoA2", "PCoA3", "PCoA4", "PCoA5"),
                               variable.name = "PCoA", value.name = "value")

all_results_species_long$color <- sapply(all_results_species_long$value, function(p) {
  if (p < 0.01) {
    return('p_value < 0.01')
  } else if ((p < 0.05) & (p >= 0.01)) {
    return('0.01 < p_value < 0.05')
  } else {
    return('p_value > 0.05')
  }
})
all_results_species_long$PCoA <- sapply(all_results_species_long$PCoA, function(string){substring(string,5,5)})

species_grid <- ggplot(all_results_species_long, aes(PCoA,feature, fill = color)) + geom_tile(color = "black") + 
  theme_bw() + theme(axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 14),
                     legend.position = 'right',
                     legend.title = element_text(size = 14),                     
                     legend.text = element_text(size = 14),
                     panel.grid.major = element_line(color = NA),
                     strip.text = element_text(face = "bold"),
                     strip.background = element_rect(color = 'black'),
                     panel.background = element_rect(color = 'black')) + # Rimuove i tick interni
  scale_fill_manual(values = c(
    'p_value < 0.01' = '#ffa32b',
    '0.01 < p_value < 0.05' = '#f0dc6c',
    'p_value > 0.05' = "white"
  ), name = "Significance Level") +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  xlab('PCoA') +  
  ylab('Feature') + 
  facet_wrap(~ Visit_number, nrow = 1,labeller = labeller(Visit_number = custom_labels)) #create the multipanel for species

#### multi-panel genus and species ####
if(correction){
  # Add a "type" column to the genus and species datasets
  all_results_genus_long$type <- "Genus"
  all_results_species_long$type <- "Species"
  
  # Combine the genus and species datasets into a single dataframe
  combined_results_long <- rbind(all_results_genus_long, all_results_species_long)
  
  # Create the combined plot with genus and species as separate panels
  combined_grid <- ggplot(combined_results_long, aes(PCoA, feature, fill = color)) + 
    geom_tile(color = "black") + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 14),
          legend.position = 'right',
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          panel.grid.major = element_line(color = NA),
          strip.text = element_text(face = "bold", size = 14),
          strip.background = element_rect(fill = "grey90", color = 'black'),
          panel.background = element_rect(color = 'black')) +
    
    # Customize fill colors based on p-value significance levels
    scale_fill_manual(values = c(
      'p_value < 0.01' = '#ffa32b',
      '0.01 < p_value < 0.05' = '#f0dc6c',
      'p_value > 0.05' = "white"
    ), name = "Significance Level") +
    
    # Remove extra spacing in x and y axes
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) +
    
    # Axis labels
    xlab('PCoA') +  
    ylab('Feature') +  
    
    # Create separate rows for "Genus" and "Species" with columns for each visit
    facet_grid(type ~ Visit_number, scales = "free_y",labeller = labeller(Visit_number = custom_labels))
  
  ggsave('results/independence_test/multipanel_genus_species.jpg', combined_grid, width = 15,height = 10)
}


####save genus and species grid####
path_folder <- 'results/independence_test/genus'
path_folder <- ifelse(correction, paste0(path_folder,'/corrected/'),paste0(path_folder,'/not_corrected/'))
ggsave(paste0(path_folder,'grid_genus.jpg'), genus_grid, width = 15,height = 10)

path_folder <- 'results/independence_test/species'
path_folder <- ifelse(correction, paste0(path_folder,'/corrected/'),paste0(path_folder,'/not_corrected/'))
ggsave(paste0(path_folder,'grid_species.jpg'), species_grid, width = 15,height = 10)


for(v in 1:4){
  make_analysis(to_check = to_check,save_results = T,visit = v,choice = 'genus',save_plots = F)
  make_analysis(to_check = to_check,save_results = T,visit = v,choice = 'species',save_plots = F)
}

####EXTRACT FINAL SIGNIFICANT VARIABLES####

df_list <- list()
for(choice in c('genus','species')){
  visit_limited_list <- list()
  for(vis in 1:4){
    path <- file.path('results/independence_test',
                      choice,'corrected',
                      paste0('results_independence_tests_', choice, paste0('_visit_', vis), '.csv'))
    df <- read.csv(path) #read the results for each visit created above
    colnames(df)[1] <- 'Variable'
    visit_limited_list[[vis]] <- df
  }
  df_list[[choice]] <- visit_limited_list
}


variables <- df_list$genus[[1]]$Variable
significances_per_class <- list()
for(choice in c('genus','species')){
  significances <- list()
  for(visit in 1:4){
    df <- df_list[[choice]][[visit]]
    significances_per_visit <- list()
    for(var in variables){
      vect <- df[which(df$Variable == var), ]
      vect <- which(vect[2:6] < 0.05)
      if(length(vect) != 0) significances_per_visit[[var]] <- paste0('PCoA',vect)
    }
    significances[[paste0('visit_',visit)]] <- significances_per_visit
  }
  significances_per_class[[choice]] <- significances
}


significant_result_genus <- Filter(function(x) length(x) > 0, significances_per_class$genus)
significant_result_genus <- as.data.frame(melt(significant_result_genus))
colnames(significant_result_genus) <- c('PCoA','Variable','Visit')
final_variables_genus <- NULL
i <- 1
for(var in unique(significant_result_genus$Variable)){
  to_check <- subset(significant_result_genus, Variable == var)
  if(length(intersect(c('PCoA1', 'PCoA2'),to_check$PCoA)) != 0){
    final_variables_genus[i] <- var
    i <- i + 1
  }else{
    if(length(unique(to_check$Visit)) != 1){
      final_variables_genus[i] <- var
      i <- i + 1
    }
  }
}
significant_result_genus$Classification <- 'genus'

significant_result_species <- Filter(function(x) length(x) > 0, significances_per_class$species)
significant_result_species <- as.data.frame(melt(significant_result_species))
colnames(significant_result_species) <- c('PCoA','Variable','Visit')
final_variables_species <- NULL
i <- 1
for(var in unique(significant_result_species$Variable)){
  to_check <- subset(significant_result_species, Variable == var)
  if(length(intersect(c('PCoA1', 'PCoA2'),to_check$PCoA)) != 0){
    final_variables_species[i] <- var
    i <- i + 1
  }else{
    if(length(unique(to_check$Visit)) != 1){
      final_variables_species[i] <- var
      i <- i + 1
    }
  }
}
significant_result_species$Classification <- 'species'

significant_result <- rbind(significant_result_genus,significant_result_species)

write.csv(significant_result, 'results/independence_test/significant_variable_independence_test.csv')

final_variables <- union(final_variables_genus,final_variables_species)

writeLines(final_variables,'data/independence_significant_variables.txt') #save in a txt file the significant variables resulting from IT

####COMPARISON BETWEEN CORRECTED AND NOT CORRECTED PCoAs####

g1 <- ggplot(pcoa(choice = 'genus',visit = 'all',correct = F)) + geom_point(aes(PCoA1,PCoA2)) + theme_bw()
g2 <- ggplot(pcoa(choice = 'genus',visit = 'all',correct = T)) + geom_point(aes(PCoA1,PCoA2)) + theme_bw()


s1 <- ggplot(pcoa(choice = 'species',visit = 'all',correct = F)) + geom_point(aes(PCoA1,PCoA2)) + theme_bw()
s2 <- ggplot(pcoa(choice = 'species',visit = 'all',correct = T)) + geom_point(aes(PCoA1,PCoA2)) + theme_bw()

#save the comparison
ggsave('results/comparison_genus_correct_not_correct.jpg', grid.arrange(g1,g2))
ggsave('results/comparison_species_correct_not_correct.jpg', grid.arrange(s1,s2))

#save the results of correlations between pcoa and technical variables
if(correction){
  write.csv(subset(results_df,P_Value < 0.05),'results/results_correlations_pcoa_technical_variables.csv')
  
  res_genus_1 <- read.csv('results/independence_test/genus/corrected/results_independence_tests_genus_visit_1.csv')
  res_genus_2 <- read.csv('results/independence_test/genus/corrected/results_independence_tests_genus_visit_2.csv')
  res_genus_3 <- read.csv('results/independence_test/genus/corrected/results_independence_tests_genus_visit_3.csv')
  res_genus_4 <- read.csv('results/independence_test/genus/corrected/results_independence_tests_genus_visit_4.csv')

  
  res_species_1 <- read.csv('results/independence_test/species/corrected/results_independence_tests_species_visit_1.csv')
  res_species_2 <- read.csv('results/independence_test/species/corrected/results_independence_tests_species_visit_2.csv')
  res_species_3 <- read.csv('results/independence_test/species/corrected/results_independence_tests_species_visit_3.csv')
  res_species_4 <- read.csv('results/independence_test/species/corrected/results_independence_tests_species_visit_4.csv')

  #wrap all results
  res_genus <- rbind(res_genus_1,res_genus_2,res_genus_3,res_genus_4)
  res_genus$Level <- 'Genus'
  res_species <- rbind(res_species_1,res_species_2,res_species_3,res_species_4)
  res_species$Level <- 'Species'
  
  res_complete <- rbind(res_genus,res_species)
  colnames(res_complete)[1] <- 'Variable'
  
  for(col in c('PCoA1','PCoA2','PCoA3','PCoA4','PCoA5')){
    res_complete[[col]] <- round(res_complete[[col]],5)
  }
  
  write.csv(x = res_complete,file = 'results/independence_test/complete_res_independence_tests.csv')
  
}






