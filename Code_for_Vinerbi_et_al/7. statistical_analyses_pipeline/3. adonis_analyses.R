# Rscript
# This script contains several functions useful to run adonis analysis with visualization
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 03/02/2025
# R version: R v4.4.1


library(vegan)
library(ggplot2)
library(stringr)
library(tidyverse)
library(reshape2)

#get and prepare the data
rel_abb_genus <- read.table('data/Relative_abundances_genus.csv', sep = '\t')
rel_abb_genus <- as.data.frame(t(rel_abb_genus))
rel_abb_genus$Code <- rownames(rel_abb_genus)

rel_abb_species <- read.table('data/Relative_abundances_species.csv', sep = '\t')
rel_abb_species <- as.data.frame(t(rel_abb_species))
rel_abb_species$Code <- rownames(rel_abb_species)

#find the dominant lactobacillus
to_find_abundance <- function(row) {
  most_abundant <- names(row)[which.max(row)]
  if (most_abundant %in% c("Lactobacillus iners", 
                           "Lactobacillus crispatus", 
                           "Lactobacillus gasseri")) {
    return(most_abundant)
  }
  return("Other")
}
most_abundant <- apply(rel_abb_species[,-ncol(rel_abb_species)], 
                                     1, to_find_abundance)
dummy_data <- as.data.frame(most_abundant)
dummy_data$Code <- rownames(dummy_data)
rownames(dummy_data) <- NULL

#get all metadata
quest <- read.csv('data/clean_questionnaire.csv')
df_feno <- read.csv('data/clean_phenotypes.csv')
df_feno <- df_feno[,c('PROG','LH','FSH','BES17','PRL','HDL','LDL','COL','TRI','AST','ALT','Code')]

df_genus <- merge(rel_abb_genus,quest,by = 'Code')
df_genus <- merge(df_genus,df_feno,by = 'Code')
df_genus$Visit_number <- as.factor(df_genus$Visit_number)

df_species <- merge(rel_abb_species,quest,by = 'Code')
df_species <- merge(df_species,df_feno,by = 'Code')
df_species$Visit_number <- as.factor(df_species$Visit_number)


#function to make adonis analysis 
do_analysis <- function(col,covariates = NULL, choice, visit,plot = F,cores = 1,rep = F,pcoa_to_plot = 1:2,
                        discretize = T){
  
  #take the piece of data relative to choice (genus or species)
  if(choice %in% c('genus','genera')){
    df <- subset(df_genus, Visit_number %in% visit) 
    to_remove <- setdiff(colnames(df_genus),colnames(rel_abb_genus))
    to_remove <- which(colnames(df_genus) %in% to_remove)
    to_remove <- c(1, to_remove)
  } else {
    df <- subset(df_species, Visit_number %in% visit)
    to_remove <- setdiff(colnames(df_species),colnames(rel_abb_species))
    to_remove <- which(colnames(df_species) %in% to_remove)
    to_remove <- c(1, to_remove)
  }
  
  all_variables <- c(col, covariates)
  df_not_na <- df[complete.cases(df[, all_variables]), ]
  #discretize the continuous data
  if ((is.numeric(df_not_na[[col]]))&(discretize)) {
    # Calculate quartile breakpoints
    quartile_breaks <<- quantile(df_not_na[[col]], probs = seq(0, 1, 0.5))
    categorized_data <- cut(df_not_na[[col]], breaks = quartile_breaks, 
                            include.lowest = TRUE,
                            labels = c(paste(quartile_breaks[[1]],'<', col, '<',quartile_breaks[[2]]), 
                                       paste(quartile_breaks[[2]],'<', col, '<',quartile_breaks[[3]])))
    df_not_na[[col]] <- categorized_data
    df_not_na[[col]] <- as.factor(df_not_na[[col]])
  }else if(!(is.numeric(df_not_na[[col]]))){
    df_not_na[[col]] <- as.factor(df_not_na[[col]])
  }
  
  # Perform the adonis2 analysis
  if(length(covariates) != 0){
    right_member <- paste(c(col,covariates),collapse = ' + ')
    formula <- as.formula(paste("df_not_na[, -to_remove] ~", right_member))
  }else{
    formula <- as.formula(paste("df_not_na[, -to_remove] ~", col))
  }
  
  
  if(rep){
    if(length(visit) != 4){
      print('If you want to take account of repetitions then you must give visit = 1:4 as input')
      return(NA)
    }else{
      res_adonis <- adonis2(formula, data = df_not_na, method = 'bray',parallel = cores,
                            strata = df_not_na$Visit_number)
    }
  }else{
    res_adonis <- adonis2(formula, data = df_not_na, method = 'bray',parallel = cores)
  }
  if(plot){
    
    # Perform PCoA analysis
    pcoa <- cmdscale(vegdist(df[, -to_remove], method = 'bray'), k = 5,eig = TRUE)
    eigenvalues <- pcoa$eig
    eigenvalues <- eigenvalues[eigenvalues > 0]
    eigenvalues <- eigenvalues/sum(eigenvalues) * 100
    eigenvalues <- eigenvalues[1:5]
    names(eigenvalues) <- c('PCoA1','PCoA2','PCoA3','PCoA4','PCoA5')
    pcoa <- as.data.frame(pcoa$points)
    pcoa$Code <- df$Code
    pcoa <- merge(pcoa,dummy_data, by = 'Code')

    # Combine PCoA results with the specified column
    d_set <- cbind(pcoa, df[[col]])
    colnames(d_set)[c(2:6,8)] <- c('PCoA1', 'PCoA2','PCoA3','PCoA4','PCoA5', col)
    
    # Calculate centroids for the selected PCoA axes
    centroid <- na.omit(d_set) %>%
      group_by(!!sym(col)) %>%
      summarize(across(starts_with("PCoA"), mean))
    
    # Set the axes for the plot
    pcoa_to_plot_1 <- paste('PCoA', pcoa_to_plot[1], sep = '')
    pcoa_to_plot_1_lab <- paste0(pcoa_to_plot_1,' (',round(eigenvalues[pcoa_to_plot_1],2),'%)')
    
    pcoa_to_plot_2 <- paste('PCoA', pcoa_to_plot[2], sep = '')
    pcoa_to_plot_2_lab <- paste0(pcoa_to_plot_2,' (',round(eigenvalues[pcoa_to_plot_2],2),'%)')   
    
    to_print <- paste('p-value = ', format(res_adonis$`Pr(>F)`[1], digits = 3), sep = '')
    
    statistics_tab <- table(d_set[[col]])
    statistics_string <- paste(names(statistics_tab), statistics_tab, sep = ": ", collapse = "\n ")
    # Create the plot
    g1 <- ggplot(d_set, aes(x = .data[[pcoa_to_plot_1]], y = .data[[pcoa_to_plot_2]], color = .data[[col]])) +
      geom_point() +
      geom_point(data = centroid, aes(x = .data[[pcoa_to_plot_1]], y = .data[[pcoa_to_plot_2]], color = .data[[col]]), 
                 shape = 15, size = 6) +
      stat_ellipse(data = d_set[!is.na(d_set[[col]]), ], type = 'norm') + 
      annotate(geom = 'text', x = Inf, y = Inf, label = to_print, hjust = 1.1, vjust = 1.1) +
      annotate(geom = 'text', x = -Inf, y = Inf, label = statistics_string, hjust = -0.05, vjust = 1.1) +
      theme(legend.position = 'bottom', panel.background = element_blank(),
            panel.grid = element_line(colour = 'lightgrey')) +
      labs(x = pcoa_to_plot_1_lab, y = pcoa_to_plot_2_lab)
    
    
    ellipse_labels <- d_set %>%
      group_by(most_abundant) %>%
      summarize(center_x = mean(!!sym(pcoa_to_plot_1)), 
                center_y = mean(!!sym(pcoa_to_plot_2)), 
                .groups = 'drop')
    
    colors <- c('Lactobacillus iners' = '#A6CEE3','Lactobacillus crispatus' = '#1F78B4','Lactobacillus gasseri' = '#B2DF8A',
                'Other' = '#559E3E')
    
    g2 <- ggplot(data = d_set, aes(x = !!sym(pcoa_to_plot_1), 
                                   y = !!sym(pcoa_to_plot_2), 
                                   color = most_abundant)) + geom_point() + 
      stat_ellipse(aes(color = most_abundant), alpha = 0.5, level = 0.95,
                   type = 'norm', linewidth = 0.6) + 
      labs(x = pcoa_to_plot_1_lab,
           y = pcoa_to_plot_2_lab) +
      theme_minimal() +
      theme(legend.position = 'none') + 
      geom_text(data = ellipse_labels, aes(x = center_x, y = center_y, 
                label = most_abundant), 
                color = 'black', 
                size = 3,
                vjust = -3)  + # Adjust size and font as needed
    scale_color_manual(values = colors)
    
    if(choice == 'species'){
      g <- gridExtra::grid.arrange(g2,g1, nrow = 1)
      print(g)
      grid::grid.text("A", x = 0.01, y = 0.99, just = c("left", "top"), gp = grid::gpar(fontsize = 15))
      grid::grid.text("B", x = 0.51, y = 0.99, just = c("left", "top"), gp = grid::gpar(fontsize = 15))
    }else{
      print(g1)
    }

  }
  return(res_adonis)
}


#adonis analyses between visit
do_visit_adonis <- function(visit_to_compare,choice,cores = 1,plot = T){
  if(choice %in% c('genus','genera')){
    
    df <- subset(df_genus, Visit_number %in% visit_to_compare)
    to_remove <- setdiff(colnames(df_genus),colnames(rel_abb_genus))
    to_remove <- which(colnames(df_genus) %in% to_remove)
    to_remove <- c(1, to_remove)
    
  } else {
    df <- subset(df_species, Visit_number %in% visit_to_compare)
    to_remove <- setdiff(colnames(df_species),colnames(rel_abb_species))
    to_remove <- which(colnames(df_species) %in% to_remove)
    to_remove <- c(1, to_remove)
  }
  formula <- as.formula(paste("df[, -to_remove] ~", 'Visit_number'))
  res_adonis <- adonis2(formula , data = df, method = 'bray',parallel = cores)

  if(plot){
    # Perform PCoA analysis
    pcoa <- as.data.frame(cmdscale(vegdist(df[, -to_remove], method = 'bray'), k = 2))

    # Combine PCoA results with the specified column
    d_set <- cbind(pcoa, df[, 'Visit_number'])
    colnames(d_set) <- c('PCoA1', 'PCoA2', 'Visit_number')

    # Calculate centroids for each group in the specified column
    centroid <- d_set %>%
      group_by(Visit_number) %>%
      summarize(PCoA1 = mean(PCoA1),
                PCoA2 = mean(PCoA2))
    centroid <- as.data.frame(centroid)
    
    to_print <- paste('p-value = ', format(res_adonis$`Pr(>F)`[1], digits = 3), sep = '')
    
    g <- ggplot(d_set, aes(PCoA1, PCoA2, color = Visit_number)) +
      geom_point() +
      geom_point(data = centroid, aes(PCoA1, PCoA2, color = Visit_number), shape = 15, size = 6) +
      stat_ellipse(type = 'norm') +
      annotate(geom = 'text', x = Inf, y = Inf, label = to_print, hjust = 1.1, vjust = 1.1) + 
      theme(legend.position = 'bottom', panel.background = element_blank(), 
            panel.grid = element_line(colour = 'lightgrey')) + 
      labs(title = paste('PCoA2','vs','PCoA1','colored by', 'Visit number',paste('(',choice,')',sep = '')))
    

    print(g)
  }
  return(res_adonis)
}

do_visit_adonis(visit_to_compare = c(1:4),choice = 'genus',cores = 4)
do_visit_adonis(visit_to_compare = c(1:4),choice = 'species',cores = 4)


#get the significance variables from independence tests
significant_result_independence <- read.csv('results/independence_test/significant_variable_independence_test.csv')[,-c(1,4)]
significant_result_independence$PCoA <- sapply(significant_result_independence$PCoA, function(x) substring(x, 5,5))

to_analyse_genus <- subset(significant_result_independence, Classification == 'genus')[,-3] 
to_analyse_species <- subset(significant_result_independence, Classification == 'species')[,-3]

#see the axes which present significances with at least one categorical column
axes_per_variable_genus <- list()
for (var in unique(to_analyse_genus$Variable)) {
  df <- subset(to_analyse_genus, Variable == var)
  unique_axes <- as.numeric(unique(df$PCoA))
  
  if(length(unique_axes) == 1){
    if(unique_axes[1] != 1){
      axes_per_variable_genus[[var]] <- as.numeric(c(1, unique_axes[1]))
    }else{
      axes_per_variable_genus[[var]] <- as.numeric(c(unique_axes[1],2))
    }
  }else{
    vect <- c()
    for (i in 1:(length(unique_axes) - 1)) {
      pair <- as.numeric(c(unique_axes[i], unique_axes[i + 1]))
    }
    
    last_pair <- as.numeric(c(unique_axes[length(unique_axes)], unique_axes[1]))
    if (!any(vect == rev(last_pair))) {
      vect <- c(vect, last_pair)
    }
    
    axes_per_variable_genus[[var]] <- vect
  }
}

axes_per_variable_species <- list()
for (var in unique(to_analyse_species$Variable)) {
  df <- subset(to_analyse_species, Variable == var)
  unique_axes <- as.numeric(unique(df$PCoA))
  
  if(length(unique_axes) == 1){
    if(unique_axes[1] != 1){
      axes_per_variable_species[[var]] <- as.numeric(c(1, unique_axes[1]))
    }else{
      axes_per_variable_species[[var]] <- as.numeric(c(unique_axes[1],2))
    }
  }else{
    vect <- c()
    for (i in 1:(length(unique_axes) - 1)) {
      pair <- as.numeric(c(unique_axes[i], unique_axes[i + 1]))
    }
    
    last_pair <- as.numeric(c(unique_axes[length(unique_axes)], unique_axes[1]))
    if (!any(vect == rev(last_pair))) {
      vect <- c(vect, last_pair)
    }
    
    axes_per_variable_species[[var]] <- vect
  }
}

#make the analyses and generate the images

for(var in unique(to_analyse_genus$Variable)){
  do_analysis(col = var,choice = 'genus',visit = 1:4,plot = T,rep = T,pcoa_to_plot = axes_per_variable_genus[[var]])
}



for(var in unique(to_analyse_species$Variable)){
  do_analysis(col = var,choice = 'species',visit = 1:4,plot = T,rep = T,pcoa_to_plot = axes_per_variable_species[[var]])
}







