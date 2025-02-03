df_feno <- read.csv('data/clean_phenotypes.csv')[,-1]
clr_genus <- read.csv('data/corrected_clr/filtered/genus/1.Tech.csv')[,-1]
clr_species <- read.csv('data/corrected_clr/filtered/species/1.Tech.csv')[,-1]

genus_complete <- merge(clr_genus,df_feno, by = 'Code')
genus_complete$Visit <- sapply(genus_complete$Code, function(x){substring(x,6,6)})

species_complete <- merge(clr_species,df_feno, by = 'Code')
species_complete$Visit <- sapply(species_complete$Code, function(x){substring(x,6,6)})

res_complete <- read.csv('results/linear_models/results_hormones.csv')[,-1]
sig_results <- subset(res_complete,p.value < 0.05)

make_scatterplot <- function(horm, bact, normalize = FALSE) {
  if (bact %in% colnames(genus_complete)) {
    df <- genus_complete
  } else if (bact %in% colnames(species_complete)) {
    df <- species_complete
  }
  
  # Create a mapping of Visit values to new labels
  df$Visit_label <- factor(df$Visit, 
                           levels = c(1, 2, 3, 4), 
                           labels = c("F", "O", "EL", "LL"))
  
  if (normalize) {
    b <- bestNormalize::orderNorm(df[[horm]])
    df[[horm]] <- b$x.t
  }
  
  ggplot(df) + 
    geom_point(aes(!!sym(bact), !!sym(horm), color = Visit)) + 
    geom_smooth(aes(!!sym(bact), !!sym(horm), color = Visit), method = 'lm', se = FALSE) +
    theme_bw() +
    facet_wrap(~Visit_label, ncol = 2)
}

list_of_plot <- list()
for(taxa in sig_results$Bacteria){
  hormones_to_consider <- subset(sig_results,Bacteria == taxa)$Hormone
  for(horm in hormones_to_consider){
    list_of_plot[[paste0(horm,'_vs_',taxa)]] <- make_scatterplot(horm, taxa)
  }
}

# path to save
output_path <- "results/linear_models/scatterplot/"

# save plots
i <- 1
for (p in list_of_plot) {
  ggsave(
    filename = paste0(output_path, "scatterplot_", names(list_of_plot)[i], ".png"),
    plot = p,
    width = 8,
    height = 6)
  i <- i + 1
}



