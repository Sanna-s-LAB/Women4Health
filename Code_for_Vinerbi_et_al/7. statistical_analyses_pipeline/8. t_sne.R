# Rscript
# This script run the tsne dimensionality reduction with visualization
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 31/01/2025
# R version: R v4.4.1




library(ggplot2)
library(gridExtra)
library(grid)
library(Rtsne)
library(vegan)
library(dplyr)
library(ggpubr)

#get the data
clr_genus <- read.csv('data/clr_genus.csv')
rownames(clr_genus) <- clr_genus$X
clr_genus <- clr_genus[,-1]
clr_species <- read.csv('data/clr_species.csv')
rownames(clr_species) <- clr_species$X
clr_species <- clr_species[,-1]

clr_genus <- clr_genus[rownames(clr_species),]

cst_class <- read.table('data/Tab_CST_SubCST.csv', header = T)
colnames(cst_class)[1] <- 'Code'
rownames(cst_class) <- cst_class$Code
cst_class <- cst_class[rownames(clr_genus),]

rel_abb_species <- t(read.table('data/Relative_abundances_species.csv'))
rel_abb_species <- as.data.frame(rel_abb_species)
rel_abb_species <- rel_abb_species[rownames(clr_genus),]

######TSNE######

#fixing the seed due to randomness of tsne
set.seed(1)

#compute the tsne for species
tsne_species <- Rtsne(clr_species)

#create che dataframe with metadata useful to plots
tsne_plot_species <- data.frame(x = tsne_species$Y[,1], 
                        y = tsne_species$Y[,2],
                        Code = cst_class$Code,
                        subCST = cst_class$subCST,
                        CST = cst_class$CST)

#####Plots like the ISALA paper#####

#FIRST PLOT
color_sub <- c('#3B75AF','#7F7F7F','#A8DD93','#B3C6E5','#D57DBE','#8D69B8','#EF8636','#F5BE82',
               '#C53A32','#84584D','#509E3E','#BCBD45','#F19D99')

names(color_sub) <- c("I-A","I-B","II","III-A","III-B","IV-A","IV-B","IV-C0","IV-C1","IV-C4","V",
                      "IV-C2","IV-C3")
  
tsne_plot_species_ISALA <- tsne_plot_species[,1:5]

tsne_plot_species_ISALA$subCST <- factor(tsne_plot_species_ISALA$subCST, 
                                              levels = names(color_sub))

tsne_plot_species_ISALA$CST <- factor(tsne_plot_species_ISALA$CST)

#Create a named color scale
color_scale <- scale_color_manual(values = color_sub)

#FIRST PLOT
g1 <- ggplot(tsne_plot_species_ISALA) + 
  geom_point(aes(x = x, y = y, color = subCST), size = 3) + 
   theme_bw() + theme(legend.position = 'bottom',axis.title.y = element_blank()) + color_scale
g1

g1_beta <- ggplot(tsne_plot_species_ISALA) + 
  geom_point(aes(x = x, y = y, color = CST),size = 4) + theme(legend.position = 'bottom') +
  scale_colour_manual(values = c('I'="#FF8C6A",'II'="#FFDA6A",'III'="#CAD94A",'IV'="#6AFF8A",'V'= "#6A9EFF")) +
  theme_bw()

# SECOND PLOT
most_abundant_genus <- apply(clr_genus,1, function(row){names(row)[which.max(row)]})

lactobacillus_species <- c('Lactobacillus.iners','Lactobacillus.crispatus','Lactobacillus.gasseri','Lactobacillus.jensenii')

most_abundant_lactobacillus <- apply(clr_species[,lactobacillus_species], 1, function(row) {names(row)[which.max(row)]})

tsne_plot_species_ISALA$most_abundant <- ifelse(most_abundant_genus == "Lactobacillus",most_abundant_lactobacillus,
                                                most_abundant_genus)

tsne_plot_species_ISALA$most_abundant <- ifelse(tsne_plot_species_ISALA$most_abundant == "Bifidobacterium",'Gardnerella',
                                                ifelse(tsne_plot_species_ISALA$most_abundant %in% c('Corynebacterium','Prevotella',
                                                                                                    'Lactobacillus.crispatus','Lactobacillus.iners',
                                                                                                    'Lactobacillus.gasseri'),
                                                       tsne_plot_species_ISALA$most_abundant,'Other'))

color_sub <- c(
  'Lactobacillus.iners' = "#A6CEE3",
  'Lactobacillus.crispatus' = "#1F78B4",
  'Other' = "#B2DF8A",
  'Lactobacillus.gasseri' = "#33A02C",
  'Lactobacillus.jensenii' = "#E31A1C",
  'Prevotella' = "#6A3D9A",
  'Corynebacterium' = "#FF7F00",  # Arancione brillante
  'Gardnerella' = "#CAB2D6"       # Lilla chiaro
)



tsne_plot_species_ISALA$most_abundant <- factor(tsne_plot_species_ISALA$most_abundant)

# Create a named color scale
color_scale <- scale_color_manual(values = color_sub)

g2 <- ggplot(tsne_plot_species_ISALA) + 
  geom_point(aes(x = x, y = y, color = most_abundant),size = 3) +
  color_scale + theme_bw() + theme(legend.position = 'bottom',axis.title.x = element_blank(),
                                   legend.text = element_text(face = 'italic'))
g2

# THIRD PLOT
# Find the second most abundant genus
second_most_abundant_genus <- apply(clr_genus, 1, function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  names(row)[sorted_indices[2]]
})

# Find the second most abundant Lactobacillus species
second_most_abundant_lactobacillus <- apply(clr_species[, lactobacillus_species], 1, 
                                            function(row) {
                                                          sorted_indices <- order(row, decreasing = TRUE)
                                                          names(row)[sorted_indices[2]]
                                                          })

# Update the tsne_plot_species_ISALA dataframe to reflect the second most abundant bacteria
tsne_plot_species_ISALA$second_most_abundant <- ifelse(second_most_abundant_genus == "Lactobacillus", 
                                                       second_most_abundant_lactobacillus,
                                                       second_most_abundant_genus)

# Modify the grouping of bacteria as before
tsne_plot_species_ISALA$second_most_abundant <- ifelse(tsne_plot_species_ISALA$second_most_abundant == "Bifidobacterium", 
                                                       'Gardnerella',
                                                       ifelse(tsne_plot_species_ISALA$second_most_abundant %in% c('Corynebacterium', 'Prevotella',
                                                                                                                  'Lactobacillus.crispatus', 'Lactobacillus.iners',
                                                                                                                  'Lactobacillus.gasseri'),
                                                              tsne_plot_species_ISALA$second_most_abundant, 'Other'))

# Assign factor levels to the second most abundant bacteria
tsne_plot_species_ISALA$second_most_abundant <- factor(tsne_plot_species_ISALA$second_most_abundant)

# Create the named color scale
color_scale <- scale_color_manual(values = color_sub)

# Generate the ggplot with the second most abundant bacteria
g3 <-  ggplot(tsne_plot_species_ISALA) + 
   geom_point(aes(x = x, y = y, color = second_most_abundant), size = 3) + 
   color_scale + theme_bw() + theme(legend.position = 'bottom', 
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_blank())
g3

#FOURTH PLOT
tsne_plot_species_ISALA$most_abundant_percentage <- apply(rel_abb_species,1, function(row){max(row)})
View(tsne_plot_species_ISALA)

color_gradient <- scale_color_gradientn(colors = c("#160883", "#74159F", "#803257",
                                                   "#E1815B","#F3ED57"))

# Generate the ggplot with the custom color scale
g4 <- ggplot(tsne_plot_species_ISALA) + 
  geom_point(aes(x = x, y = y, color = most_abundant_percentage), size = 3) + 
  color_gradient + 
  theme_bw() + theme(legend.position = 'bottom')
g4

#Plot with arrows indicating CST shifts
tsne_plot_species_ISALA$Patient <- sapply(tsne_plot_species_ISALA$Code,function(x) substring(x,1,4))
tsne_plot_species_ISALA$Visit_number <- sapply(tsne_plot_species_ISALA$Code,function(x) substring(x,6,6))

#Doing the math
changing_women <- NULL #find the changing women
i <- 1
for(ID in tsne_plot_species_ISALA$Patient){
  df <- subset(tsne_plot_species_ISALA,Patient == ID)
  if(length(unique(df$subCST)) != 1){
    changing_women[[i]] <- ID
    i <- i + 1
  }
}
 
list_changing <- list() #find where the change is happened
i <- 1
for(ID in changing_women){
  df <- subset(tsne_plot_species_ISALA, Patient == ID)
  changing <- NULL  # Initialize as an empty vector to accumulate changes for each ID
  
  for(visit in df$Visit_number){
    visit <- as.numeric(visit)
    if(visit %in% 1:3){
      next_visit <- visit + 1
      
      # Ensure that the next visit exists in the data
      if(next_visit %in% df$Visit_number){
        starting_value <- as.vector(subset(df, Visit_number == visit)$subCST)
        final_value <- as.vector(subset(df, Visit_number == next_visit)$subCST)
        
        # Compare the subCST values
        if(starting_value != final_value){
          change <- paste(visit, '-', next_visit, ': ', starting_value, '->', final_value, sep = '')
          changing <- c(changing, change)  # Append to the vector
        }
      }
    }
  }
  
  # Store the changes in the list
  list_changing[[ID]] <- changing
}

g_special <- g1 #let's create the plot with arrows
text_data <- data.frame()
for(ID in names(list_changing)){
  changes <- list_changing[[ID]]
  for(change in changes){
    start_visit <- substring(change,1,1)
    end_visit <- substring(change,3,3)
    start_point <- subset(tsne_plot_species_ISALA,(Patient == ID) & (Visit_number == start_visit))
    end_point <- subset(tsne_plot_species_ISALA,(Patient == ID) & (Visit_number == end_visit))
    g_special <- g_special + geom_segment(x = start_point$x, y = start_point$y, 
                                              xend = end_point$x, yend = end_point$y,
                                          arrow = arrow(length = unit(1.5, "mm")), color = "black",linewidth = 0.1)
    
    text_data <- rbind(text_data, data.frame(
      x = (start_point$x + end_point$x) / 2,
      y = (start_point$y + end_point$y) / 2,
      label = paste(start_visit,'->',end_visit)
    ))
      
      # Add text annotations to the plot
      g_special <- g_special + geom_text(
        data = text_data,
        aes(x = x, y = y, label = label),
        color = "black",
        size = 3,
        vjust = -1,  # Adjust vertical position
        hjust = 0.5  # Adjust horizontal position
      )
    
    
  }
}
g_special

# Create a data frame for visit changes
change_counts <- data.frame()
for(ID in names(list_changing)){
  changes <- list_changing[[ID]]
  for(change in changes){
    start_visit <- substring(change,1,1)
    end_visit <- substring(change,3,3)
    change_counts <- rbind(change_counts, data.frame(
      change = paste(start_visit, '->', end_visit)
    ))
  }
}

# Count the occurrences of each type of change
change_summary <- as.data.frame(table(change_counts$change))
names(change_summary) <- c("change", "count")

# Calculate proportions
change_summary$prop <- change_summary$count / sum(change_summary$count)
data <- change_summary %>% 
  arrange(desc(change)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
pastel_colors <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6")

#let's create a pie chart to see the proportion of visit between which the changes are happened
pie_chart <- ggplot(data, aes(x = "", y = prop, fill = change)) +
  geom_bar(width = 1, stat = "identity",color = 'white') +
  coord_polar(theta = "y",start = 0) +
  theme_void() +
  theme(legend.position = "none") + 
  geom_text(aes(y = ypos, label = change), color = "black", size=2) +
  scale_fill_manual(values = pastel_colors)  # Use the custom pastel colors
pie_chart_grob <- ggplotGrob(pie_chart)

#wrap the two plots
g_final <- g_special + annotation_custom(grob = pie_chart_grob,xmin = -10,xmax = -5,ymin = -15,ymax = -10)

g_final


#INNOVATIVE PLOT (changing between CST)
list_changing <- list()
i <- 1
for(ID in changing_women){
  df <- subset(tsne_plot_species_ISALA, Patient == ID)
  changing <- NULL  # Initialize as an empty vector to accumulate changes for each ID
  
  for(visit in df$Visit_number){
    visit <- as.numeric(visit)
    if(visit %in% 1:3){
      next_visit <- visit + 1
      
      # Ensure that the next visit exists in the data
      if(next_visit %in% df$Visit_number){
        starting_value <- as.vector(subset(df, Visit_number == visit)$CST)
        final_value <- as.vector(subset(df, Visit_number == next_visit)$CST)
        
        # Compare the subCST values
        if(starting_value != final_value){
          change <- paste(visit, '-', next_visit, ': ', starting_value, '->', final_value, sep = '')
          changing <- c(changing, change)  # Append to the vector
        }
      }
    }
  }
  
  # Store the changes in the list
  list_changing[[ID]] <- changing
}

#map visit_number to the phases
visit_labels <- c("1" = "F", "2" = "O", "3" = "EL", "4" = "LL")

g_special_beta <- g1_beta + theme(legend.text = element_text(size =7), axis.title.x = element_blank(),
                                  axis.title.y = element_blank())
text_data <- data.frame()

for (ID in names(list_changing)) {
  changes <- list_changing[[ID]]
  
  for (change in changes) {
    start_visit <- substring(change, 1, 1)
    end_visit <- substring(change, 3, 3)
    
    #substitute the number with the phases' name
    start_label <- visit_labels[start_visit]
    end_label <- visit_labels[end_visit]
    
    start_point <- subset(tsne_plot_species_ISALA, (Patient == ID) & (Visit_number == start_visit))
    end_point <- subset(tsne_plot_species_ISALA, (Patient == ID) & (Visit_number == end_visit))
    
    #add arrow
    g_special_beta <- g_special_beta + geom_segment(
      x = start_point$x, y = start_point$y, 
      xend = end_point$x, yend = end_point$y,
      arrow = arrow(length = unit(1.5, "mm")), color = "black", linewidth = 0.5
    )
    
    #create the text to print above the arrows
    text_data <- rbind(text_data, data.frame(
      x = (start_point$x + end_point$x) / 2,
      y = (start_point$y + end_point$y) / 2,
      label = paste(start_label, '->', end_label)
    ))
  }
}

#add annotations
g_special_beta <- g_special_beta + geom_text(
  data = text_data,
  aes(x = x, y = y, label = label),
  color = "black",
  size = 5,
  vjust = -1.4,
  hjust = 0.5
)

#show the plot
g_special_beta

# Create a mapping between visit numbers and the desired labels
visit_labels <- c("1" = "F", "2" = "O", "3" = "EL", "4" = "LL")

# Create a data frame for visit changes
change_counts <- data.frame()
for(ID in names(list_changing)){
  changes <- list_changing[[ID]]
  
  for(change in changes){
    # Extract start and end visit numbers
    start_visit <- substring(change, 1, 1)
    end_visit <- substring(change, 3, 3)
    
    # Replace visit numbers with their corresponding labels
    start_label <- visit_labels[start_visit]
    end_label <- visit_labels[end_visit]
    
    # Create a new change label using the mapped visit labels
    change_counts <- rbind(change_counts, data.frame(
      change = paste(start_label, '->', end_label)
    ))
  }
}

# Count the occurrences of each type of change
change_summary <- as.data.frame(table(change_counts$change))
names(change_summary) <- c("change", "count")

# Calculate proportions
change_summary$prop <- change_summary$count / sum(change_summary$count)

# Arrange the data and calculate label positions for the pie chart
data <- change_summary %>% 
  arrange(desc(change)) %>%
  mutate(ypos = cumsum(prop) - 0.5 * prop)

# Define custom pastel colors
pastel_colors <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6")

# Create the pie chart
pie_chart <- ggplot(data, aes(x = "", y = prop, fill = change)) +
  geom_bar(width = 1, stat = "identity", color = 'white') +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(legend.position = "none",text = element_text(family = "sans")) +
  geom_text(aes(y = ypos, label = change), color = "black", size = 4) +
  scale_fill_manual(values = pastel_colors)  # Use custom pastel colors

# Convert the pie chart to a grob
pie_chart_grob <- ggplotGrob(pie_chart)

# Add the pie chart to the existing plot (g_special_beta)
g_final_beta <- g_special_beta + annotation_custom(
  grob = pie_chart_grob,
  xmin = -30,  # Increase horizontal space
  xmax = 9,   
  ymin = -15,  # Increase vertical space
  ymax = -5
) + theme(text = element_text(family = "sans"))

# Display the final plot
g_final_beta

ggsave('results/tsne/tsne_subcst.jpg',g1,width = 7,height = 7)
ggsave('results/tsne/tsne_most_abundant.jpg',g2,width = 7,height = 7)
ggsave('results/tsne/tsne_second_most_abundant.jpg',g3,width = 7,height = 7)
ggsave('results/tsne/tsne_most_abudant_percentage.jpg',g4,width = 7,height = 7)


g2 <- g2 + theme(legend.title = element_blank(),text = element_text(family = "sans"))
g3 <- g3 + theme(legend.title = element_blank(),text = element_text(family = "sans"))
g4 <- g4 + theme(legend.title = element_blank(),text = element_text(family = "sans"))
g1 <- g1 + theme(legend.title = element_blank(),text = element_text(family = "sans"))

combined_plot_A_B <- ggarrange(g2, g3, ncol = 2, nrow = 1, 
                               labels = c("A", "B"), font.label = list(face = "plain"),
                               common.legend = TRUE, legend = "bottom",
                               label.x = c(0, 0), label.y = c(1, 1))

combined_plot_C_D <- ggarrange(g4, g1, ncol = 2, nrow = 1, 
                               labels = c("C", "D"), font.label = list(face = "plain"),
                                legend = "bottom",
                               label.x = c(0, 0), label.y = c(1.02,1.02))


combined_plot <- ggarrange(combined_plot_A_B, combined_plot_C_D, ncol = 1,nrow = 2)

#save
png("results/tsne/t_SNE_complete.png",width = 3000,height = 1500,res = 200)
plot(combined_plot)
dev.off()


