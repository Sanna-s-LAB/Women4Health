# Rscript
# This script runs independence tests between women with changing CSTs and those with stable CSTs.
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 31/01/2025
# R version: R v4.4.1

library(ggplot2)
library(lme4)
library(reshape2)

############## CHECK IF THERE ARE DIFFERENCES BETWEEN CHANGING WOMAN AND NOT ############## 

#features to check
col_to_check <- c('Age_category','Smoking_habits','Swab_morning_or_not','Swab_after_feces',
                  'sex_less_than_two_days','Pregnancy_category','BMI_ranges',
                  'WHR_ranges','Pill_use',
                  'Med_or_supp_check','Bristol_stool_scale')
#get the data
quest <- read.csv('~/complete_pipeline/data/clean_questionnaire.csv')

cst_class <- read.csv('~/complete_pipeline/data/Valencia_results.csv', sep = '\t')[,c('sampleID','subCST','CST')]
colnames(cst_class)[1] <- 'Code'
cst_class$Woman <- sapply(cst_class$Code, function(string){substring(string,1,4)})
cst_class$Visit_number <- sapply(cst_class$Code, function(string){substring(string,6,6)})

its_counts <- read.csv('~/complete_pipeline/data/Counts_ITS.csv', sep = '\t')[,1:3]

changing_women <- c("X004", "X008", "X012", "X016", "X030", "X058", "X061", "X069", "X083", "X091", "X092")
quest$CST_change <- ifelse(quest$Patient_id %in% changing_women,'Change','No_change')

#we shall not consider the women at the enrollment
quest <- subset(quest, Visit_number != 0)

p_values <- list()  #initialize empty list
for(col in col_to_check){
  df <- quest[,c(col,'Patient_id','CST_change')] #get the data
  if(col %in% c('Age_category','Smoker','BMI_ranges','WHR_ranges','Pill_use',
                'Ex_smoker','Pregnancy_category')){ #Test on steady variables (Fisher)
    df <- aggregate(. ~ Patient_id, data = df, FUN = function(vect){vect[1]})
    tab <- table(df[[col]],df[['CST_change']])
    #do the Fisher test with or without simulation depending on the variables's number of levels
    if(length(unique(quest[[col]])) > 2){ 
      test <- fisher.test(tab, simulate.p.value = TRUE)
    } else {
      test <- fisher.test(tab)
    }
    p_values[[col]] <- test$p.value
  }else{ #Test on weekly changing variables (Generalized logistic model)
    if(length(unique(df[[col]])) > 1) { 
      df$CST_change <- as.factor(ifelse(df$CST_change == "Change", 1, 0))
      formula <- as.formula(paste0('CST_change ~ ', col, ' + (1|Patient_id)'))
      model <- glmer(formula, data = df, family = binomial)
      p_value <- summary(model)$coefficients[2, "Pr(>|z|)"]
      p_values[[col]] <- p_value
    }
  }
}

#store results in a dataframe
p_values <- as.data.frame(t(as.data.frame(p_values)))
p_values$Feature <- rownames(p_values)
rownames(p_values) <- NULL
colnames(p_values)[1] <- 'p_value'
p_values <- p_values[,c(2,1)]
p_values$'-log10 p_value' <- -log10(p_values$p_value)
p_values$order <- 1:nrow(p_values)

#plot the results
plot <- ggplot(p_values) + geom_bar(aes(x = reorder(Feature,order),
                                        y = !!sym('-log10 p_value')),stat = 'identity',
                                    fill = 'lightgreen',color = 'black') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.45),
        panel.background = element_rect(fill = 'white'), panel.grid = element_line(colour = 'grey'),
        axis.text = element_text(size = 13), axis.title = element_text(size = 13)) + 
  geom_abline(slope = 0,intercept = -log10(0.05),color = 'red') +
  annotate(geom = 'text',x = 4,y = 1.25,label = '**significance_line**',size = 4) + xlab('Feature')

plot

