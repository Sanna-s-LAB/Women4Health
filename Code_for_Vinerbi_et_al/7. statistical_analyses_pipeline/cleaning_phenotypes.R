df_feno <- read.csv('/mnt/sannaLAB-Data/PRJ_WOMEN4HEALTH/BLOOD_BIOCHEMISTRY_PHENOTYPES/W4H_Trieste_Dati_volontarie_allsamples_noduplicates13122024.csv')

#fix the id
fix_id <- function(stringa){
  stringa_number <- substring(stringa,4,6)
  stringa_finale <- paste('X',stringa_number,sep = '')
  return(stringa_finale)
}

df_feno$Code <- sapply(df_feno$Record.ID, fix_id) 
df_feno$Code <- paste(df_feno$Code, df_feno$Numero.Visita, sep = '_')

df_feno <- df_feno[,-which(colnames(df_feno) %in% c("Record.ID","Numero.Visita","Age","sesso",
                                                    "Note1","Note2","Note3","Note.4",
                                                    'stringa', 'Codice.accettazione',
                                                    'data.accettazione.',"data.analisi.profilo.lipidico",
                                                    "data.analisi.profilo.ormonale",
                                                    "data.spedizione","data.spedizione"))]
colnames(df_feno)
colnames(df_feno) <- c('GL','COL','HDL','LDL','TRI','AST','ALT','INS','HOMA_IR',
                       'HOMA_B','PROG','TSH','FT4','LH','FSH','BES17','TST','PRL','Bacth','Code')
#write.csv(df_feno,'~/complete_pipeline/data/sintetic_phenotypes.csv')
