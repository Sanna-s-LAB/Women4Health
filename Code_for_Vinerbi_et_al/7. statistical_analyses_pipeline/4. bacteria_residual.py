# Python script
# This script compute the residuals of bacterial abundances by significant results by indipendence tests
# Author: Fabio Chillotti (fabiochillotti@cnr.it)
# Last change: 13/01/2025
# Python version: v3.10.12

from utilities import GeneralizedModel #see the script utilities.py

#import the necessary packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import sys
remove_by_prevalence = sys.argv[1] #from the command line it must be provide whether filtering the bacteria for prevalence and mean abundances (True or False)

if remove_by_prevalence == 'True': 
  remove_by_prevalence = True
else:
  remove_by_prevalence = False

#get the abundances
clr_genus = pd.read_csv('data/clr_genus.csv') 
clr_genus.columns = ['Code'] + list(clr_genus.columns[1::])

clr_species = pd.read_csv('data/clr_species.csv')
clr_species.columns = ['Code'] + list(clr_species.columns[1::])


#get the lab data
df_lab = pd.read_csv('data/Laboratory_data.csv').drop('Unnamed: 0',axis = 1)

#get the phenotype and questionnaire
df_feno = pd.read_csv('data/clean_phenotypes.csv')
df_feno = df_feno.drop('Unnamed: 0', axis = 1)
quest = pd.read_csv('data/clean_questionnaire.csv')

#get the relative abundance for genus and species
rel_abb_genus = pd.read_csv('data/Relative_abundances_genus.csv',
                           sep = '\t') 

rel_abb_genus = rel_abb_genus.transpose()

rel_abb_species = pd.read_csv('data/Relative_abundances_species.csv',
                             sep = '\t')

rel_abb_species = rel_abb_species.transpose()

#get the taxa that occurs more than 20% of the samples
right_genera = (((rel_abb_genus != 0).sum(axis = 0) / len(rel_abb_genus)) > 0.2)[(((rel_abb_genus != 0).sum(axis = 0) / len(rel_abb_genus)) > 0.2)].index
right_species = (((rel_abb_species != 0).sum(axis = 0) / len(rel_abb_species)) > 0.2)[(((rel_abb_species != 0).sum(axis = 0) / len(rel_abb_species)) > 0.2)].index

right_genera = list(right_genera)
rel_abb_genus = rel_abb_genus[right_genera]
right_species = list(right_species)
rel_abb_species = rel_abb_species[right_species]

#get the taxa with mean abundances greater than 0.1
right_genera = list(rel_abb_genus.apply(np.mean)[rel_abb_genus.apply(np.mean) > 0.1].index)
right_genera.append('Code')

right_species = list(rel_abb_species.apply(np.mean)[rel_abb_species.apply(np.mean) > 0.1].index)
right_species.append('Code')

if remove_by_prevalence: #filtering for prevalence and mean abundance if True is provided
  print('Filtering by prevalence and mean abundance...')
  clr_genus = clr_genus[right_genera]
  clr_species = clr_species[right_species]
else:
  pass

#create the complete datasets with bacteria and metadata
df_complete_genus = pd.merge(clr_genus,quest,on = 'Code')
df_complete_genus = pd.merge(df_complete_genus,df_feno, on = 'Code')
df_complete_genus = pd.merge(df_complete_genus,df_lab, on = 'Code')
df_complete_genus['Visit'] = df_complete_genus['Code'].apply(lambda x: int(x.split('_')[1]))

df_complete_species = pd.merge(clr_species,quest,on = 'Code')
df_complete_species = pd.merge(df_complete_species,df_feno, on = 'Code')
df_complete_species = pd.merge(df_complete_species,df_lab)
df_complete_species['Visit'] = df_complete_species['Code'].apply(lambda x: int(x.split('_')[1]))

bacteria_genus = list(clr_genus.drop('Code',axis = 1).columns)
#bacteria_genus = [bact for bact in bacteria_genus if not bact.startswith('Unclassified')]

bacteria_species = list(clr_species.drop('Code',axis = 1).columns)
#bacteria_species = [bact for bact in bacteria_species if not bact.startswith('Unclassified')]

hormones = ['PRL','PROG','BES17','FSH','LH']

from sklearn.linear_model import LinearRegression

def grab_residuals(model, df_complete):       #function to create datasets from the attribute residuals of fitted models  
    
    res_to_save = {}

    bacteria = np.unique(model.residuals.Outcome.values)

    for bact in bacteria:
        res_to_save[bact] = model.residuals[model.residuals.Outcome == bact].Residual
        res_to_save[bact].index = range(212)
    to_save = pd.DataFrame.from_dict(res_to_save)
    to_save['Code'] = df_complete['Code']

    return(to_save)

f = open("data/independence_significant_variables.txt", "r")
variables = [el.split('\n')[0] for el in f.readlines()]

#create a short name (depending on the results of independence tests)
short_names = {
    'tech' : 'Tech',
    'tech_Age_category' : 'Tech_Age',
    'tech_Age_category_Pregnancy_category' : 'Tech_Age_Preg',
    'tech_Age_category_Pregnancy_category_Pill_use': 'Tech_Age_Preg_Pill',
    'tech_Age_category_Pregnancy_category_Pill_use_Swab_morning_or_not' : 'Tech_Age_Preg_Pill_SwabM',
    'tech_Age_category_Pregnancy_category_Pill_use_Swab_morning_or_not_Swab_after_feces' : 'Tech_Age_Preg_Pill_SwabM_A',
    'tech_Age_category_Pregnancy_category_Pill_use_Swab_morning_or_not_Swab_after_feces_WHR_ranges' : 'Tech_Age_Preg_Pill_SwabM_A_WHR',
    'tech_Age_category_Pregnancy_category_Pill_use_Swab_morning_or_not_Swab_after_feces_WHR_ranges_Bristol_stool_scale' : 'Tech_Age_Preg_Pill_SwabM_A_WHR_Bristol',
    'tech_Age_category_Pregnancy_category_Pill_use_Swab_morning_or_not_Swab_after_feces_WHR_ranges_Bristol_stool_scale_sex_less_than_two_days' : 'Tech_Age_Preg_Pill_SwabM_A_WHR_Bristol_coitus'
} #IT'S IMPORTANT TO ORDER THE FILE "independence_significant_variables.txt" as the short names above (they must be coherent)
#Be sure that in the txt file there are not extra spaces or \n as it would raise an error in this script

#to store the residuals
dict_of_data_genus = {}
dict_of_data_species = {}

print('Computing residuals...')
for cov in [['Qubit_DNA','Total_counts','Qubit_Library']]+ [['Qubit_DNA','Total_counts','Qubit_Library'] + variables[:i] for i in range(1, len(variables) + 1)]:
    
    print('Bacterium' + ' ~ ' + ' + '.join(cov))
    cov_key = "_".join(['tech'] + cov[3:])

    # Process genus
    model_g = GeneralizedModel(df_complete_genus, outcomes=bacteria_genus, covariates=cov,
                               model=LinearRegression())
    model_g.compute_resid()
    dict_of_data_genus[cov_key] = grab_residuals(model_g, df_complete_genus)

    # Process species
    model_s = GeneralizedModel(df_complete_species, outcomes=bacteria_species, covariates=cov,
                               model=LinearRegression())
    model_s.compute_resid()
    dict_of_data_species[cov_key] = grab_residuals(model_s, df_complete_species)
    
dict_of_data_genus = {short_names[key]: value for key, value in dict_of_data_genus.items() if key in short_names}
dict_of_data_species = {short_names[key]: value for key, value in dict_of_data_species.items() if key in short_names}

#save all residuals
i = 1
for corr in dict_of_data_genus.keys():
    name_file = str(i) + '.' + corr + '.csv'
    if remove_by_prevalence:
      dict_of_data_genus[corr].to_csv('data/corrected_clr/filtered/genus/' + name_file)
    else:
      dict_of_data_genus[corr].to_csv('data/corrected_clr/genus/' + name_file)
    i = i + 1

i = 1
for corr in dict_of_data_species.keys():
    name_file = str(i) + '.' + corr + '.csv'
    
    if remove_by_prevalence:
      dict_of_data_species[corr].to_csv('data/corrected_clr/filtered/species/' + name_file)
    else:
      dict_of_data_species[corr].to_csv('data/corrected_clr/species/' + name_file)
    i = i + 1
