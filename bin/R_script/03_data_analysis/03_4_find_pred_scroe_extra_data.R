## find pred_score extra data 
library(tidyverse)
## TME included
load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")

pred_score_extra_data_update_TME <- read.table("results/extra_data/04_1_Not_validated_Validated_neoepitopes_extra_predicted_score_tme_included.txt", sep = "\t", header = T)

pred_score_extra_data_update_TME  <- pred_score_extra_data_update_TME %>% 
  rename( "prediction_rf_tme" = 'mean_prediction_rf') %>% 
  mutate(identity = paste(Sample ,HLA_allele, Mut_peptide, sep = "_"))

# 

## TME excluded
pred_score_extra_data_update_TME_excluded <- read.table("results/extra_data/04_1_Not_validated_Validated_neoepitopes_extra_predicted_score_tme_excluded.txt", sep = "\t",header = T)

Pred_Modelling_extra <- pred_score_extra_data_update_TME_excluded %>% 
  rename("prediction_rf" = 'mean_prediction_rf') %>% 
  mutate(identity = paste(Sample ,HLA_allele, Mut_peptide, sep = "_")) %>% 
  select(identity, prediction_rf) %>% 
  right_join(.,pred_score_extra_data_update_TME )
  

# join df 
# --------------------------------

Pred_Modelling_extra$Patient <- gsub("Neye","MM",Pred_Modelling_extra$Patient)
Pred_Modelling_extra$Patient <- gsub("BC","mUC",Pred_Modelling_extra$Patient)


cols_selected = c('Patient','HLA_allele','Mut_peptide','Aro', 'Inst', 'CysRed','RankEL','RankBA','NetMHCExp',
                    'Expression','SelfSim','Prime','PropHydroAro','HydroCore','pI',
                    'PropSmall','PropAro','PropBasic','PropAcidic','DAI','Stability','Foreigness',
                    'CelPrev','PrioScore','CYT','HLAexp','Monocytes',
                     "prediction_rf","prediction_rf_tme")

Pred_Modelling_extra <- Pred_Modelling_extra %>% select(all_of(cols_selected))
Pred_Modelling_to_join <- Pred_Modelling %>% select(all_of(cols_selected))

Pred_Modelling_total_all_peptide <- Pred_Modelling_to_join %>% 
  bind_rows(.,Pred_Modelling_extra) 


save(Pred_Modelling_total_all_peptide,file = "data/04_plotting/Survival_data/All_pep_for_survival.Rdata")

