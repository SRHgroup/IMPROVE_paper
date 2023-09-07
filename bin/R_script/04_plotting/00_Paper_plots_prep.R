# Paper plots 20221112
# -----------------------------

# -------------------------------------------------------------
####                     Plots for paper #######
# -------------------------------------------------------------
# plot in generel new 
# load libraries 
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggbeeswarm)
library(ggrepel)
library(cowplot)
library(ggpubr)
# for roc curves 

# core plot 
#install.packages("devtools")
#devtools::install_github("taiyun/corrplot", build_vignettes = TRUE)
library(GGally)
library(corrplot)
library(openxlsx)

# fore roc curves 
# Loading package
library(caTools)
library(randomForest)
library(pROC)
library(ROCit)
library(caret)
library(ROCR)
library(verification)

# survival packages 
library(survminer)
library(survival)

library(rstatix)

library(data.table)


# load data 
load('data/04_plotting/Rdata/04_pepides_for_plotting.Rdata')
source("bin/R_script/99_functions.R")

colnames(all_peptides)


# load data from RF model
# ----------------------------
#insert data 
# -------------------------------------------------
# NNalign
NNalign_results <- read.table("results/NNalign/all.test_pred_cmb")
NNalign_results <- NNalign_results %>% select(-V10)

colnames(NNalign_results) <- NNalign_results[1,]
NNalign_results <- NNalign_results[-1,]
NNalign_results <- as.tibble(NNalign_results)
NNalign_results <-  NNalign_results %>% 
  filter(!Measure=="Measure") %>% 
  mutate(response = case_when(Measure=="0.001000" ~ 0,
                              Measure=="1.000000" ~ 1 )) %>% 
  rename(.,   "HLA_allele" = HLA ) %>% 
  rename(.,  "Patient" =donor ) 

NNalign_results <- NNalign_results %>% 
  mutate(ID = paste(Patient,HLA_allele,Peptide, sep = "_")) %>% 
  select(ID,Prediction)


# load prediction results 
# ------------------------------------
Feature_round = 5
# IMPROVE 5 cv 
Feature_importance <- read.table(file = "results/5_fold_CV/TME_excluded/Feature_importance_TME_excluded.txt", header = T)
pred_df <- read.table(file = "results/5_fold_CV/TME_excluded/pred_df_TME_excludedprime_1.txt", header = T)
pred_df <- pred_df %>% mutate(ID = paste(Patient,HLA_allele,Mut_peptide, sep = "_")) %>% select(ID,prediction_rf)
# IMPROVE TME 5 cv 
Feature_importance_TME <- read.table("results/5_fold_CV/TME_included/Feature_importance_TME_included.txt", sep = "\t")
pred_df_TME_include <- read.table(file = "results/5_fold_CV/TME_included/pred_df_TME_includedprime_1.txt", header = T)
pred_df_TME_include <- pred_df_TME_include %>% 
  rename(.,  "prediction_rf_tme" = "prediction_rf" ) %>% 
  mutate(ID = paste(Patient,HLA_allele,Mut_peptide, sep = "_")) 

colnames(pred_df)
table(pred_df$response)
# -------------------------------


Pred_Modelling <- pred_df_TME_include %>% 
  right_join(.,pred_df, by = "ID") %>% 
  left_join(., NNalign_results, by = "ID")
nrow(Pred_Modelling)
# data handling 
# ---------------------------------
Pred_Modelling <- Pred_Modelling %>% 
  mutate(response_lab = case_when(response=="1" ~ "yes",
                                  TRUE ~ "no")) %>% 
  mutate(cohort = case_when(cohort == "melanoma" ~ "Melanoma", 
                            cohort == "bladder" ~ "mUC",
                            cohort == "Basket" ~ "Basket")) %>% 
  mutate(pep_length = str_length(Mut_peptide)) %>% 
  mutate(HLA_type = str_sub(HLA_allele, start = 1L, end = 5L)) %>% 
  mutate(HLA_num = str_sub(HLA_allele, start = 6L, end = 10L))

Pred_Modelling$Prediction <- as.numeric(Pred_Modelling$Prediction )
Pred_Modelling$Patient <- gsub("Neye","MM",Pred_Modelling$Patient)
Pred_Modelling$Patient <- gsub("BC","mUC",Pred_Modelling$Patient)

all_peptides$Patient <- gsub("BC","mUC",all_peptides$Patient)


length(unique(all_peptides$Patient))

# data handling 
# ---------------------------------
all_peptides <- all_peptides %>% 
  mutate(response_lab = case_when(response=="1" ~ "yes",
                                  TRUE ~ "no")) %>% 
  mutate(cohort = case_when(cohort == "melanoma" ~ "Melanoma", 
                            cohort == "bladder" ~ "mUC",
                            cohort == "Basket" ~ "Basket")) %>% 
  mutate(pep_length = str_length(Mut_peptide)) %>% 
  mutate(HLA_type = str_sub(HLA_allele, start = 1L, end = 5L)) %>% 
  mutate(HLA_num = str_sub(HLA_allele, start = 6L, end = 10L))

all_peptides$Patient <- gsub("Neye","MM",all_peptides$Patient)


# remove space from colnames
names(all_peptides) <- gsub(" ", "_", names(all_peptides))


#### new forginess score 


Pred_Modelling <- Pred_Modelling %>% 
  mutate(RankEL_minus = -RankEL) %>%
  mutate(RankBA_minus = -RankBA) %>% 
  mutate(Stability_minus = -Stability) %>% 
  mutate(NetMHCExp_minus = -NetMHCExp) 
colnames(all_peptides)

all_peptides<- all_peptides %>% 
  mutate(RankEL_minus = -RankEL) %>% 
  mutate(RankBA_minus = -RankBA) %>% 
  mutate(Stability_minus = -Stability) %>% 
  mutate(NetMHCExp_minus = -NetMHCExp) %>% 
  mutate(minus_PropAcidic = -PropAcidic) %>% 
  mutate(minus_PropSmall = -PropSmall) %>% 
  mutate(minus_inst = -Inst) %>% 
  mutate(minus_pI = -pI) 


save(Pred_Modelling,all_peptides, file = "data/04_plotting/Rdata/Final_prep_plotting.Rdata")
write.table(Pred_Modelling, file = "data/04_plotting/Pred_modelling.txt", sep = "\t", col.names = T,quote = F)

