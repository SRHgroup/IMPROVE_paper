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
#load necessary packages
library(ggplot2)

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
load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
load("data/04_plotting/Survival_data/All_pep_for_survival.Rdata")
source("bin/R_script/99_functions.R")
#cutoff data 
load("data/04_plotting/Survival_data/05_cut_offs.Rdata")




# ------------------------------------------------------------------------------
#                 Figure 5
# ------------------------------------------------------------------------------

# figure 5 A and B
# ------------------------------------------------------------------------------
top_figure_percent(top_number = 20)
top_figure_percent(top_number = 50)
# figure 5 C
# ------------------------------------------------------------------------------

Fig5C_cutoff_figure_TME <- Fig5C_cutoff_figure_TME +   
  geom_vline(xintercept = cut_off_tme_cross) +
  geom_vline(xintercept = cut_off_tme_90) 

Fig5C_cutoff_figure_legend <- get_legend(Fig5C_cutoff_figure)
Fig5C_cutoff_figure <- Fig5C_cutoff_figure + 
  theme(legend.position =  "none") +
  geom_vline(xintercept = cut_off_cross) +
  geom_vline(xintercept = cut_off_90) 


# figure 5 D confusion matrix 
# ------------------------------------------------------------------------------

#model with TME
# ---------------------------

Pred_Modelling$predicted_cat_cross <- ifelse(Pred_Modelling$prediction_rf_tme >cut_off_tme_cross, 1,0)
#Pred_Modelling$predicted_cat <- ifelse(Pred_Modelling$prediction_rf_tme >0.482, 1,0)
Pred_Modelling$predicted_cat_spec <- ifelse(Pred_Modelling$prediction_rf_tme >cut_off_tme_90, 1,0)
Pred_Modelling$predicted_cat_cross <- factor(Pred_Modelling$predicted_cat_cross, levels = c(0,1))
Pred_Modelling$predicted_cat_spec <- factor(Pred_Modelling$predicted_cat_spec, levels = c(0,1))
Pred_Modelling$response <- factor(Pred_Modelling$response, levels = c(0,1))

CM_TME_cross <- confusionMatrix(Pred_Modelling$response, Pred_Modelling$predicted_cat_cross, positive = "1")
CM_TME_spec <- confusionMatrix(Pred_Modelling$response, Pred_Modelling$predicted_cat_spec, positive = "1")


Pred_Modelling <- Pred_Modelling %>% 
  mutate(predition_class_cross = case_when(
    response==1 & predicted_cat_cross==1 ~ "TP",
    response==1 & predicted_cat_cross==0 ~ "FN",
    response==0 & predicted_cat_cross==0 ~ "TN",
    response==0 & predicted_cat_cross==1 ~ "FP"
  )) %>% 
  mutate(predition_class_spec = case_when(
    response==1 & predicted_cat_spec==1 ~ "TP",
    response==1 & predicted_cat_spec==0 ~ "FN",
    response==0 & predicted_cat_spec==0 ~ "TN",
    response==0 & predicted_cat_spec==1 ~ "FP"
  ))

CM_plotting_TME_cross <- Pred_Modelling %>% group_by(response,predicted_cat_cross,predition_class_cross) %>% tally()

CM_plotting_TME_spec <- Pred_Modelling %>% group_by(response,predicted_cat_spec,predition_class_spec) %>% tally()


CM_figure_TME_cross <- ggplot(data =  CM_plotting_TME_cross, 
                              mapping = aes(x = response, y =  predicted_cat_cross)) +
  geom_tile(fill = "white", colour = "white") +
  geom_hline(yintercept = 1.5)+
  geom_vline(xintercept = 1.5)+
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 8) +
  scale_x_discrete(labels = c("1" = "P","0" = "N") ,expand=c(0,0)) +
  scale_y_discrete(labels = c("0" = "N","1" = "P"),expand=c(0,0) ) +
  geom_text(aes(label = predition_class_cross), vjust = -1, size = 8) +
  labs(x = "Actual",y="IMPROVE TME") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 12)) 

CM_figure_TME_spec <- ggplot(data =  CM_plotting_TME_spec, 
                             mapping = aes(x = response, y =  predicted_cat_spec)) +
  geom_tile(fill = "white", colour = "white") +
  geom_hline(yintercept = 1.5)+
  geom_vline(xintercept = 1.5)+
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 8) +
  scale_x_discrete(labels = c("1" = "P","0" = "N") ,expand=c(0,0)) +
  scale_y_discrete(labels = c("0" = "N","1" = "P"),expand=c(0,0) ) +
  geom_text(aes(label = predition_class_spec), vjust = -1, size = 8) +
  labs(x = "Actual",y="IMPROVE TME \n 90% speceficty") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 12)) 




# model without TME
# ---------------------------

Pred_Modelling$predicted_cat_cross <- ifelse(Pred_Modelling$prediction_rf>cut_off_cross, 1,0)
Pred_Modelling$predicted_cat_spec <- ifelse(Pred_Modelling$prediction_rf>cut_off_90, 1,0)
Pred_Modelling$predicted_cat_cross <- factor(Pred_Modelling$predicted_cat_cross, levels = c(0,1))
Pred_Modelling$predicted_cat_spec <- factor(Pred_Modelling$predicted_cat_spec, levels = c(0,1))
Pred_Modelling$response <- factor(Pred_Modelling$response, levels = c(0,1))
CM_rf_cross <- confusionMatrix(Pred_Modelling$response, Pred_Modelling$predicted_cat_cross, positive = "1")
CM_rf_spec <- confusionMatrix(Pred_Modelling$response, Pred_Modelling$predicted_cat_spec, positive = "1")


Pred_Modelling <- Pred_Modelling %>% 
  mutate(predition_class_cross = case_when(
    response==1 & predicted_cat_cross==1 ~ "TP",
    response==1 & predicted_cat_cross==0 ~ "FN",
    response==0 & predicted_cat_cross==0 ~ "TN",
    response==0 & predicted_cat_cross==1 ~ "FP"
  )) %>% 
  mutate(predition_class_spec = case_when(
    response==1 & predicted_cat_spec==1 ~ "TP",
    response==1 & predicted_cat_spec==0 ~ "FN",
    response==0 & predicted_cat_spec==0 ~ "TN",
    response==0 & predicted_cat_spec==1 ~ "FP"
  ))


CM_plotting_rf_cross <- Pred_Modelling %>% group_by(response,predicted_cat_cross,predition_class_cross) %>% tally()
CM_plotting_rf_spec <- Pred_Modelling %>% group_by(response,predicted_cat_spec,predition_class_spec) %>% tally()


CM_figure_rf_cross <- ggplot(data =  CM_plotting_rf_cross, 
                             mapping = aes(x = response, y =  predicted_cat_cross)) +
  geom_tile(fill = "white", colour = "white") +
  geom_hline(yintercept = 1.5)+
  geom_vline(xintercept = 1.5)+
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 8) +
  geom_text(aes(label = predition_class_cross), vjust = -1, size = 8) +
  scale_x_discrete(labels = c("1" = "P","0" = "N") ,expand=c(0,0)) +
  scale_y_discrete(labels = c("0" = "N","1" = "P") ,expand=c(0,0)) +
  labs(x = "Actual",y="IMPROVE") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 12)) 

CM_figure_rf_spec <- ggplot(data =  CM_plotting_rf_spec, 
                            mapping = aes(x = response, y =  predicted_cat_spec)) +
  geom_tile(fill = "white", colour = "white") +
  geom_hline(yintercept = 1.5)+
  geom_vline(xintercept = 1.5)+
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 8) +
  geom_text(aes(label = predition_class_spec), vjust = -1, size = 8) +
  scale_x_discrete(labels = c("1" = "P","0" = "N") ,expand=c(0,0)) +
  scale_y_discrete(labels = c("0" = "N","1" = "P") ,expand=c(0,0)) +
  labs(x = "Actual",y="IMPROVE \n 90% speceficty") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 12)) 




# EL rank model 
# ---------------------------

Pred_Modelling$EL_cat_spec <- ifelse(Pred_Modelling$RankEL < 0.05, 1,0)
Pred_Modelling$EL_cat_cross <- ifelse(Pred_Modelling$RankEL < 0.45, 1,0)
Pred_Modelling$EL_cat_spec <- factor(Pred_Modelling$EL_cat_spec, levels = c(0,1))
Pred_Modelling$EL_cat_cross <- factor(Pred_Modelling$EL_cat_cross, levels = c(0,1))
Pred_Modelling$response <- factor(Pred_Modelling$response, levels = c(0,1))


Pred_Modelling <- Pred_Modelling %>% 
  mutate(predition_class_EL_cross = case_when(
    response==1 & EL_cat_cross==1 ~ "TP",
    response==1 & EL_cat_cross==0 ~ "FN",
    response==0 & EL_cat_cross==0 ~ "TN",
    response==0 & EL_cat_cross==1 ~ "FP"
  )) %>% 
  mutate(predition_class_EL_spec = case_when(
    response==1 & EL_cat_spec==1 ~ "TP",
    response==1 & EL_cat_spec==0 ~ "FN",
    response==0 & EL_cat_spec==0 ~ "TN",
    response==0 & EL_cat_spec==1 ~ "FP"
  ))


CM_plotting_EL_cross <- Pred_Modelling %>% group_by(response,EL_cat_cross,predition_class_EL_cross) %>% tally()

CM_plotting_EL_spec <- Pred_Modelling %>% group_by(response,EL_cat_spec,predition_class_EL_spec) %>% tally()

CM_figure_EL_cross <- ggplot(data =  CM_plotting_EL_cross, 
                             mapping = aes(x = response, y =  EL_cat_cross)) +
  geom_tile(fill = "white", colour = "white") +
  geom_hline(yintercept = 1.5)+
  geom_vline(xintercept = 1.5)+
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 8) +
  geom_text(aes(label = predition_class_EL_cross), vjust = -1, size =8) +
  scale_x_discrete(labels = c("1" = "P","0" = "N") ,expand=c(0,0)) +
  scale_y_discrete(labels = c("0" = "N","1" = "P") ,expand=c(0,0)) +
  labs(x = "Actual",y="RankEL") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 12)) 

CM_figure_EL_spec <- ggplot(data =  CM_plotting_EL_spec, 
                            mapping = aes(x = response, y =  EL_cat_spec)) +
  geom_tile(fill = "white", colour = "white") +
  geom_hline(yintercept = 1.5)+
  geom_vline(xintercept = 1.5)+
  geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, size = 8) +
  geom_text(aes(label = predition_class_EL_spec), vjust = -1, size =8) +
  scale_x_discrete(labels = c("1" = "P","0" = "N") ,expand=c(0,0)) +
  scale_y_discrete(labels = c("0" = "N","1" = "P") ,expand=c(0,0)) +
  labs(x = "Actual",y="RankEL predicted \n 90% speceficty") +
  theme_bw() + theme(legend.position = "none", 
                     axis.title = element_text(size = 18), 
                     axis.text = element_text(size = 12)) 





# calculate MCC 
# ---------------------------------------------------------------
# from CM_plotting_EL
TN = CM_plotting_EL_cross$n[CM_plotting_EL_cross$predition_class_EL_cross=="TN"]
FP = CM_plotting_EL_cross$n[CM_plotting_EL_cross$predition_class_EL_cross=="FP"]
TP = CM_plotting_EL_cross$n[CM_plotting_EL_cross$predition_class_EL_cross=="TP"]
FN = CM_plotting_EL_cross$n[CM_plotting_EL_cross$predition_class_EL_cross=="FN"]

num1 = as.numeric((TP+FP)*(TP+FN))
num2 = as.numeric((TN+FP)*(TN+FN))
MCC_EL_cross = (TP*TN-FP*FN)/sqrt(num1*num2)
MCC_EL_cross = round(MCC_EL_cross,2)

ACC_EL_cross = (TP+TN)/(TP+TN+FP+FN)
ACC_EL_cross = round(ACC_EL_cross,2)

TN = CM_plotting_EL_spec$n[CM_plotting_EL_spec$predition_class_EL_spec=="TN"]
FP = CM_plotting_EL_spec$n[CM_plotting_EL_spec$predition_class_EL_spec=="FP"]
TP = CM_plotting_EL_spec$n[CM_plotting_EL_spec$predition_class_EL_spec=="TP"]
FN = CM_plotting_EL_spec$n[CM_plotting_EL_spec$predition_class_EL_spec=="FN"]

num1 = as.numeric((TP+FP)*(TP+FN))
num2 = as.numeric((TN+FP)*(TN+FN))
MCC_EL_spec = (TP*TN-FP*FN)/sqrt(num1*num2)
MCC_EL_spec = round(MCC_EL_spec,2)

ACC_EL_spec = (TP+TN)/(TP+TN+FP+FN)
ACC_EL_spec = round(ACC_EL_spec,2)

CM_plotting_EL_cross$model <- "RankEL"
colnames(CM_plotting_EL_cross) <- c("response", "predicted","class","number","model")

CM_plotting_EL_cross$spec <- "RankEL"
colnames(CM_plotting_EL_spec) <- c("response", "predicted","class","number","model")


# from CM_plotting_rf
# ---------------------------------
TN = CM_plotting_rf_cross$n[CM_plotting_rf_cross$predition_class_cross=="TN"]
FP = CM_plotting_rf_cross$n[CM_plotting_rf_cross$predition_class_cross=="FP"]
TP = CM_plotting_rf_cross$n[CM_plotting_rf_cross$predition_class_cross=="TP"]
FN = CM_plotting_rf_cross$n[CM_plotting_rf_cross$predition_class_cross=="FN"]

FP+TP
num1 = as.numeric((TP+FP)*(TP+FN))
num2 = as.numeric((TN+FP)*(TN+FN))
MCC_rf_cross = (TP*TN-FP*FN)/sqrt(num1*num2)
MCC_rf_cross = round(MCC_rf_cross,2)
ACC_rf_cross = (TP+TN)/(TP+TN+FP+FN)
ACC_rf_cross = round(ACC_rf_cross,2)


TN = CM_plotting_rf_spec$n[CM_plotting_rf_spec$predition_class_spec=="TN"]
FP = CM_plotting_rf_spec$n[CM_plotting_rf_spec$predition_class_spec=="FP"]
TP = CM_plotting_rf_spec$n[CM_plotting_rf_spec$predition_class_spec=="TP"]
FN = CM_plotting_rf_spec$n[CM_plotting_rf_spec$predition_class_spec=="FN"]

FP+TP
num1 = as.numeric((TP+FP)*(TP+FN))
num2 = as.numeric((TN+FP)*(TN+FN))
MCC_rf_spec = (TP*TN-FP*FN)/sqrt(num1*num2)
MCC_rf_spec= round(MCC_rf_spec,2)
ACC_rf_spec = (TP+TN)/(TP+TN+FP+FN)
ACC_rf_spec = round(ACC_rf_spec,2)

CM_plotting_rf_cross$model <- "RF"
colnames(CM_plotting_rf_cross) <- c("response", "predicted","class","number","model")

CM_plotting_rf_spec$model <- "RF"
colnames(CM_plotting_rf_spec) <- c("response", "predicted","class","number","model")


# from CM_plotting_TME 
# -------------------------------
TN = CM_plotting_TME_cross$n[CM_plotting_TME_cross$predition_class_cross=="TN"]
FP = CM_plotting_TME_cross$n[CM_plotting_TME_cross$predition_class_cross=="FP"]
TP = CM_plotting_TME_cross$n[CM_plotting_TME_cross$predition_class_cross=="TP"]
FN = CM_plotting_TME_cross$n[CM_plotting_TME_cross$predition_class_cross=="FN"]

num1 = as.numeric((TP+FP)*(TP+FN))
num2 = as.numeric((TN+FP)*(TN+FN))
MCC_rf_tme_cross = (TP*TN-FP*FN)/sqrt(num1*num2)
MCC_rf_tme_cross = round(MCC_rf_tme_cross,2)
FP+TP
ACC_rf_tme_cross = (TP+TN)/(TP+TN+FP+FN)
ACC_rf_tme_cross = round(ACC_rf_tme_cross ,2)


TN = CM_plotting_TME_spec$n[CM_plotting_TME_spec$predition_class_spec=="TN"]
FP = CM_plotting_TME_spec$n[CM_plotting_TME_spec$predition_class_spec=="FP"]
TP = CM_plotting_TME_spec$n[CM_plotting_TME_spec$predition_class_spec=="TP"]
FN = CM_plotting_TME_spec$n[CM_plotting_TME_spec$predition_class_spec=="FN"]

num1 = as.numeric((TP+FP)*(TP+FN))
num2 = as.numeric((TN+FP)*(TN+FN))
MCC_rf_tme_spec = (TP*TN-FP*FN)/sqrt(num1*num2)
MCC_rf_tme_spec = round(MCC_rf_tme_spec,2)
FP+TP
ACC_rf_tme_spec = (TP+TN)/(TP+TN+FP+FN)
ACC_rf_tme_spec = round(ACC_rf_tme_spec ,2)


CM_plotting_TME_cross$model <- "RF_TME"
colnames(CM_plotting_TME_cross) <- c("response", "predicted","class","number","model")

CM_plotting_TME_spec$model <- "RF_TME"
colnames(CM_plotting_TME_spec) <- c("response", "predicted","class","number","model")



CM_matrix_all <- bind_rows(CM_plotting_EL, CM_plotting_rf, CM_plotting_TME)
write.table(CM_matrix_all, file = "results/tabels/CM_matrix_all.txt", sep = "\t", quote = F, col.names = T, row.names = F)
#



# ---------------------------------------------------------------------
# Survival figures 
# ---------------------------------------------------------------------

load("data/04_plotting/Survival_data/All_pep_for_survival.Rdata")
All_samples_clinical<- read.csv('data/04_plotting/Survival_data/All_samples_clinical.txt', sep = "\t", header = T)
All_samples_clinical$Patient <- gsub("Neye","MM",All_samples_clinical$Patient)
All_samples_clinical$Patient <- gsub("BC","mUC",All_samples_clinical$Patient)

# add all data together 
All_samples_clinical_cohort <- all_peptides %>% dplyr::select(Patient ,cohort) %>% distinct(Patient, .keep_all = T)
All_samples_clinical <- All_samples_clinical %>% right_join(.,All_samples_clinical_cohort, by = "Patient")
colnames(Pred_Modelling_total_all_peptide)

# add sample in neodata 
sample_vec <- all_peptides$Sample
names(sample_vec) <- all_peptides$Patient
Pred_Modelling_total_all_peptide$Sample <- sample_vec[Pred_Modelling_total_all_peptide$Patient]
# #timepoint 1

n2 <- all_peptides %>% select(Sample,Patient) %>%
  distinct(Sample, .keep_all = T)  %>% 
  group_by(Patient) %>%
  add_tally() %>%
  filter(n==2)

n2$timepoint <- str_sub(n2$Sample,-1,-1)
n2 <- n2 %>% filter(timepoint %in% c("A","1","2"))

n1 <- all_peptides %>% select(Sample,Patient) %>%
  group_by(Patient) %>%
  distinct(Sample, .keep_all = T)  %>% add_tally() %>%
  filter(n==1) %>% mutate(timepoint = "1")


Sample_to_incude <- bind_rows(n1,n2) %>% filter(!Sample %in% c("PATIEN08","MM909_15_1","MM909_24_1","MM909_22_1")) %>% pull(Sample) #"MM909_24_1","MM909_22_1",
Sample_to_incude <- c(Sample_to_incude,"MM909_15_2") 
All_samples_clinical <- All_samples_clinical %>% filter(!Patient %in% c("RH-08","MM-22","MM-24"))  #,"MM-22" "MM-24",

Pred_Modelling_total_all_peptide <- Pred_Modelling_total_all_peptide %>% 
  filter(Sample %in% Sample_to_incude)  %>% 
  filter(Expression >= 0.01) %>%  # 
  mutate(ID = paste(Patient, HLA_allele,Mut_peptide)) %>% 
  distinct(ID, .keep_all = T)


All_samples_clinical <- Pred_Modelling_total_all_peptide %>%
  filter(RankEL< rankel_cutoff_cross) %>% 
  group_by(Patient) %>%
  tally(n = "prediction_EL_05") %>%
  right_join(.,All_samples_clinical)



All_samples_clinical <- Pred_Modelling_total_all_peptide %>%
  filter(prediction_rf > as.numeric(cut_off_cross)) %>%
  group_by(Patient) %>%
  tally(n = "prediction_rf") %>%
  right_join(.,All_samples_clinical)


All_samples_clinical <- Pred_Modelling_total_all_peptide %>%
  filter(prediction_rf_tme > as.numeric(cut_off_tme_cross)) %>%
  #  filter(prediction_rf_tme > as.numeric(0.536)) %>%
  group_by(Patient) %>%
  tally(n = "prediction_rf_TME") %>%
  right_join(.,All_samples_clinical)


All_samples_clinical <- Pred_Modelling %>% 
  filter(response==1) %>% 
  group_by(Patient) %>% 
  tally(n = "number_response") %>% 
  right_join(.,All_samples_clinical) 

All_samples_clinical <- Pred_Modelling %>% 
  group_by(Patient) %>% 
  tally(n = "number_screened") %>% 
  right_join(.,All_samples_clinical)

All_samples_clinical  <- All_samples_clinical %>% mutate(fraction_response = number_response/number_screened)



All_samples_clinical$prediction_EL_05[is.na(All_samples_clinical$prediction_EL_05)==T] <- 0
All_samples_clinical$prediction_rf_TME[is.na(All_samples_clinical$prediction_rf_TME)==T] <- 0
All_samples_clinical$prediction_rf[is.na(All_samples_clinical$prediction_rf)==T] <- 0
All_samples_clinical$number_response[is.na(All_samples_clinical$number_response)==T] <- 0
#All_samples_clinical$mean_score_tme[is.na(All_samples_clinical$mean_score_tme_pos)==T] <- 0
All_samples_clinical$fraction_response[is.na(All_samples_clinical$fraction_response)==T] <- 0

All_samples_clinical$PFS.Time <- gsub(",",".", All_samples_clinical$PFS.Time)
All_samples_clinical$PFS.Time <- as.numeric(All_samples_clinical$PFS.Time)
All_samples_clinical$OS.Time <- gsub(",",".", All_samples_clinical$OS.Time)
All_samples_clinical$OS.Time <- as.numeric(All_samples_clinical$OS.Time)


quantile_prediction_rf_TME <-as.numeric(quantile(All_samples_clinical$prediction_rf_TME))
quantile_prediction_rf <-as.numeric(quantile(All_samples_clinical$prediction_rf))
quantile_prediction_EL_05 <-as.numeric(quantile(All_samples_clinical$prediction_EL_05))
quantile_number_responses <- as.numeric(quantile(All_samples_clinical$number_response))
quantile_mut_load <- as.numeric(quantile(All_samples_clinical$mut_load))


All_samples_clinical_subset <- All_samples_clinical %>% 
  mutate(pred_neo_model_cut_TME = 
           case_when(prediction_rf_TME <= quantile_prediction_rf_TME[2] ~ "low",
                     prediction_rf_TME > quantile_prediction_rf_TME[2] & prediction_rf_TME < quantile_prediction_rf_TME[3]  ~"medium_low",
                     prediction_rf_TME > quantile_prediction_rf_TME[3] & prediction_rf_TME < quantile_prediction_rf_TME[4]  ~"medium_high",
                     prediction_rf_TME >= quantile_prediction_rf_TME[4] ~ "high")) %>% 
  mutate(pred_neo_model_cut = 
           case_when(prediction_rf <= quantile_prediction_rf[2] ~ "low",
                     prediction_rf > quantile_prediction_rf[2] & prediction_rf < quantile_prediction_rf[3]  ~"medium_low",
                     prediction_rf > quantile_prediction_rf[3] & prediction_rf < quantile_prediction_rf[4]  ~"medium_high",
                     prediction_rf >= quantile_prediction_rf[4] ~ "high")) %>% 
  # mutate(Mut_load_cut = 
  #          case_when(mut_load > median(mut_load) ~ "high",
  #                    TRUE ~"low")) %>% 
  mutate(pred_EL_05 = 
           case_when(prediction_EL_05 <= quantile_prediction_EL_05[2] ~ "low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[3]  ~"medium_low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[4]  ~"medium_high",
                     prediction_EL_05 >= quantile_prediction_EL_05[4] ~ "high")) %>% 
  mutate(number_response_cut = 
           case_when(number_response <= quantile_number_responses[2] ~ "low",
                     number_response  > quantile_number_responses[2] & number_response  < quantile_number_responses[3]  ~"medium_low",
                     number_response  >quantile_number_responses[2] & number_response  < quantile_number_responses[4]  ~"medium_high",
                     number_response  >= quantile_number_responses[4] ~ "high")) %>% 
  mutate(mut_load_cut = 
           case_when(mut_load <= quantile_mut_load[2] ~ "low",
                     mut_load  > quantile_mut_load[2] & mut_load  < quantile_mut_load[3]  ~"medium_low",
                     mut_load > quantile_mut_load[2] & mut_load  < quantile_mut_load[4]  ~"medium_high",
                     mut_load >= quantile_mut_load[4] ~ "high")) 


All_samples_clinical_subset$pred_neo_model_cut_TME <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_neo_model_cut <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_EL_05 <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$number_response_cut <- factor(All_samples_clinical_subset$number_response_cut, levels = c("high", "medium_high","medium_low","low")) 
All_samples_clinical_subset$mut_load_cut <- factor(All_samples_clinical_subset$mut_load_cut, levels = c("high", "medium_high","medium_low","low")) 


save(All_samples_clinical_subset, file = "data/04_plotting/All_samples_clinical_subset.Rdata")

pfs <- Surv(time =  All_samples_clinical_subset$PFS.Time, event =  All_samples_clinical_subset$PFS.Event)
os <- Surv(time =  All_samples_clinical_subset$OS.Time, event =  All_samples_clinical_subset$OS.Event) 

All_samples_clinical_subset$pred_neo_model_cut_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("low", "high"))
HR_and_pval(col = All_samples_clinical_subset$pred_neo_model_cut_high_low,
            surv_object = pfs)
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                        surv_object = pfs)


PFS_pred_neo_model_cut_survfit <- survfit(pfs ~ pred_neo_model_cut, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_survplot <- surv_plot(PFS_pred_neo_model_cut_survfit, 
                                             legend_lab_names =c("high", "medium high ","medium low","low"),
                                             legend_position = c(.9, .95),
                                             plot_title = "IMPROVE without TME",
                                             # xlimits = c(0,30),
                                             pcoord = c(26, 1),
                                             ylab_name = 'PFS',
                                             p_value = surv_list$pval_survival,
                                             palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                             HR = surv_list$HR_survival)

surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                        surv_object = os)
OS_pred_neo_model_cut_survfit<- survfit(os ~ pred_neo_model_cut, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_survplot <- surv_plot(OS_pred_neo_model_cut_survfit, 
                                            legend_lab_names =c("high", "medium high ","medium low","low"),
                                            legend_position = c(.9, .95),
                                            plot_title = "IMPROVE without TME",
                                            #  xlimits = c(0,30),
                                            pcoord = c(26, 1),
                                            p_value = surv_list$pval_survival,
                                            palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                            HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_neo_model_cut_TME_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = pfs)
PFS_pred_neo_model_cut_TME_survfit <- survfit(pfs ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_TME_survplot <- surv_plot(PFS_pred_neo_model_cut_TME_survfit, 
                                                 legend_lab_names =c("high", "medium high ","medium low","low"),
                                                 legend_position = c(.9, .95),
                                                 plot_title = "IMPROVE with TME",
                                                 #    xlimits = c(0,30),
                                                 pcoord = c(26, 1),
                                                 p_value = surv_list$pval_survival,
                                                 ylab_name = 'PFS',
                                                 palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                 HR = surv_list$HR_survival)

surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = os)
OS_pred_neo_model_cut_TME_survfit <- survfit(os ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_TME_survplot <- surv_plot(OS_pred_neo_model_cut_TME_survfit, 
                                                legend_lab_names =c("high", "medium high ","medium low","low"),
                                                legend_position = c(.9, .95),
                                                plot_title = "IMPROVE with TME",
                                                p_value = surv_list$pval_survival,
                                                #  xlimits = c(0,30),
                                                pcoord = c(27, 1),
                                                palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_EL_05_high_low <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("low", "high"))
surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                        surv_object = pfs)
PFS_pred_EL_05_cut_survfit <- survfit(pfs ~ pred_EL_05, data = All_samples_clinical_subset)
PFS_pred_EL_05_cut_survplot <- surv_plot(PFS_pred_EL_05_cut_survfit,
                                         legend_lab_names =c("high", "medium high ","medium low","low"),
                                         legend_position = c(.9, .95),
                                         p_value = surv_list$pval_survival,
                                         plot_title = "RankEL",
                                         #  xlimits = c(0,30),
                                         pcoord = c(26, 1),
                                         ylab_name = 'PFS',
                                         palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                         HR = surv_list$HR_survival)
surv_list <- HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                         surv_object = os)
OS_pred_EL_05_cut_survfit <- survfit(os ~ pred_EL_05, data = All_samples_clinical_subset)
OS_pred_EL_05_cut_survplot <- surv_plot(OS_pred_EL_05_cut_survfit,
                                        legend_lab_names =c("high", "medium high ","medium low","low"),
                                        legend_position = c(.9, .95),
                                        p_value = surv_list$pval_survival,
                                        plot_title = "RankEL",
                                        #  xlimits = c(0,30),
                                        pcoord = c(26, 1),
                                        palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                        HR = surv_list$HR_survival)




Fig5F <- ggarrange(PFS_pred_EL_05_cut_survplot$plot, PFS_pred_neo_model_cut_survplot$plot, PFS_pred_neo_model_cut_TME_survplot$plot,
                   PFS_pred_EL_05_cut_survplot$table, PFS_pred_neo_model_cut_survplot$table, PFS_pred_neo_model_cut_TME_survplot$table,
                   OS_pred_EL_05_cut_survplot$plot, OS_pred_neo_model_cut_survplot$plot, OS_pred_neo_model_cut_TME_survplot$plot,
                   OS_pred_EL_05_cut_survplot$table, OS_pred_neo_model_cut_survplot$table,OS_pred_neo_model_cut_TME_survplot$table,
                   heights = c(1, 0.4),
                   ncol = 3, nrow = 4,
                   align = "v")

## save figure 

pdf(file = "results/PaperPlots/Fig5/Fig5_ALL.pdf", width = 18, height =20 )
ggdraw() +
  draw_plot(p_rankel,        .0, .66, .31, 0.12)+
  draw_plot(Fig5C_cutoff_figure,        .335, .66, .31, 0.12)+
  draw_plot(Fig5C_cutoff_figure_TME,    .67, .66, .31, 0.12) +
  draw_plot(Fig5C_cutoff_figure_legend, .46, 0.82, 0.1, 0.12) +
  
  draw_plot(Fig5F,         .0, .0, 1, 0.52) +
  draw_plot_label(c(paste0("MCC:",MCC_EL_cross),paste0("MCC:",MCC_EL_spec),
                    paste0("MCC:",MCC_rf_cross),paste0("MCC:",MCC_rf_spec),
                    paste0("MCC:",MCC_rf_tme_cross),paste0("MCC:",MCC_rf_tme_spec)), 
                  x = c(0.025,0.2,0.37,0.55,.7,0.88), 
                  y = c(0.66,0.66,0.66,0.66,0.66,0.66), size = 16) + 
  draw_plot(CM_figure_EL_cross, .0, .5, .156, 0.14)+
  draw_plot(CM_figure_EL_spec, .16, .5, .16, 0.14)+
  draw_plot(CM_figure_rf_cross, .335, .5, .155, 0.14) +
  draw_plot(CM_figure_rf_spec, .5, .5, .16, 0.14) +
  draw_plot(CM_figure_TME_cross, .67, .5, .155, 0.14) +
  draw_plot(CM_figure_TME_spec, .84, .5, .16, 0.14) +
  draw_plot_label(c("A","B","C","D","E"), x = c(0,0.5,.0,.0,.0), y = c(1,1,0.80,0.65,0.5), size = 22)  


dev.off()




#---------------------------------------------------------
#       melanoma cohort
#---------------------------------------------------------

All_samples_clinical_subset <- All_samples_clinical %>% filter(cohort=="Melanoma")
quantile_prediction_rf_TME <-as.numeric(quantile(All_samples_clinical_subset$prediction_rf_TME))
quantile_prediction_rf <-as.numeric(quantile(All_samples_clinical_subset$prediction_rf))
quantile_prediction_EL_05 <-as.numeric(quantile(All_samples_clinical_subset$prediction_EL_05))
quantile_number_responses <- as.numeric(quantile(All_samples_clinical_subset$number_response))
quantile_mut_load <- as.numeric(quantile(All_samples_clinical$mut_load))


All_samples_clinical_subset <- All_samples_clinical_subset %>% 
  mutate(pred_neo_model_cut_TME = 
           case_when(prediction_rf_TME <= quantile_prediction_rf_TME[2] ~ "low",
                     prediction_rf_TME > quantile_prediction_rf_TME[2] & prediction_rf_TME < quantile_prediction_rf_TME[3]  ~"medium_low",
                     prediction_rf_TME > quantile_prediction_rf_TME[3] & prediction_rf_TME < quantile_prediction_rf_TME[4]  ~"medium_high",
                     prediction_rf_TME >= quantile_prediction_rf_TME[4] ~ "high")) %>% 
  mutate(pred_neo_model_cut = 
           case_when(prediction_rf <= quantile_prediction_rf[2] ~ "low",
                     prediction_rf > quantile_prediction_rf[2] & prediction_rf < quantile_prediction_rf[3]  ~"medium_low",
                     prediction_rf > quantile_prediction_rf[3] & prediction_rf < quantile_prediction_rf[4]  ~"medium_high",
                     prediction_rf >= quantile_prediction_rf[4] ~ "high")) %>% 
  # mutate(Mut_load_cut = 
  #          case_when(mut_load > median(mut_load) ~ "high",
  #                    TRUE ~"low")) %>% 
  mutate(pred_EL_05 = 
           case_when(prediction_EL_05 <= quantile_prediction_EL_05[2] ~ "low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[3]  ~"medium_low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[4]  ~"medium_high",
                     prediction_EL_05 >= quantile_prediction_EL_05[4] ~ "high")) %>% 
  mutate(number_response_cut = 
           case_when(number_response <= quantile_number_responses[2] ~ "low",
                     number_response  > quantile_number_responses[2] & number_response  < quantile_number_responses[3]  ~"medium_low",
                     number_response  >quantile_number_responses[2] & number_response  < quantile_number_responses[4]  ~"medium_high",
                     number_response  >= quantile_number_responses[4] ~ "high")) %>% 
  mutate(mut_load_cut = 
           case_when(mut_load  <= quantile_mut_load[2] ~ "low",
                     mut_load   > quantile_mut_load[2] & mut_load  < quantile_mut_load[3]  ~"medium_low",
                     mut_load   >quantile_mut_load[2] & mut_load  < quantile_mut_load[4]  ~"medium_high",
                     mut_load   >= quantile_mut_load[4] ~ "high")) 

All_samples_clinical_subset$pred_neo_model_cut_TME <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_neo_model_cut <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_EL_05 <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$mut_load <- factor(All_samples_clinical_subset$mut_load, levels = c("high", "medium_high","medium_low","low")) 


pfs <- Surv(time =  All_samples_clinical_subset$PFS.Time, event =  All_samples_clinical_subset$PFS.Event)
os <- Surv(time =  All_samples_clinical_subset$OS.Time, event =  All_samples_clinical_subset$OS.Event) 

All_samples_clinical_subset$pred_neo_model_cut_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                        surv_object = pfs)

PFS_pred_neo_model_cut_survfit <- survfit(pfs ~ pred_neo_model_cut, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_survplot <- surv_plot(PFS_pred_neo_model_cut_survfit, 
                                             legend_lab_names =c("high", "medium high ","medium low","low"),
                                             legend_position = c(.9, .95),
                                             plot_title = "IMPROVE without TME",
                                             # xlimits = c(0,30),
                                             pcoord = c(26, 1),
                                             ylab_name = 'PFS',
                                             p_value = surv_list$pval_survival,
                                             palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                             HR = surv_list$HR_survival)

surv_list <- HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                         surv_object = os)

OS_pred_neo_model_cut_survfit<- survfit(os ~ pred_neo_model_cut, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_survplot <- surv_plot(OS_pred_neo_model_cut_survfit, 
                                            legend_lab_names =c("high", "medium high ","medium low","low"),
                                            legend_position = c(.9, .95),
                                            plot_title = "IMPROVE without TME",
                                            #  xlimits = c(0,30),
                                            pcoord = c(36, 1),
                                            p_value = surv_list$pval_survival,
                                            palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                            HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_neo_model_cut_TME_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = pfs)
PFS_pred_neo_model_cut_TME_survfit <- survfit(pfs ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_TME_survplot <- surv_plot(PFS_pred_neo_model_cut_TME_survfit, 
                                                 legend_lab_names =c("high", "medium high ","medium low","low"),
                                                 legend_position = c(.9, .95),
                                                 plot_title = "IMPROVE with TME",
                                                 #    xlimits = c(0,30),
                                                 pcoord = c(26, 1),
                                                 p_value = surv_list$pval_survival,
                                                 ylab_name = 'PFS',
                                                 palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                 HR = surv_list$HR_survival)

surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = os)
OS_pred_neo_model_cut_TME_survfit <- survfit(os ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_TME_survplot <- surv_plot(OS_pred_neo_model_cut_TME_survfit, 
                                                legend_lab_names =c("high", "medium high ","medium low","low"),
                                                legend_position = c(.9, .95),
                                                plot_title = "IMPROVE with TME",
                                                p_value = surv_list$pval_survival,
                                                #  xlimits = c(0,30),
                                                pcoord = c(36, 1),
                                                palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_EL_05_high_low <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("low", "high"))
surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low  ,
                        surv_object = pfs)
PFS_pred_EL_05_cut_survfit <- survfit(pfs ~ pred_EL_05, data = All_samples_clinical_subset)
PFS_pred_EL_05_cut_survplot <- surv_plot(PFS_pred_EL_05_cut_survfit, 
                                         legend_lab_names =c("high", "medium high ","medium low","low"),
                                         legend_position = c(.9, .95),
                                         p_value = surv_list$pval_survival,
                                         plot_title = "RankEL",
                                         #  xlimits = c(0,30),
                                         pcoord = c(26, 1),
                                         ylab_name = 'PFS',
                                         palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                         HR = surv_list$HR_survival)
surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                        surv_object = os)
OS_pred_EL_05_cut_survfit <- survfit(os ~ pred_EL_05, data = All_samples_clinical_subset)
OS_pred_EL_05_cut_survplot <- surv_plot(OS_pred_EL_05_cut_survfit, 
                                        legend_lab_names =c("high", "medium high ","medium low","low"),
                                        legend_position = c(.9, .95),
                                        p_value = surv_list$pval_survival,
                                        plot_title = "RankEL",
                                        #  xlimits = c(0,30),
                                        pcoord = c(36, 1),
                                        palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                        HR = surv_list$HR_survival)


plot <- ggarrange(PFS_pred_EL_05_cut_survplot$plot, PFS_pred_neo_model_cut_survplot$plot, PFS_pred_neo_model_cut_TME_survplot$plot,
                  PFS_pred_EL_05_cut_survplot$table, PFS_pred_neo_model_cut_survplot$table, PFS_pred_neo_model_cut_TME_survplot$table,
                  OS_pred_EL_05_cut_survplot$plot, OS_pred_neo_model_cut_survplot$plot, OS_pred_neo_model_cut_TME_survplot$plot,
                  OS_pred_EL_05_cut_survplot$table, OS_pred_neo_model_cut_survplot$table,OS_pred_neo_model_cut_TME_survplot$table,
                  heights = c(1, 0.4),
                  ncol = 3, nrow = 4,
                  align = "v") 


plot_melanoma <- annotate_figure(plot, top = text_grob("Melanoma cohort", 
                                                       color = "black",  size = 24))
ggsave(plot, file = "results/PaperPlots/Fig5/SupplementaryFig5D_melanoma.pdf",width = 22, height = 11)




# All_samples_clinical_subset$pred_neo_model_cut_TME <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("low","high") )
# All_samples_clinical_subset$pred_neo_model_cut <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("low","high"))
# All_samples_clinical_subset$pred_EL_05 <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("low","high") )
#  
# ggforest(coxph(pfs ~ pred_neo_model_cut, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# ggforest(coxph(os ~ pred_neo_model_cut, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# 
# ggforest(coxph(pfs ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'pfs Hazard Ratio')
# ggforest(coxph(os ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# 
# ggforest(coxph(pfs ~ pred_EL_05, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# ggforest(coxph(os ~ pred_EL_05, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# 



#  ---------------------------------------------------------
#       bladder cohort
# ---------------------------------------------------------

All_samples_clinical_subset <- All_samples_clinical %>% filter(cohort=="mUC")
quantile_prediction_rf_TME <-as.numeric(quantile(All_samples_clinical_subset$prediction_rf_TME))
quantile_prediction_rf <-as.numeric(quantile(All_samples_clinical_subset$prediction_rf))
quantile_prediction_EL_05 <-as.numeric(quantile(All_samples_clinical_subset$prediction_EL_05))
quantile_number_responses <- as.numeric(quantile(All_samples_clinical_subset$number_response))
quantile_mut_load <- as.numeric(quantile(All_samples_clinical$mut_load))

All_samples_clinical_subset <- All_samples_clinical_subset %>% 
  mutate(pred_neo_model_cut_TME = 
           case_when(prediction_rf_TME <= quantile_prediction_rf_TME[2] ~ "low",
                     prediction_rf_TME > quantile_prediction_rf_TME[2] & prediction_rf_TME <= quantile_prediction_rf_TME[3]  ~"medium_low",
                     prediction_rf_TME > quantile_prediction_rf_TME[3] & prediction_rf_TME <= quantile_prediction_rf_TME[4]  ~"medium_high",
                     prediction_rf_TME > quantile_prediction_rf_TME[4] ~ "high")) %>% 
  mutate(pred_neo_model_cut = 
           case_when(prediction_rf <= quantile_prediction_rf[2] ~ "low",
                     prediction_rf > quantile_prediction_rf[2] & prediction_rf <= quantile_prediction_rf[3]  ~"medium_low",
                     prediction_rf > quantile_prediction_rf[3] & prediction_rf <= quantile_prediction_rf[4]  ~"medium_high",
                     prediction_rf > quantile_prediction_rf[4] ~ "high")) %>% 
  # mutate(Mut_load_cut = 
  #          case_when(mut_load > median(mut_load) ~ "high",
  #                    TRUE ~"low")) %>% 
  mutate(pred_EL_05 = 
           case_when(prediction_EL_05 <= quantile_prediction_EL_05[2] ~ "low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[3]  ~"medium_low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[4]  ~"medium_high",
                     prediction_EL_05 >= quantile_prediction_EL_05[4] ~ "high")) %>% 
  mutate(number_response_cut = 
           case_when(number_response <= quantile_number_responses[2] ~ "low",
                     number_response  > quantile_number_responses[2] & number_response  < quantile_number_responses[3]  ~"medium_low",
                     number_response  >quantile_number_responses[2] & number_response  < quantile_number_responses[4]  ~"medium_high",
                     number_response  >= quantile_number_responses[4] ~ "high")) %>% 
  mutate(mut_load_cut = 
           case_when(mut_load  < quantile_mut_load[2] ~ "low",
                     mut_load   > quantile_mut_load[2] & mut_load  < quantile_mut_load[3]  ~"medium_low",
                     mut_load   >quantile_mut_load[2] & mut_load  < quantile_mut_load[4]  ~"medium_high",
                     mut_load   >= quantile_mut_load[4] ~ "high"))

All_samples_clinical_subset$pred_neo_model_cut_TME <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_neo_model_cut <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_EL_05 <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$mut_load <- factor(All_samples_clinical_subset$mut_load, levels = c("high", "medium_high","medium_low","low")) 


pfs <- Surv(time =  All_samples_clinical_subset$PFS.Time, event =  All_samples_clinical_subset$PFS.Event)
os <- Surv(time =  All_samples_clinical_subset$OS.Time, event =  All_samples_clinical_subset$OS.Event) 

All_samples_clinical_subset$pred_neo_model_cut_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                        surv_object = pfs)

PFS_pred_neo_model_cut_survfit <- survfit(pfs ~ pred_neo_model_cut, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_survplot <- surv_plot(PFS_pred_neo_model_cut_survfit, 
                                             legend_lab_names =c("high", "medium high ","medium low","low"),
                                             legend_position = c(.9, .95),
                                             plot_title = "IMPROVE without TME",
                                             # xlimits = c(0,30),
                                             pcoord = c(26, 1),
                                             ylab_name = 'PFS',
                                             p_value = surv_list$pval_survival,
                                             palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                             HR = surv_list$HR_survival)

surv_list <- HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                         surv_object = os)

OS_pred_neo_model_cut_survfit<- survfit(os ~ pred_neo_model_cut, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_survplot <- surv_plot(OS_pred_neo_model_cut_survfit, 
                                            legend_lab_names =c("high", "medium high ","medium low","low"),
                                            legend_position = c(.9, .95),
                                            plot_title = "IMPROVE without TME",
                                            #  xlimits = c(0,30),
                                            pcoord = c(26, 1),
                                            p_value = surv_list$pval_survival,
                                            palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                            HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_neo_model_cut_TME_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = pfs)
PFS_pred_neo_model_cut_TME_survfit <- survfit(pfs ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_TME_survplot <- surv_plot(PFS_pred_neo_model_cut_TME_survfit, 
                                                 legend_lab_names =c("high", "medium high ","medium low","low"),
                                                 legend_position = c(.9, .95),
                                                 plot_title = "IMPROVE with TME",
                                                 #    xlimits = c(0,30),
                                                 pcoord = c(26, 1),
                                                 p_value = surv_list$pval_survival,
                                                 ylab_name = 'PFS',
                                                 palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                 HR = surv_list$HR_survival)

surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = os)
OS_pred_neo_model_cut_TME_survfit <- survfit(os ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_TME_survplot <- surv_plot(OS_pred_neo_model_cut_TME_survfit, 
                                                legend_lab_names =c("high", "medium high ","medium low","low"),
                                                legend_position = c(.9, .95),
                                                plot_title = "IMPROVE with TME",
                                                p_value = surv_list$pval_survival,
                                                #  xlimits = c(0,30),
                                                pcoord = c(27, 1),
                                                palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_EL_05_high_low <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("low", "high"))
surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                        surv_object = pfs)
PFS_pred_EL_05_cut_survfit <- survfit(pfs ~ pred_EL_05, data = All_samples_clinical_subset)
PFS_pred_EL_05_cut_survplot <- surv_plot(PFS_pred_EL_05_cut_survfit, 
                                         legend_lab_names =c("high", "medium high ","medium low","low"),
                                         legend_position = c(.9, .95),
                                         p_value = surv_list$pval_survival,
                                         plot_title = "RankEL",
                                         #  xlimits = c(0,30),
                                         pcoord = c(26, 1),
                                         ylab_name = 'PFS',
                                         palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                         HR = surv_list$HR_survival)
surv_list <- HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                         surv_object = os)
OS_pred_EL_05_cut_survfit <- survfit(os ~ pred_EL_05, data = All_samples_clinical_subset)
OS_pred_EL_05_cut_survplot <- surv_plot(OS_pred_EL_05_cut_survfit, 
                                        legend_lab_names =c("high", "medium high ","medium low","low"),
                                        legend_position = c(.9, .95),
                                        p_value = surv_list$pval_survival,
                                        plot_title = "RankEL",
                                        #  xlimits = c(0,30),
                                        pcoord = c(26, 1),
                                        palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                        HR = surv_list$HR_survival)

plot <- ggarrange(PFS_pred_EL_05_cut_survplot$plot, PFS_pred_neo_model_cut_survplot$plot, PFS_pred_neo_model_cut_TME_survplot$plot,
                  PFS_pred_EL_05_cut_survplot$table, PFS_pred_neo_model_cut_survplot$table, PFS_pred_neo_model_cut_TME_survplot$table,
                  OS_pred_EL_05_cut_survplot$plot, OS_pred_neo_model_cut_survplot$plot, OS_pred_neo_model_cut_TME_survplot$plot,
                  OS_pred_EL_05_cut_survplot$table, OS_pred_neo_model_cut_survplot$table,OS_pred_neo_model_cut_TME_survplot$table,
                  heights = c(1, 0.4),
                  ncol = 3, nrow = 4,
                  align = "v") 

plot_mUC <- annotate_figure(plot, top = text_grob("mUC cohort", 
                                                  color = "black",  size = 24))
ggsave(plot, file = "results/PaperPlots/Fig5/SupplementaryFig5E_bladder.pdf",width = 22, height = 11)


#  # 
# All_samples_clinical_subset$pred_neo_model_cut_TME <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("low","high") )
# All_samples_clinical_subset$pred_neo_model_cut <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("low","high"))
# All_samples_clinical_subset$pred_EL_05 <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("low","high") )
#  
# ggforest(coxph(pfs ~ pred_neo_model_cut, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# ggforest(coxph(os ~ pred_neo_model_cut, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
#  
# ggforest(coxph(pfs ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'pfs Hazard Ratio')
# ggforest(coxph(os ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
#  
# ggforest(coxph(pfs ~ pred_EL_05, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
# ggforest(coxph(os ~ pred_EL_05, data = All_samples_clinical_subset), data = All_samples_clinical_subset, main = 'OS Hazard Ratio')
#  
#  

#  ---------------------------------------------------------
#       Basket cohort
# ---------------------------------------------------------

All_samples_clinical_subset <- All_samples_clinical %>% filter(cohort=="Basket")
quantile_prediction_rf_TME <-as.numeric(quantile(All_samples_clinical_subset$prediction_rf_TME))
quantile_prediction_rf <-as.numeric(quantile(All_samples_clinical_subset$prediction_rf))
quantile_prediction_EL_05 <-as.numeric(quantile(All_samples_clinical_subset$prediction_EL_05))
quantile_number_responses <- as.numeric(quantile(All_samples_clinical_subset$number_response))
quantile_mut_load <- as.numeric(quantile(All_samples_clinical$mut_load))

All_samples_clinical_subset <- All_samples_clinical_subset %>% 
  mutate(pred_neo_model_cut_TME = 
           case_when(prediction_rf_TME <= quantile_prediction_rf_TME[2] ~ "low",
                     prediction_rf_TME > quantile_prediction_rf_TME[2] & prediction_rf_TME < quantile_prediction_rf_TME[3]  ~"medium_low",
                     prediction_rf_TME > quantile_prediction_rf_TME[3] & prediction_rf_TME < quantile_prediction_rf_TME[4]  ~"medium_high",
                     prediction_rf_TME >= quantile_prediction_rf_TME[4] ~ "high")) %>% 
  mutate(pred_neo_model_cut = 
           case_when(prediction_rf <= quantile_prediction_rf[2] ~ "low",
                     prediction_rf > quantile_prediction_rf[2] & prediction_rf < quantile_prediction_rf[3]  ~"medium_low",
                     prediction_rf > quantile_prediction_rf[3] & prediction_rf < quantile_prediction_rf[4]  ~"medium_high",
                     prediction_rf >= quantile_prediction_rf[4] ~ "high")) %>% 
  # mutate(Mut_load_cut = 
  #          case_when(mut_load > median(mut_load) ~ "high",
  #                    TRUE ~"low")) %>% 
  mutate(pred_EL_05 = 
           case_when(prediction_EL_05 <= quantile_prediction_EL_05[2] ~ "low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[3]  ~"medium_low",
                     prediction_EL_05 > quantile_prediction_EL_05[2] & prediction_EL_05 < quantile_prediction_EL_05[4]  ~"medium_high",
                     prediction_EL_05 >= quantile_prediction_EL_05[4] ~ "high")) %>% 
  mutate(number_response_cut = 
           case_when(number_response <= quantile_number_responses[2] ~ "low",
                     number_response  > quantile_number_responses[2] & number_response  < quantile_number_responses[3]  ~"medium_low",
                     number_response  >quantile_number_responses[2] & number_response  < quantile_number_responses[4]  ~"medium_high",
                     number_response  > quantile_number_responses[4] ~ "high")) %>% 
  mutate(mut_load_cut = 
           case_when(mut_load  <= quantile_mut_load[2] ~ "low",
                     mut_load   > quantile_mut_load[2] & mut_load  < quantile_mut_load[3]  ~"medium_low",
                     mut_load   >quantile_mut_load[2] & mut_load  < quantile_mut_load[4]  ~"medium_high",
                     mut_load   >= quantile_mut_load[4] ~ "high"))

All_samples_clinical_subset$pred_neo_model_cut_TME <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_neo_model_cut <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("high", "medium_high","medium_low","low"))
All_samples_clinical_subset$pred_EL_05 <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("high", "medium_high","medium_low","low")) 
All_samples_clinical_subset$mut_load <- factor(All_samples_clinical_subset$mut_load, levels = c("high", "medium_high","medium_low","low")) 


pfs <- Surv(time =  All_samples_clinical_subset$PFS.Time, event =  All_samples_clinical_subset$PFS.Event)
os <- Surv(time =  All_samples_clinical_subset$OS.Time, event =  All_samples_clinical_subset$OS.Event) 


All_samples_clinical_subset$pred_neo_model_cut_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                        surv_object = pfs)

PFS_pred_neo_model_cut_survfit <- survfit(pfs ~ pred_neo_model_cut, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_survplot <- surv_plot(PFS_pred_neo_model_cut_survfit, 
                                             legend_lab_names =c("high", "medium high ","medium low","low"),
                                             legend_position = c(.9, .95),
                                             plot_title = "IMPROVE without TME",
                                             xlimits = c(0,30),
                                             pcoord = c(9, 1),
                                             ylab_name = 'PFS',
                                             p_value = surv_list$pval_survival,
                                             palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                             HR = surv_list$HR_survival)

surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_high_low,
                        surv_object = os)

OS_pred_neo_model_cut_survfit<- survfit(os ~ pred_neo_model_cut, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_survplot <- surv_plot(OS_pred_neo_model_cut_survfit, 
                                            legend_lab_names =c("high", "medium high ","medium low","low"),
                                            legend_position = c(.9, .95),
                                            plot_title = "IMPROVE without TME",
                                            xlimits = c(0,30),
                                            pcoord = c(9, 1),
                                            p_value = surv_list$pval_survival,
                                            palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                            HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_neo_model_cut_TME_high_low <- factor(All_samples_clinical_subset$pred_neo_model_cut_TME, levels = c("low", "high"))
surv_list <-HR_and_pval(col =All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = pfs)
PFS_pred_neo_model_cut_TME_survfit <- survfit(pfs ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
PFS_pred_neo_model_cut_TME_survplot <- surv_plot(PFS_pred_neo_model_cut_TME_survfit, 
                                                 legend_lab_names =c("high", "medium high ","medium low","low"),
                                                 legend_position = c(.9, .95),
                                                 plot_title = "IMPROVE with TME",
                                                 xlimits = c(0,30),
                                                 pcoord = c(9, 1),
                                                 p_value = surv_list$pval_survival,
                                                 ylab_name = 'PFS',
                                                 palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                 HR = surv_list$HR_survival)

surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
                        surv_object = os)
OS_pred_neo_model_cut_TME_survfit <- survfit(os ~ pred_neo_model_cut_TME, data = All_samples_clinical_subset)
OS_pred_neo_model_cut_TME_survplot <- surv_plot(OS_pred_neo_model_cut_TME_survfit, 
                                                legend_lab_names =c("high", "medium high ","medium low","low"),
                                                legend_position = c(.9, .95),
                                                plot_title = "IMPROVE with TME",
                                                p_value = surv_list$pval_survival,
                                                xlimits = c(0,30),
                                                pcoord = c(9, 1),
                                                palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                                HR = surv_list$HR_survival)

All_samples_clinical_subset$pred_EL_05_high_low <- factor(All_samples_clinical_subset$pred_EL_05, levels = c("low", "high"))
surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                        surv_object = pfs)
PFS_pred_EL_05_cut_survfit <- survfit(pfs ~ pred_EL_05, data = All_samples_clinical_subset)
PFS_pred_EL_05_cut_survplot <- surv_plot(PFS_pred_EL_05_cut_survfit, 
                                         legend_lab_names =c("high", "medium high ","medium low","low"),
                                         legend_position = c(.9, .95),
                                         p_value = surv_list$pval_survival,
                                         plot_title = "RankEL",
                                         xlimits = c(0,30),
                                         pcoord = c(9, 1),
                                         ylab_name = 'PFS',
                                         palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                         HR = surv_list$HR_survival)
surv_list <-HR_and_pval(col = All_samples_clinical_subset$pred_EL_05_high_low ,
                        surv_object = os)
OS_pred_EL_05_cut_survfit <- survfit(os ~ pred_EL_05, data = All_samples_clinical_subset)
OS_pred_EL_05_cut_survplot <- surv_plot(OS_pred_EL_05_cut_survfit, 
                                        legend_lab_names =c("high", "medium high ","medium low","low"),
                                        legend_position = c(.9, .95),
                                        p_value = surv_list$pval_survival,
                                        plot_title = "RankEL",
                                        xlimits = c(0,30),
                                        pcoord = c(9, 1),
                                        palette = c('#6a3d9a', '#e31a1c','#1f78b4', '#33a02c'),
                                        HR = surv_list$HR_survival)




plot <- ggarrange(PFS_pred_EL_05_cut_survplot$plot, PFS_pred_neo_model_cut_survplot$plot, PFS_pred_neo_model_cut_TME_survplot$plot,
                  PFS_pred_EL_05_cut_survplot$table, PFS_pred_neo_model_cut_survplot$table, PFS_pred_neo_model_cut_TME_survplot$table,
                  OS_pred_EL_05_cut_survplot$plot, OS_pred_neo_model_cut_survplot$plot, OS_pred_neo_model_cut_TME_survplot$plot,
                  OS_pred_EL_05_cut_survplot$table, OS_pred_neo_model_cut_survplot$table,OS_pred_neo_model_cut_TME_survplot$table,
                  heights = c(1, 0.4),
                  ncol = 3, nrow = 4,
                  align = "v") 

plot_Basket <- annotate_figure(plot, top = text_grob("Basket cohort", 
                                                     color = "black",  size = 24))
ggsave(plot, file = "results/PaperPlots/Fig5/SupplementaryFig5F_basket.pdf",width = 22, height = 11)






pdf(file = "results/PaperPlots/Fig5/SupplementaryFig5_ALL.pdf", width = 20, height =31 )
ggdraw() +
  draw_plot(plot_Basket ,  .0, .66, 1, 0.33)+
  draw_plot(plot_melanoma, .0, .33, 1, 0.33) +
  draw_plot(plot_mUC,      .0, 0,   1, 0.33) +
  draw_plot_label(c("A","B","C"), x = c(0,0,0), y = c(0.98,0.65,.32), size = 26)  


dev.off()


