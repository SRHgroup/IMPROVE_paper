#### scripts for adding features to data 

# load library
library(tidyverse)
library(Peptides)
library(dplyr)
library(tidyr)
library(stringr)

# -----------------------------------------------
# Run python feature calculation script first
# python3 bin/python_script/feature_calculations.py --file data/01_data/01_Validated_neoepitopes.txt --dataset "Validated_neoepitopes" --ProgramDir "/Users/annieborch/Documents/programs/" --outfile "data/02_feature_data/02_Validated_neoepitopes_calculated_features.tsv" --TmpDir "/Users/annieborch/Documents/programs/"
# -----------------------------------------------
table(all_peptides$response,all_peptides$Sample)
# load data 
all_peptides <- read.table("data/02_feature_data/02_Validated_neoepitopes_calculated_features.tsv", sep = "\t")
colnames(all_peptides) <- all_peptides[1,]
all_peptides <- all_peptides[-1,]

# source function script 
source("bin/R_script/99_functions.R")


# --------------------------------------------------------------------------------
##                      ## rename HLA binding features and calculate DAI #### 
# ---------------------------------------------------------------------------------
all_peptides <- all_peptides %>% 
  rename("RankEL" = "Mut_MHCrank_EL", 
         "RankBA" = "Mut_MHCrank_BA", 
         "RankEL_wt" = "Norm_MHCrank_EL")

all_peptides$RankEL <- as.numeric(all_peptides$RankEL)
all_peptides$RankEL_wt <- as.numeric(all_peptides$RankEL_wt)

all_peptides <- all_peptides %>% mutate(DAI = RankEL/RankEL_wt )
# --------------------------------------------------------------------------------
##                      ## Add features #### 
# ---------------------------------------------------------------------------------
# make everything numeric 
all_peptides <- 
  all_peptides %>% 
  mutate(Cancer_Driver_Gene = case_when(Cancer_Driver_Gene=="Yes" ~ 1,
                                        Cancer_Driver_Gene=="YES" ~ 1,
                                        Cancer_Driver_Gene=="NO" ~ 0,
                                        Cancer_Driver_Gene=="No" ~ 0)) %>% 
  mutate(Misense_mutation = case_when(Mutation_Consequence =="M" ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(Frameshift_mutation = case_when(Mutation_Consequence =="F" ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(Inframe_deletion_mutation = case_when(Mutation_Consequence =="D" ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(Inframe_insertion = case_when(Mutation_Consequence =="I" ~ 1,
                                      TRUE ~ 0)) %>% 
  mutate(response = case_when(response=="yes" ~ 1,
                              response=="no" ~ 0,
                              response=="1" ~ 1,
                              response=="0" ~ 0))

unique(all_peptides$Cancer_Driver_Gene)


# -------------------------------
# DAI
# --------------------------------





# define Improve and conserved binders 
# ---------------------------------------
# make numeric 
all_peptides[c("RankEL_wt","RankEL")] <- sapply(all_peptides[c("RankEL_wt","RankEL")],as.numeric)
all_peptides$IB_CB<-  all_peptides$RankEL_wt/all_peptides$RankEL
all_peptides <- all_peptides %>% mutate(IB_CB_cat = case_when(IB_CB > 1.2 ~ "CB",
                                              TRUE ~ "IB"))

#-------------------------------------------------------
### add foreignness_score  score 
#-------------------------------------------------------

load("data/02_feature_data/foreignness_score/foreignness_score_output_all.Rdata")

all_peptides <-  all_peptides %>% right_join(.,fs_data)

# -----------------------------------------------------------------------------
# MCP counter features
# -----------------------------------------------------------------------------
all_peptides$Sample_TME <- all_peptides$Sample
#unique(all_peptides$Sample)
#all_peptides$Sample_TME <- gsub("MM909_15_2","MM909_15_1", all_peptides$Sample_TME)
# load MCP counter data
load("data/02_feature_data/TME_and_mcp_counter_raw.Rdata")
unique(all_peptides$Sample)
# peptides are predicted with 22_2

all_peptides <- all_peptides %>% 
  left_join(.,MCP_counter_matrix, by = c("Sample_TME"="Sample")) 

all_peptides$MCPmean <- MCP_counter_mean[all_peptides$Sample_TME]



# -----------------------------------------------------------------------------
# Add Expression for HLA alleles 
# -----------------------------------------------------------------------------
# get HLA type 
all_peptides <- all_peptides %>% 
mutate(HLA_type = str_sub(HLA_allele, start = 1L, end = 5L)) %>% 
  mutate(HLA_num = str_sub(HLA_allele, start = 6L, end = 10L)) 



hla_exp <- HLA_feature %>% 
  mutate(sample_hla = paste(Sample,hugo_symbol,sep = ".")) %>% 
  mutate(HLAexp  = mean_exp) %>% 
  ungroup() %>% 
  dplyr::select(sample_hla,HLAexp) 



all_peptides <- all_peptides %>% 
  mutate(sample_hla = paste(Sample_TME,HLA_type,sep = ".")) %>% 
  left_join(.,hla_exp, by = "sample_hla")

# -----------------------------------------------------------------------------
# CYT 
# --------------------------------------------------
# ----- geometric mean ----- # 
all_peptides <- all_peptides %>% 
  left_join(TME_marker %>% select(Sample,CYT), by = c("Sample_TME"="Sample"))


# RNA data confirm mutations from ibel 
# -----------------------------------------

mutaion_val <- read.csv(file = "data/02_feature_data/Validated_RNA_seq/confirm_rnaseq_all_cohorts.csv", sep = "\t", header = T)


mutaion_val <- mutaion_val %>% 
  mutate(identi_pep_patient = paste(Mut_peptide,HLA_allele,Patient, sep = "_")) %>% 
  distinct(identi_pep_patient, .keep_all = T)


all_peptides <- all_peptides %>% 
  mutate(identi_pep_patient = paste(Mut_peptide,HLA_allele,Patient, sep = "_")) %>% 
  distinct(identi_pep_patient, .keep_all = T)

mutaion_val <- mutaion_val %>% 
  select(identi_pep_patient,rna_confirm,
         rna_var,rna_total,rna_af,rna_coef_100,rna_bin) %>% 
        rename( "ValMutRNACoef"= rna_coef_100)

all_peptides <- all_peptides %>% left_join(.,mutaion_val)


# -----------------------
## Netmhc exp 

net_mhc_exp <- read.table("data/02_feature_data/NetMHCExp/norm_peps.netmhcpanexp.out.form", header = T)
colnames(net_mhc_exp)

net_mhc_exp$cohort <- gsub("basket","Basket",net_mhc_exp$Cohort) 
net_mhc_exp <- net_mhc_exp %>% rename( "NetMHCExp" = "Prediction" ) %>% 
  rename("HLA_allele" = "MHC_mol" ) %>% 
  rename( "response" = "Measure") %>% 
  mutate(norm_pMHC = paste(HLA_allele,Peptide, sep = "_"))
all_peptides <- all_peptides %>%  mutate(norm_pMHC = paste(HLA_allele,Norm_peptide, sep = "_"))


# insert in all_peptide
vec <- net_mhc_exp$NetMHCExp
names(vec) <- net_mhc_exp$norm_pMHC
all_peptides$NetMHCExp <- vec[all_peptides$norm_pMHC]
table(is.na(all_peptides$NetMHCExp))

# get NetExp values where expression is na
# --------------------------------------------
New_exp_values <- read.table("data/02_feature_data/NetMHCExp/New_exp_values.txt")
all_peptides$Expression_Level <- as.numeric(all_peptides$Expression_Level)

New_exp_values <- New_exp_values %>% select(pMHC,gexp_annotated_tcs) %>% distinct(pMHC, .keep_all = T)

all_peptides <- all_peptides %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_")) %>% 
  left_join(.,New_exp_values, by = "pMHC" ) %>% 
  mutate(Expression_Level =
           case_when(is.na(Expression_Level)==T ~ gexp_annotated_tcs, 
                     TRUE ~ Expression_Level)) %>% 
  select(-gexp_annotated_tcs)


# rename mupexi columns:
# -------------------------------------
all_peptides <- all_peptides %>% 
  rename("PrioScore" = "priority_Score", 
         "Expression" = "Expression_Level", 
         "CelPrev" = "cellular_prevalence",
         "VarAlFreq"= "Allele_Frequency")

# make everything numeric 
# --------------------------------------
cols_to_numeric = c('Aro', 'Inst', 'CysRed','RankEL','RankBA','NetMHCExp',
'Expression','SelfSim','Prime','PropHydroAro','HydroCore','HydroAll','pI',
'PropSmall','PropAro','PropBasic','PropAcidic','DAI','Stability','Foreigness',
'CelPrev','PrioScore','VarAlFreq','CYT','HLAexp','Monocytes',
'Tcells','TcellsCD8', 'CytoxLympho','Blinage','NKcells',
'MyeloidDC','Neutrophils','Endothelial' ,'Fibroblasts','MCPmean')


colnames(all_peptides)
all_peptides[cols_to_numeric] <- sapply(all_peptides[cols_to_numeric],as.numeric)

all_peptides$response <- as.factor(all_peptides$response)
which(is.na(all_peptides[cols_to_numeric]))




# if NA put in the mean of the column as RF do not take NA as input
# ----------------------------------------
for(i in 1:ncol(all_peptides)){
  all_peptides[is.na(all_peptides[,i]), i] <- mean(all_peptides[,i], na.rm = TRUE)
}


which(is.na(all_peptides[cols_to_numeric]))

# save for plotting
save(all_peptides, file = "data/04_plotting/Rdata/04_pepides_for_plotting.Rdata")

# Save data for 5-fold CV
# -------------------------------------
write.table(all_peptides, file = "data/02_feature_data/02_2_Validated_neoepitopes_calculated_features_TME.txt", quote = F, sep = '\t')


