#### scripts for adding features to data 

# load library
library(tidyverse)
library(Peptides)
library(dplyr)
library(tidyr)
library(stringr)

# -----------------------------------------------
# Run python feature calculation script first
# python3 bin/python_script/feature_calculations.py --file data/01_data/01_Not_validated_neoepitopes_extra.txt --dataset "Not_Validated_neoepitopes_extra" --ProgramDir "/Users/annieborch/Documents/programs/" --outfile "data/02_feature_data/02_Not_validated_neoepitopes_extra_calculated_features.tsv" --TmpDir "/Users/annieborch/Documents/programs/"
# -----------------------------------------------

# load data 
all_peptides_extra_data <- read.table("data/02_feature_data/02_Not_validated_neoepitopes_extra_calculated_features.tsv", sep = "\t")
colnames(all_peptides_extra_data) <- all_peptides_extra_data[1,]
all_peptides_extra_data <- all_peptides_extra_data[-1,]

unique(all_peptides_extra_data$Patient)

colnames(all_peptides_extra_data)
unique(all_peptides_extra_data$Sample)

# source function script 
source("bin/R_script/99_functions.R")


# --------------------------------------------------------------------------------
##                      ## rename HLA binding features and calculate DAI #### 
# ---------------------------------------------------------------------------------
all_peptides_extra_data <- all_peptides_extra_data %>%
  rename("RankEL" = "Mut_MHCrank_EL",
         "RankBA" = "Mut_MHCrank_BA",
         "RankEL_wt" = "Norm_MHCrank_EL")

all_peptides_extra_data$RankEL <- as.numeric(all_peptides_extra_data$RankEL)
all_peptides_extra_data$RankEL_wt <- as.numeric(all_peptides_extra_data$RankEL_wt)

# -------------------------------
# DAI
# --------------------------------

all_peptides_extra_data <- all_peptides_extra_data %>% mutate(DAI = RankEL/RankEL_wt )
# --------------------------------------------------------------------------------
##                      ## Add features #### 
# ---------------------------------------------------------------------------------
# make everything numeric 
all_peptides_extra_data <- 
  all_peptides_extra_data %>% 
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
                                       TRUE ~ 0))

unique(all_peptides_extra_data$Cancer_Driver_Gene)







# define Improve and conserved binders 
# ---------------------------------------
# make numeric 
all_peptides_extra_data[c("RankEL_wt","RankEL")] <- sapply(all_peptides_extra_data[c("RankEL_wt","RankEL")],as.numeric)
all_peptides_extra_data$IB_CB<-  all_peptides_extra_data$RankEL_wt/all_peptides_extra_data$RankEL
all_peptides_extra_data <- all_peptides_extra_data %>% mutate(IB_CB_cat = case_when(IB_CB > 1.2 ~ "CB",
                                                              TRUE ~ "IB"))


#-------------------------------------------------------
### add foreignness_score  score 
#-------------------------------------------------------

load("data/02_feature_data/foreignness_score/foreignness_score_output_EXTRA_DATA.Rdata")

fs_EXTRA_data <- fs_EXTRA_data %>% 
  select(-nmer) %>% 
  distinct(Mut_peptide, .keep_all = T) %>% 
  rename("Foreigness" = "foreignness_score")

all_peptides_extra_data <-  all_peptides_extra_data %>% left_join(.,fs_EXTRA_data)

# -----------------------------------------------------------------------------
# MCP counter features
# -----------------------------------------------------------------------------
# load MCP counter data
load("data/02_feature_data/TME_and_mcp_counter.Rdata")
unique(all_peptides_extra_data$Sample)

all_peptides_extra_data <- all_peptides_extra_data %>% 
  left_join(.,MCP_counter_matrix) 

all_peptides_extra_data$MCPmean <- MCP_counter_mean[all_peptides_extra_data$Sample]

# -----------------------------------------------------------------------------
# Add Expression for HLA alleles 
# -----------------------------------------------------------------------------
# get HLA type 
all_peptides_extra_data <- all_peptides_extra_data %>% 
  mutate(HLA_type = str_sub(HLA_allele, start = 1L, end = 5L)) %>% 
  mutate(HLA_num = str_sub(HLA_allele, start = 6L, end = 10L)) 


hla_exp <- HLA_feature %>% 
  mutate(sample_hla = paste(Sample,hugo_symbol,sep = ".")) %>% 
  mutate(HLAexp  = mean_exp) %>% 
  dplyr::select(Sample,sample_hla,HLAexp)



all_peptides_extra_data <- all_peptides_extra_data %>% 
  mutate(sample_hla = paste(Sample,HLA_type,sep = ".")) %>% 
  left_join(.,hla_exp)





# -----------------------------------------------------------------------------
# CYT 
# --------------------------------------------------
# ----- geometric mean ----- # 
all_peptides_extra_data <- all_peptides_extra_data %>% 
  left_join(TME_marker %>% select(Sample,CYT), by = "Sample")



# -----------------------
## Netmhc exp 

net_mhc_exp <- read.table("data/02_feature_data/NetMHCExp/norm_peps.netmhcpanexp.out.extra.data.form")
nrow(net_mhc_exp)
colnames(net_mhc_exp) <- net_mhc_exp[1,]
net_mhc_exp <- net_mhc_exp[-1,]
net_mhc_exp <- net_mhc_exp %>% rename( "NetMHCExp" = "%rank_EL" ) %>% 
  mutate(Norm_pMHC = paste(HLA_allele,Norm_peptide, sep = "_"))
all_peptides_extra_data <- all_peptides_extra_data %>% 
  mutate(Norm_pMHC = paste(HLA_allele,Norm_peptide, sep = "_"))

table(net_mhc_exp$Mut_peptide %in% all_peptides_extra_data$Mut_peptide)

# insert in all_peptide
vec <- net_mhc_exp$NetMHCExp
names(vec) <- net_mhc_exp$Norm_pMHC
all_peptides_extra_data$NetMHCExp <- vec[all_peptides_extra_data$Norm_pMHC]
table(is.na(all_peptides_extra_data$NetMHCExp))

# get NetExp values where expression is na
# --------------------------------------------
New_exp_values <- read.table("data/02_feature_data/NetMHCExp/New_exp_values.txt")
all_peptides_extra_data$Expression_Level <- as.numeric(all_peptides_extra_data$Expression_Level)



all_peptides_extra_data <- all_peptides_extra_data %>% 
  left_join(.,New_exp_values) %>% 
  mutate(Expression_Level =
           case_when(is.na(Expression_Level)==T ~ gexp_annotated_tcs, 
                     TRUE ~ Expression_Level)) %>% 
  select(-gexp_annotated_tcs)


# rename mupexi columns:
# -------------------------------------
all_peptides_extra_data <- all_peptides_extra_data %>% 
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


cols_to_numeric %in%  colnames(all_peptides_extra_data) 
all_peptides_extra_data[cols_to_numeric] <- sapply(all_peptides_extra_data[cols_to_numeric],as.numeric)

all_peptides_extra_data$response <- as.factor(all_peptides_extra_data$response)
which(is.na(all_peptides_extra_data[cols_to_numeric]))




# if NA put in the mean of the column as RF do not take NA as input
# ----------------------------------------
# 
for(i in 1:ncol(all_peptides_extra_data)){
  all_peptides_extra_data[is.na(all_peptides_extra_data[,i]), i] <- mean(all_peptides_extra_data[,i], na.rm = TRUE)
}


which(is.na(all_peptides_extra_data[cols_to_numeric]))
# Save data for 5-fold CV
# -------------------------------------
write.table(all_peptides_extra_data, file = "data/02_feature_data/02_2_Not_validated_Validated_neoepitopes_extra_calculated_features_TME.txt", quote = F, sep = '\t')
# ----------------------------------
## Now run predict tool 


#python3 bin/python_script/Predict_immunogenicity.py --file "data/02_feature_data/02_2_Not_validated_Validated_neoepitopes_extra_calculated_features_TME.txt" --model TME_included  --outfile "extra_data/04_1_Not_validated_Validated_neoepitopes_extra_predicted_score_tme_included.txt"

#python3 bin/python_script/Predict_immunogenicity.py --file "data/02_feature_data/02_2_Not_validated_Validated_neoepitopes_extra_calculated_features_TME.txt" --model TME_excluded  --outfile "extra_data/04_1_Not_validated_Validated_neoepitopes_extra_predicted_score_tme_excluded.txt" 



