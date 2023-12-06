# sort peptides for partitioning 
# ---  data without patient data 
library(tidyverse)

# ---------------------------------------------------------
###### Validated data 
# --------------------------------------------------------
load("/Users/annieborch/Documents/GitHub/Immunugenicity/data/01_preprosessing_data/Rdata/01_3_merge_peptides_corrected_prior_score.Rdata")

raw <- all_peptides
load("/Users/annieborch/Documents/GitHub/Immunugenicity/data/03_partitioning_data/Rdata/03_2_filtered_data_model_plotting.Rdata")


# select included pep from raw 
Validated_neoepitopes <-  raw %>% 
  mutate(validation = case_when(raw$identity %in% all_peptides$identity ~ "Sufficient", 
                                 TRUE ~ "Insufficient"))
                                 
colnames(Validated_neoepitopes)
#Validated_neoepitopes <- raw[raw$identity %in% all_peptides$identity,] 
Validated_neoepitopes <- Validated_neoepitopes %>% select(Sample,Patient,HLA_allele,Norm_peptide,Mut_peptide,Mut_MHCrank_EL,Mut_MHCrank_BA,Norm_MHCrank_EL, 
                                 Gene_ID,Transcript_ID,Genomic_Position,Protein_position,Mutation_Consequence,Allele_Frequency,Amino_Acid_Change,
                                 Gene_Symbol,Cancer_Driver_Gene,Expression_Level,priority_Score,cellular_prevalence,
                                 cohort,identity,response,validation)
unique(Validated_neoepitopes$cohort)
Validated_neoepitopes$Norm_peptide <- as.character(Validated_neoepitopes$Norm_peptide)
Validated_neoepitopes$Norm_peptide[Validated_neoepitopes$Norm_peptide=="-----------"] <- NA

Validated_neoepitopes <- Validated_neoepitopes %>% filter(validation=="Sufficient")

write.table(Validated_neoepitopes, file = "data/01_data/01_Validated_neoepitopes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
nrow(Validated_neoepitopes)
# ---------------------------------------------------------
###### Extra data - Not validated
# --------------------------------------------------------

All_mupexi_data_extra <- read.table(file = "/Users/annieborch/Documents/GitHub/Immunugenicity/data/00_raw_data/01_All_mupexi_data_extra.txt.sim", sep = "\t")
colnames(All_mupexi_data_extra) <- All_mupexi_data_extra[1,]
All_mupexi_data_extra <- All_mupexi_data_extra[-1,]
All_mupexi_data_extra <- All_mupexi_data_extra %>% select(Sample,Patient,HLA_allele,Norm_peptide,Mut_peptide,,Mut_MHCrank_EL,Mut_MHCrank_BA,Norm_MHCrank_EL, 
                                  Gene_ID,Transcript_ID,Genomic_Position,Protein_position,Mutation_Consequence,Allele_Frequency,Amino_Acid_Change,
                                 Gene_Symbol,Cancer_Driver_Gene,Expression_Level,priority_Score,cellular_prevalence)
All_mupexi_data_extra$Norm_peptide <- as.character(All_mupexi_data_extra$Norm_peptide)
All_mupexi_data_extra$Norm_peptide[All_mupexi_data_extra$Norm_peptide=="-----------"] <- NA

write.table(All_mupexi_data_extra, "data/01_data/01_Not_validated_neoepitopes_extra.txt", sep = "\t",quote = F)


# ---------------------------------------------------------
###### becnmark data 
# --------------------------------------------------------
benchmark_data_old <- benchmark_data
benchmark_data <- read.csv("/Users/annieborch/Documents/GitHub/Immunugenicity/data/benchmark_data/cedar_benchmark_features_pepX.tsv", header = T, sep = "\t")
benchmark_data <- benchmark_data %>% 
  select(HLA_allele,Norm_peptide,Mut_peptide,Target,NetMHCExp,Foreignness,Expression) %>% 
  rename( "response" = "Target")
benchmark_data$HLA_allele <- sub("(.{7})(:*)", "\\1:\\2", benchmark_data$HLA_allele)


write.table(benchmark_data, "data/01_data/01_cedar_benchmark_neoepitopes.tsv", sep = "\t", quote = F)

# ---------------------------------------------------------
###### forginess score  
# --------------------------------------------------------
load("/Users/annieborch/Documents/GitHub/Immunugenicity/data/02_feature_data/Rdata/foreignness_score_output_all.Rdata")
fs_data <- fs_data %>% 
  rename("Mut_peptide" = "Peptide") %>% 
  select(-nmer) %>% 
  distinct(Mut_peptide, .keep_all = T) %>% 
  rename("Foreigness" = "foreignness_score")
fs_data <- fs_data[fs_data$Mut_peptide %in% Validated_neoepitopes$Mut_peptide,]
save(fs_data,file ="data/02_feature_data/foreignness_score/foreignness_score_output_all.Rdata")
# ---------------------------------------------------------
###### NetMHCexppan  
# --------------------------------------------------------
net_mhc_exp <- read.table("/Users/annieborch/Documents/GitHub/Immunugenicity/data/02_feature_data/txt/norm_peps.netmhcpanexp.out.form",header = T)

net_mhc_exp  <- net_mhc_exp %>%  mutate(norm_pMHC = paste(MHC_mol,Peptide,sep = "_"))
Validated_neoepitopes <- Validated_neoepitopes %>%  mutate(norm_pMHC = paste(HLA_allele,Norm_peptide, sep = "_"))

net_mhc_exp <- net_mhc_exp[net_mhc_exp$norm_pMHC %in% Validated_neoepitopes$norm_pMHC,]
write.table(net_mhc_exp, file = "data/02_feature_data/NetMHCExp/norm_peps.netmhcpanexp.out.form", sep = "\t")
# ---------------------------------------------------------
###### mut_val  
# --------------------------------------------------------
mutaion_val <- read.csv(file = "/Users/annieborch/Documents/GitHub/Immunugenicity/data/00_raw_data/validation_mutation_RNA/confirm_rnaseq_all_cohorts.csv", sep = " ")
nrow(mutaion_val)
mutaion_val$Sample <- sprintf("%02d", as.numeric(mutaion_val$Sample))
mutaion_val <- mutaion_val %>% 
  mutate(Patient = 
           case_when(Cohort=="Bladder" ~ paste("BC",Sample, sep = "-"),
                     Cohort=="Melanoma" ~ paste("Neye",Sample, sep = "-"),
                     Cohort=="Basket" ~ paste("RH",Sample, sep = "-")))

mutaion_val <- mutaion_val %>% 
  mutate(identi_pep_patient = paste(Mut_peptide,HLA_allele,Patient, sep = "_")) %>% 
  distinct(identi_pep_patient, .keep_all = T)


Validated_neoepitopes <- Validated_neoepitopes %>% 
  mutate(identi_pep_patient = paste(Mut_peptide,HLA_allele,Patient, sep = "_")) %>% 
  distinct(identi_pep_patient, .keep_all = T)

mutaion_val <- mutaion_val[mutaion_val$identi_pep_patient %in% Validated_neoepitopes$identi_pep_patient,]

write.table(mutaion_val, file = "data/02_feature_data/Validated_RNA_seq/confirm_rnaseq_all_cohorts.csv", sep ="\t")

# ----------
# new exp values 
# -------------

New_exp_values <- read.csv('/Users/annieborch/Documents/GitHub/Immunugenicity/data/01_preprosessing_data/txt/01_3_merge_peptides_Extra.txt', sep = "\t", header = T)
New_exp_values <-New_exp_values %>% select(HLA_allele,Mut_peptide,gexp_annotated_tcs) %>% mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_"))
Validated_neoepitopes  <- Validated_neoepitopes %>% mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_"))
New_exp_values <- New_exp_values[New_exp_values$pMHC %in% Validated_neoepitopes$pMHC,]

write.table(New_exp_values, "data/02_feature_data/NetMHCExp/New_exp_values.txt")






# Write table to ACT-TIL500 project 

screened_epitopes_melanoma <- Validated_neoepitopes %>% 
  filter(cohort=="melanoma") %>% 
  select(Patient,HLA_allele,Mut_peptide,Norm_peptide,Gene_ID,Gene_Symbol,Transcript_ID,Genomic_Position, Mutation_Consequence, Amino_Acid_Change,response)

write.table(screened_epitopes_melanoma, file = "screened_epitopes_melanoma.txt", quote = F, col.names = T, row.names = F, sep = "\t")



