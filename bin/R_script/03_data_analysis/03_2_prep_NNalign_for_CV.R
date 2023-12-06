# sort peptides for partitioning 
# ---  data without patient data 
library(tidyverse)
library(caret)
# ---------------------------------------------------------
###### Sort peptide for partioning with patient data 
# ---------------------------------------------------------
# make NNalign input

all_peptides <- read.csv('data/02_feature_data/02_2_Validated_neoepitopes_calculated_features_TME.txt', sep ="\t")
colnames(all_peptides)
table(all_peptides$Partition)
table(all_peptides$Cancer_Driver_Gene, all_peptides$response)
# filter out insufficient peptides: 


# add partition 
#  -------------
Partition_dat <- read.csv('data/03_data_for_CV/IMPROVE/03_final_peptide_features_Partition.txt', sep ="\t")
Partition_vec <- Partition_dat$Partition
names(Partition_vec) <- Partition_dat$Patient
all_peptides$Partition <- Partition_vec[all_peptides$Patient]

# ---------- save data
table(all_peptides$response)
table(all_peptides$Partition)
colnames(all_peptides)
# write table 
write.table(all_peptides, file= "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# write a file to tools folder:
#test_file_for_tool_dir <- all_peptides %>% sample_n(.,100) %>% select(Patient,HLA_allele,Mut_peptide,Norm_peptide, Expression,PrioScore,CelPrev,NetMHCExp,Foreigness,CYT,HLAexp,MCPmean)
#write.table(test_file_for_tool_dir, file = "../IMPROVE_git/IMPROVE_tool/data/test_file_for_feature_calculation.tsv", quote = F, sep = "\t")

# Save data for plotting
# -------------------------------------
## NNalign
NNalign_input <- all_peptides %>%
  arrange(Partition) %>% 
  dplyr::select(Mut_peptide,response, Partition,HLA_allele, Patient)
write.table(NNalign_input, file= "data/03_data_for_CV/NNalign/03_NNalign_input_with_p_inf.txt", sep = "\t", quote = F, col.names = F, row.names = F)



