# sort peptides for partitioning 
# ---  data without patient data 
library(tidyverse)
library(caret)
# ---------------------------------------------------------
###### Sort peptide for partioning with patient data 
# ---------------------------------------------------------
# make NNalign input

all_peptides <- read.csv('data/03_1_feature_partition_data/03_2_final_peptide_features_Partition_V2.txt', sep ="\t")
colnames(all_peptides)
table(all_peptides$Partition)
table(all_peptides$Cancer_Driver_Gene, all_peptides$response)
# filter out insufficient peptides: 
#all_peptides <- all_peptides %>% filter(validation=="Sufficient")

# ## Prime 1.0 
# # ----------------
# prime_s <- read.table("/Users/annieborch/Documents/GitHub/Immunugenicity/results/RandomForrest/tabels/pred_df_TME_includedadvance_final_updated_forgines_and_DAI.txt", header = T)
# prime_s <- prime_s %>%
#   mutate(pMHC = paste(HLA_allele,Peptide, sep = "_")) %>%
#   select(pMHC,Score_PRIME,Expression_Level, pI)
# # 
# vec <- prime_s$Score_PRIME
# names(vec) <- prime_s$pMHC
# 
# vec_ex <- prime_s$Expression_Level
# names(vec_ex) <- prime_s$pMHC
# 
# all_peptides <- all_peptides %>%
#   mutate(pMHC = paste(HLA_allele,Mut_peptide)) %>%
#   rename("Prime_1.0" = "Prime" ) %>%
#   rename("Expression_new" = "Expression" )
# 
# all_peptides$Prime <- vec[all_peptides$pMHC]
# all_peptides$Expression <- vec_ex[all_peptides$pMHC]
# 



# all_peptides <- all_peptides %>%
#   mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_")) %>%
#   rename("pI_new" = "pI" )
# table(prime_s$pMHC %in% all_peptides$pMHC)
# # pi 
# vec_pi <- prime_s$pI
# names(vec_pi) <- prime_s$pMHC
# 
# 
# 
# 
# 
# all_peptides$pI <- vec_pi[all_peptides$pMHC]

# # ----------------END
table(all_peptides$response)
# write table 
write.table(all_peptides, file= "data/03_2_data_for_CV/03_3_final_peptide_features_Partition.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Save data for plotting
# -------------------------------------
#save(all_peptides, file = "data/04_plotting/Rdata/04_pepides_for_plotting.Rdata")

## NNalign
NNalign_input <- all_peptides %>%
  arrange(Partition) %>% 
  dplyr::select(Mut_peptide,response, Partition,HLA_allele, Patient)
write.table(NNalign_input, file= "data/03_data_for_CV/NNalign/03_NNalign_input_with_p_inf.txt", sep = "\t", quote = F, col.names = F, row.names = F)



