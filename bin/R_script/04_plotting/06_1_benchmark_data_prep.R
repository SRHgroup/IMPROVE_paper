## benchmark data 
# ---
library(ggplot2)
library(openxlsx)
library(tidyverse)
library(ggsignif)
library(ggbeeswarm)
library(ggrepel)
library(cowplot)
library(ggpubr)

# fore roc curves 
# Loading package
library(caTools)
library(randomForest)
library(pROC)
library(ROCit)
library(caret)
library(ROCR)
library(verification)

# feature calclulations 
# -----------------------------------------------
# Run python feature calculation script first
# python3 bin/python_script/feature_calculations.py --file data/01_data/01_cedar_benchmark_neoepitopes.tsv --dataset "CEDAR_neoepitopes" --ProgramDir "/Users/annieborch/Documents/programs/" --outfile "data/02_feature_data/02_Benchmark_neoepitopes_calculated_features.tsv" --TmpDir "/Users/annieborch/Documents/programs/"
# -----------------------------------------------

CEDAR_all <- read.csv("data/01_data/01_cedar_benchmark_neoepitopes.tsv")

benchmark_data <- read.csv("data/02_feature_data/02_Benchmark_neoepitopes_calculated_features.tsv", header = TRUE, sep = '\t')
#benchmark_data_40 <- read.csv("data/02_feature_data/02_Benchmark_neoepitopes_calculated_features_netmhc40.tsv", header = TRUE, sep = '\t')

benchmark_data$Mut_peptide

max(all_peptides$RankEL)
colnames(benchmark_data)

benchmark_data_CEDAR <- benchmark_data %>% 
  dplyr::rename("RankEL" = "RankEL_4.1")  %>%
  dplyr::rename("RankBA" = "RankBA_4.1")  %>%
  dplyr::rename("DAI" = "DAI_4.1")  %>%
  dplyr::rename("Foreigness" = "Foreignness")  %>% 
   mutate(Patient = "Unknown") %>% 
  filter(RankEL<2) %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_"), .after = "Mut_peptide") %>% 
  distinct(pMHC , .keep_all = T)
table(benchmark_data_CEDAR$response)

write.table(benchmark_data_CEDAR, file = "data/02_feature_data/02_2_Benchmark_neoepitopes_calculated_features.tsv", sep = "\t", quote = F, row.names=F)

benchmark_data_CEDAR_neoepi <- benchmark_data_CEDAR %>% select(HLA_allele,Mut_peptide,response)
write.table(benchmark_data_CEDAR_neoepi, file = "neoepitope_tabels/Neoepitopes_CEDAR_benchmark_data.tsv", sep = "\t", quote = F, row.names=F)

# Predictions
# -----------------------------------------------
# Run python feature calculation script first
#python3 bin/python_script/Predict_immunogenicity.py --file "data/02_feature_data/02_2_Benchmark_neoepitopes_calculated_features.tsv" --model Simple  --outfile "benchmark_data/04_Benchmark_neoepitopes_prediction.tsv" 
# -----------------------------------------------


# -------------------------------------------------

# Prepare for validation tools 


# -------------------------------------------------

# HLA athena
HLA_athena <- benchmark_data_CEDAR %>% mutate(pep = Mut_peptide) %>% select(HLA_allele,pep,Expression)
unique(HLA_athena$HLA_allele)
for (hla in unique(HLA_athena$HLA_allele)) {
  print(hla)
  HLA_athena_selcted <- HLA_athena %>% filter(HLA_allele == hla) %>% select(-HLA_allele)
  peps <- HLA_athena_selcted %>% select(pep) %>% pull() 
  peps <- paste(peps, collapse = ",")
write.table(peps, file = paste0("data/benchmark_comparison/cedar_bench_results/hla_athena_input/",hla,".txt"), quote = F,row.names = F, col.names = F)
}

write.table(HLA_athena, file = "data/benchmark_comparison/HLA_athena_tab.txt", sep = "\t", quote = F, row.names = F)

unique(benchmark_data_prediction$HLA_allele)


# MHCflurry 
Peptides <- benchmark_data_CEDAR %>% select(Mut_peptide) 
write.table(Peptides, file = "results/PaperPlots/Fig6/tables_for_tools/MHC_flurry_peptides.txt", sep = "\n", quote = F, row.names = F)
HLA_alleles_flurry <- benchmark_data_CEDAR %>% select(HLA_allele) %>% pull()
HLA_alleles_flurry
write.table(HLA_alleles_flurry, file = "data/benchmark_comparison/tables_for_tools/HLA_alleles_flurry.txt", sep = "", quote = F, row.names = F)


# for MHCflurry and antigen garnish 

table_cedar <- benchmark_data_CEDAR %>% select(Mut_peptide, HLA_allele)
write.table(table_cedar, file = "data/benchmark_comparison/tables_for_tools/table_cedar.txt", sep = "\t", quote = F, row.names = F)


## correct MHCflurry output:
# ................................
# MHC_flurry_results <- read.table("data/benchmark_comparison/MHCflurry_out.txt", fill = T, sep = ",", header = T)
# MHC_flurry_results <- MHC_flurry_results[!MHC_flurry_results$peptide=="peptide",]
# 
# 
# MHC_flurry_results <- MHC_flurry_results %>%
#   mutate(pMHC = paste(allele,peptide, sep = "_"))
# 
# table_cedar <- table_cedar %>% mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_"))
# table(MHC_flurry_results$pMHC %in% table_cedar$pMHC)
# Missing_mhc_flurry_results <- table_cedar[!table_cedar$pMHC %in% MHC_flurry_results$pMHC,]
# 
# Missing_mhc_flurry_results <- Missing_mhc_flurry_results %>%  select(Mut_peptide, HLA_allele)
# write.table(Missing_mhc_flurry_results, file = "data/benchmark_comparison/tables_for_tools/Missing_mhc_flurry_results.txt", sep = "\t", quote = F, row.names = F)
# 





# for deep net bim 9-mers 
library(stringr)
benchmark_data_CEDAR$pep_length <- str_length(benchmark_data_CEDAR$Mut_peptide)

benchmark_data_CEDAR_net_deep <- benchmark_data_CEDAR %>%
  mutate(mhc = HLA_allele) %>% mutate(sequence = Mut_peptide) %>% 
  filter(pep_length == 9) %>% 
  select(mhc,sequence)
write.table(benchmark_data_CEDAR_net_deep, file = "data/benchmark_data/tables_for_tools/benchmark_data_CEDAR_net_deep.txt", sep = "\t", quote = F, row.names = F)



# deep immu
# ---------
benchmark_data_CEDAR_deep_immu<- benchmark_data_CEDAR %>%
 select(Mut_peptide,HLA_allele)
write.table(benchmark_data_CEDAR_deep_immu, file = "data/benchmark_data/tables_for_tools/benchmark_data_CEDAR_deep_immu.txt", sep = "\t", quote = F, row.names = F, col.names = F)


# --- ittcarf


benchmark_data_CEDAR_ittca <- benchmark_data_CEDAR %>% mutate(fata_name = paste0(">",HLA_allele,"_",Mut_peptide))
new_fasta <- data.frame()
for (row in 1:nrow(benchmark_data_CEDAR_ittca)) {
  fasta_line <- rbind(benchmark_data_CEDAR_ittca$fata_name[row],benchmark_data_CEDAR_ittca$Mut_peptide[row])
  new_fasta <- rbind(new_fasta,fasta_line)
}

write.table(new_fasta, file = "data/benchmark_comparison/tables_for_tools/benchmark_data_CEDAR_ittca.txt", quote = F, row.names = F, col.names = F)



## mixMHCpred HLA list 
# ----------------


unique_hla_formixmhcpred <- table_cedar %>% distinct(HLA_allele)
write.table(unique_hla_formixmhcpred,file = "data/benchmark_comparison/tables_for_tools/unique_hla_formixmhcpred.txt", quote = F)



