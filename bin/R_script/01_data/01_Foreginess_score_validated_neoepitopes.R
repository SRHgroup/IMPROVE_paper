# module load tools anaconda3/2021.11 gcc/8.2.0 intel/perflibs/2019_update5 R/4.0.0 perl/5.20.1 ncbi-blast/2.12.0+ parallel/20220422

library(magrittr)
library(data.table)
library(antigen.garnish)
library(Biostrings)
library(tidyverse)

load("/home/projects/SRHgroup/projects/MuPeXI_project/scripts/antigen_garnish/03_1_filtered_data.Rdata")


Sys.setenv(AG_DATA_DIR = "/home/projects/SRHgroup/apps/antigen.garnish")
Sys.setenv(ANTIGEN_GARNISH_DIR = "/home/projects/SRHgroup/apps/antigen.garnish/netMHC")
Sys.setenv(HOME = "/home/projects/SRHgroup/apps")



fs_data <- data.frame()

test <- all_peptides$Peptide

for (pep in test) {
fs <- foreignness_score(pep, db = "human") 
fs_data <- rbind(fs_data, c(pep,fs$foreignness_score,fs$nmer))
}
colnames(fs_data) <- c("Peptide","foreignness_score","nmer")

save(fs_data, file = "foreignness_score_output_all.Rdata")

