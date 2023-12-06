# RUN ON C2

# module load tools anaconda3/2021.11 gcc/8.2.0 intel/perflibs/2019_update5 R/4.0.0 perl/5.20.1 ncbi-blast/2.12.0+ parallel/20220422

library(magrittr)
library(data.table)
library(antigen.garnish)
library(Biostrings)
library(tidyverse)

inpu_data <- read.table("/home/projects/SRHgroup/projects/MuPeXI_project/scripts/antigen_garnish/01_All_mupexi_data_extra.txt.sim", sep ="\t", header = TRUE)


Sys.setenv(AG_DATA_DIR = "/home/projects/SRHgroup/apps/antigen.garnish")
Sys.setenv(ANTIGEN_GARNISH_DIR = "/home/projects/SRHgroup/apps/antigen.garnish/netMHC")
Sys.setenv(HOME = "/home/projects/SRHgroup/apps")


fs_EXTRA_data <- data.frame()

test <- inpu_data$Mut_peptide

for (pep in test) {
fs <- foreignness_score(pep, db = "human") 
fs_EXTRA_data <- rbind(fs_EXTRA_data, c(pep,fs$foreignness_score,fs$nmer))
}
colnames(fs_EXTRA_data) <- c("Mut_peptide","foreignness_score","nmer")

save(fs_EXTRA_data, file = "foreignness_score_output_EXTRA_DATA.Rdata")

