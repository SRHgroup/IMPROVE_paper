## MCP counter 

library(tidyverse)
#install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
library(devtools)
#install_github("ebecht/MCPcounter",ref="master", subdir="Source")
source("bin/R_script/99_functions.R")
# ----------------------------
RNA_all_data_raw <- read.table("data/01_data/RNA_data/01_2_RNA_all_data.txt")
# samples from cell lines
RNA_all_data  <- RNA_all_data_raw 

unique(RNA_all_data$Sample)
RNA_all_data$mean_exp <- as.numeric(RNA_all_data$mean_exp )
RNA_all_data <- RNA_all_data %>% group_by(Sample,hugo_symbol) %>% 
  summarise_at(vars(mean_exp), list(mean_exp = mean))

#Making the RNA_data in correct format;
RNA_data_wide <- spread(hugo_symbol, mean_exp, data=RNA_all_data) 
#RNA_data_wide <- RNA_data_wide[!colnames(RNA_data_wide) %in% "hugo_symbol"]
RNA_data_wide <-as.data.frame(RNA_data_wide)
rownames(RNA_data_wide) <- RNA_data_wide$Sample

RNA_data_wide$Sample <-NULL
RNA_data_wide <- t(RNA_data_wide)

#Loading genes and running MCPcounter
require(curl)
genes=read.table(curl:::curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
intersect(rownames(RNA_data_wide),genes[,"HUGO symbols"])

MCP_counter_matrix <- MCPcounter::MCPcounter.estimate(RNA_data_wide, featuresType="HUGO_symbols")

# try with mean --
# ------------------------------------------
MCP_counter_mean <- colMeans(MCP_counter_matrix)

# ------end -----------------------------------------

MCP_counter_matrix <- t(MCP_counter_matrix)
MCP_counter_matrix <- as.data.frame(MCP_counter_matrix)
MCP_counter_matrix$Sample <- rownames(MCP_counter_matrix)
unique(all_peptides$Sample)

MCP_counter_matrix <- MCP_counter_matrix %>% 
  rename("Tcells" = `T cells`,"TcellsCD8" =`CD8 T cells`, "CytoxLympho" = `Cytotoxic lymphocytes`,
         "Blinage"= `B lineage`,"NKcells" = `NK cells`,"Monocytes" =`Monocytic lineage`,
         "MyeloidDC" = `Myeloid dendritic cells`,"Endothelial"=`Endothelial cells` )
unique(all_peptides$Sample)


##################### 

# CYT
#----------
TME_feature <- RNA_all_data %>% 
  filter(hugo_symbol %in% c("GZMA","PRF1")) 

TME_marker <- TME_feature %>%
  #  filter(Sample %in% Sampel_to_select ) %>% 
  mutate(Sampe_hugo = paste(Sample,hugo_symbol, sep = "_")) %>% 
  distinct(Sampe_hugo, .keep_all = T) %>% 
  dplyr::select(Sample,hugo_symbol,mean_exp) %>% 
  filter(hugo_symbol %in% c("PRF1","GZMA")) %>% 
  spread(., key = hugo_symbol, value = mean_exp) 

# calculate geometric mean 
TME_marker$CYT <- apply(TME_marker[ ,c('PRF1', 'GZMA')], 1, gm_mean)

## HLA expression 
# -----------------------------------
HLA_feature <- RNA_all_data %>% 
  filter(hugo_symbol %in% c("HLA-A","HLA-B","HLA-C")) 

  
save(MCP_counter_matrix, HLA_feature,TME_marker,MCP_counter_mean, file = "data/02_feature_data/TME_and_mcp_counter.Rdata")


unique(all_peptides$Sample)

