## MCP counter 

library(tidyverse)
#install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
library(devtools)
#install_github("ebecht/MCPcounter",ref="master", subdir="Source")
# library(openxlsx)
# library(devtools)
# library(ggsignif)
# library(ggplot2)
# library(ggpubr)
# library(ggbeeswarm)
# library(readxl)

#load("data/01_preprosessing_data/Rdata/01_3_merge_peptides.Rdata")
# ----------------------------
RNA_all_data <- read.table("data/01_data/RNA_data/01_2_RNA_all_data.txt")

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
MCP_counter_matrix <- t(MCP_counter_matrix)
MCP_counter_matrix <- as.data.frame(MCP_counter_matrix)
MCP_counter_matrix$Sample <- rownames(MCP_counter_matrix)



MCP_counter_matrix <- MCP_counter_matrix %>% 
  rename("Tcells" = `T cells`,"TcellsCD8" =`CD8 T cells`, "CytoxLympho" = `Cytotoxic lymphocytes`,
         "Blinage"= `B lineage`,"NKcells" = `NK cells`,"Monocytes" =`Monocytic lineage`,
         "MyeloidDC" = `Myeloid dendritic cells`,"Endothelial"=`Endothelial cells` )

# M22_MCP <- MCP_counter_matrix %>% 
#   gather(., key = "feature", value = "count" , -Sample) %>% 
#   filter(Sample %in% c("MM909_22_1","MM909_22_2")) %>% 
#   group_by(feature) %>% summarize(count = mean(count, na.rm=TRUE)) %>% 
#   mutate(Sample = "MM909_22_1") %>% 
#   spread(., key = feature , value = count)



#MCP_counter_matrix <- MCP_counter_matrix %>% filter(!Sample %in%c("MM909_22_1","MM909_22_2")) %>% bind_rows(.,M22_MCP)

##################### 

TME_feature <- RNA_all_data %>% 
  filter(hugo_symbol %in% c("GZMA","PRF1","HLA-A","HLA-B","HLA-C")) 


# # Sample 22 mean 
# MM22_TME <- TME_exp_dat %>% 
#   filter(Sample %in% c("MM909_22_1","MM909_22_2")) %>% 
#   group_by(hugo_symbol) %>% summarize(mean_exp = mean(mean_exp, na.rm=TRUE)) %>% 
#   mutate(Sample = "MM909_22_1") 
#   
# TME_feature <- TME_feature %>% filter(!Sample %in% c("MM909_22_1","MM909_22_2")) %>% 
#   bind_rows(., MM22_TME)


save(MCP_counter_matrix, TME_feature, file = "data/02_feature_data/TME_and_mcp_counter.Rdata")
