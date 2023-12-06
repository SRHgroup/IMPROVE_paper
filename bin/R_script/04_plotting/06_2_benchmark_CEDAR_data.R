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


load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")

# result simple_model 
# ------------------------------

benchmark_data_prediction <- read.table(file = "results/benchmark_data/04_Benchmark_neoepitopes_prediction.tsv", sep = "\t", header = T)
nrow(benchmark_data_prediction)
benchmark_data_prediction <-benchmark_data_prediction %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_")) %>% 
  distinct(pMHC, .keep_all = T)

# merge bencmark results 
# ------------------------------

# IEDB pred
IEDB_pred <- read.table("data/benchmark_comparison/cedar_bench_results/IEDB_prediction.txt", sep = ",", fill = T)

colnames(IEDB_pred) <- IEDB_pred[3,]
IEDB_pred <- IEDB_pred[-c(1:3),]

IEDB_pred <- IEDB_pred %>% distinct(peptide, .keep_all=T)
IEDB_pred <- IEDB_pred %>% rename("IEDB_score" = score)

# join with banchdata 
nrow(benchmark_data_prediction %>% left_join(IEDB_pred, by = c( "Mut_peptide" = "peptide")))
benchmark_data_prediction <- benchmark_data_prediction %>% left_join(IEDB_pred, by = c( "Mut_peptide" = "peptide"))


## iTTCA-RF_results.txt
iTTCA_pred <- read.table("data/benchmark_comparison/cedar_bench_results/iTTCA-RF_results.txt", sep = "\t", fill = T)
iTTCA_pred$iTTCA_score <- gsub("\\%","",iTTCA_pred$V3)
iTTCA_pred$iTTCA_score <- as.numeric(iTTCA_pred$iTTCA_score)
iTTCA_pred <- iTTCA_pred %>% separate(V1, c("HLA_Allele","Mut_peptide"), "_")

iTTCA_pred  <- iTTCA_pred %>% 
  mutate(pMHC = paste(HLA_Allele,Mut_peptide, sep = "_")) %>% 
  distinct(pMHC, .keep_all = T) %>% 
  select(pMHC, iTTCA_score) 

iTTCA_pred$pMHC
benchmark_data_prediction$pMHC
benchmark_data_prediction <- benchmark_data_prediction %>% left_join(iTTCA_pred, by = "pMHC") 



## micMHCpred  mixMHCpred_out.txt
mixMHCpred <- read.table("data/benchmark_comparison/cedar_bench_results/mixMHCpred_out.txt", sep = "\t", fill = T)

colnames(mixMHCpred) <- mixMHCpred[1,]
mixMHCpred <- mixMHCpred[-1,]

mixMHCpred <- mixMHCpred %>% gather(., key = "type", value = "score", -Peptide)

mixMHCpred$HLA_allele <- str_sub(mixMHCpred$type,-5,-1)
library(stringr)
mixMHCpred$HLA_allele_new <- str_replace(mixMHCpred$HLA_allele, pattern = "(.{3})(.*)", replacement = "\\1:\\2")
mixMHCpred$HLA_allele_new <- paste0("HLA-", mixMHCpred$HLA_allele_new)

mixMHCpred <- mixMHCpred %>%
  filter(!type %in% c("Score_bestAllele","BestAllele","%Rank_bestAllele")) %>% 
  mutate(score_type = 
           case_when(grepl("Score",type) ~ "mixMHCpred_Score",
                     grepl("Rank",type) ~ "mixMHCpred_Rank",
                     TRUE ~ "other"))

table(mixMHCpred$pMHC %in% benchmark_data_prediction$pMHC)

mixMHCpred_test <- mixMHCpred %>%
  mutate(pMHC = paste(HLA_allele_new,Peptide, sep = "_")) %>% 
  mutate(ID = paste(pMHC,score_type, sep = "_")) %>% 
  distinct(ID, .keep_all = T)  %>% 
  select(score_type,score,pMHC) %>% 
  spread(., key = score_type, value = score)

mixMHCpred <- mixMHCpred_test[mixMHCpred_test$pMHC %in% benchmark_data_prediction$pMHC,]

benchmark_data_prediction <- benchmark_data_prediction %>% right_join(mixMHCpred, by = "pMHC") 

# netdeepbim
# --
net_deep_bim_results <- read.table(file = "data/benchmark_comparison/cedar_bench_results/net_deep_bim_results_prediction.txt", header = T)

net_deep_bim_results <- net_deep_bim_results %>%
  mutate(pMHC = paste(mhc,sequence, sep = "_")) %>% 
  rename("pred_affinity_netdeepbim" = pred_affinity) %>% 
  rename("pred_immuno_netdeepbim" = pred_immuno) %>% 
  rename("immuno_probability_netdeepbim" = immuno_probability) %>% 
  select(-c(sequence,mhc))


table(net_deep_bim_results$pMHC %in% benchmark_data_prediction$pMHC)
benchmark_data_prediction <- benchmark_data_prediction %>% left_join(.,net_deep_bim_results, by = "pMHC")



# MHCflurry
# ---------------
MHCflurry_results <- read.table(file = "data/benchmark_comparison/cedar_bench_results/MHCflurry_out.txt", header = T, sep = ",")
MHCflurry_results <- MHCflurry_results[!MHCflurry_results$peptide=="peptide",]


MHCflurry_results <- MHCflurry_results %>%
  mutate(pMHC = paste(allele,peptide, sep = "_")) %>% 
  select(-c(allele,peptide)) 


benchmark_data_prediction <- benchmark_data_prediction %>% left_join(., MHCflurry_results, by ="pMHC") 

duplicated(benchmark_data_prediction$pMHC)

benchmark_data_prediction[duplicated(benchmark_data_prediction$pMHC)==T,]

benchmark_data_prediction_filter <- benchmark_data_prediction %>% 
  mutate(IMPROVE = mean_prediction_rf) %>% 
  distinct(pMHC, .keep_all = T)

table(benchmark_data_prediction_filter$response)
# any common peptides

prime_train <- read.xlsx("data/benchmark_comparison/train_set_other_tools/mmc6.xlsx", colNames = T)
colnames(prime_train) <- prime_train[1,]
prime_train <- prime_train[-1,]
prime_train$HLA_allele_new <- str_replace(prime_train$Allele, pattern = "(.{3})(.*)", replacement = "\\1:\\2")
prime_train$HLA_allele_new <- paste0("HLA-", prime_train$HLA_allele_new)


prime_train <- prime_train %>% mutate(pMHC = paste(HLA_allele_new,Mutant, sep = "_"))

table(benchmark_data_prediction_filter$pMHC %in% prime_train$pMHC)
benchmark_data_prediction_filter <- benchmark_data_prediction_filter %>%
  mutate(Prime_new = case_when(benchmark_data_prediction_filter$Mut_peptide %in% prime_train$Mutant ~ NA,
                               TRUE ~ benchmark_data_prediction_filter$Prime))# %>% 

table(benchmark_data_prediction_filter$Mut_peptide %in% prime_train$Mutant)
#mutate(IMPROVE_new = case_when(benchmark_data_prediction_filter$Mut_peptide %in% prime_train$Mutant ~ NA,
 #                            TRUE ~ benchmark_data_prediction_filter$IMPROVE))

mixmhcpred_train_dat <- read.table("data/benchmark_comparison/train_set_other_tools/mmc4.txt", fill = T)
mixmhcpred_train_dat <- mixmhcpred_train_dat[,1:2]
colnames(mixmhcpred_train_dat) <- mixmhcpred_train_dat[2,]
mixmhcpred_train_dat <- mixmhcpred_train_dat[-c(1:2),]


#mhc_flurry_train <- read.xlsx("data/benchmark_comparison/train_set_other_tools/1-s2.0-S2405471220302398-mmc5.xlsx")
#table(benchmark_data_prediction$Mut_peptide %in% mhc_flurry_train$peptide)

write.table(benchmark_data_prediction_filter, file = "data/benchmark_comparison/benchmark_data_prediction_filter.txt", sep = "\t", quote = F)
colnames(benchmark_data_prediction_filter)

benchmark_data_prediction_filter <- benchmark_data_prediction_filter %>% 
  select(IMPROVE,Foreigness,Prime_new,Prime,RankEL,
         IEDB_score,iTTCA_score,immuno_probability_netdeepbim,
         mhcflurry_prediction,mhcflurry_prediction_percentile,
         mixMHCpred_Rank,mixMHCpred_Score,response)

#benchmark_data_prediction_filter <- benchmark_data_prediction_filter[is.na(benchmark_data_prediction_filter$Prime_new)==F,]
benchmark_data_prediction_filter[1:ncol(benchmark_data_prediction_filter)] <- sapply(benchmark_data_prediction_filter[1:ncol(benchmark_data_prediction_filter)],as.numeric)


table(benchmark_data_prediction_filter$response)


auc_df_bench_V2 <- data.frame()

for (col in colnames(benchmark_data_prediction_filter)[1:ncol(benchmark_data_prediction_filter)-1]) {
  print(col)
  df = benchmark_data_prediction_filter[,c(col,"response")]
  df =  na.omit(df)
  print(table(df$response))
  number_pep <- nrow(df)
  y = df["response"] %>% pull()
  x = df[col]  %>% pull()
  rocobj <- roc(y,x, direction = "<")
  #define object to plot and calculate AUC
  AUC <- round(auc(y, x),3)
  # auc 0.1
  pred <- prediction(x,y)
  auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
  auc_01 <- round(auc_01@y.values[[1]],4)
  auc_df_bench_V2 <- rbind(auc_df_bench_V2,c(col,AUC,auc_01,number_pep))
}
colnames(auc_df_bench_V2) <- c("model","AUC","AUC01","number_pep")

#auc_df_bench_V2 <- rbind(auc_df_bench_V2,c("IMPROVE Simple", auc_becnhmark_V2 ,auc_benchmark_0.1_V2))


auc_df_bench_V2$AUC <- as.numeric(auc_df_bench_V2$AUC)
auc_df_bench_V2$AUC01 <- as.numeric(auc_df_bench_V2$AUC01)

table(is.na(benchmark_data_prediction_filter$Prime_new))

save(benchmark_data_prediction_filter,auc_df_bench_V2, file = "results/benchmark_data/plotting_data/CEDAR_data_peformance.Rdata")
