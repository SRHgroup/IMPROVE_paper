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

# result simple_model 
# ------------------------------
benchmark_data_prediction <- read.table(file = "results/benchmark_data/04_Benchmark_neoepitopes_prediction.tsv", sep = "\t", header = T)
#benchmark_data <-  read.table(file = "data/02_feature_data/02_2_Benchmark_neoepitopes_calculated_features.tsv", sep = "\t", header = T)

# merge data 
# benchmark_data_prediction <- benchmark_data %>% 
#   mutate(identity = paste(Patient,HLA_allele,Mut_peptide, sep = "_")) %>% 
#   select(-c(Patient,HLA_allele,Mut_peptide))%>% 
#   left_join(.,benchmark_data_prediction, by = "identity" ) %>% 
#   rename("prediction_rf" = "mean_prediction_rf")

benchmark_data_prediction <- benchmark_data_prediction %>% 
  rename("prediction_rf" = "mean_prediction_rf") %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide , sep = "_")) %>% 
distinct(pMHC, .keep_all = T)

colnames(benchmark_data_prediction)


### find pred score bench mark data 
# -----------------------------------------
pred_df_simple <- read.table(file = "results/5_fold_CV/Simple/pred_df_Simple.txt", sep = " ", header = T)

### Evaluate benchmark performance 
# -----------------------------------------
# Simple model 
# -------------

rocobj_simple <- roc(pred_df_simple$response, pred_df_simple$prediction_rf, direction = "<")
#define object to plot and calculate AUC
auc_simple <- round(auc(pred_df_simple$response, pred_df_simple$prediction_rf),3)
# auc 0.1
pred <- prediction(pred_df_simple$prediction_rf,pred_df_simple$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_simple_0.1 <- round(auc_01@y.values[[1]],4)


# filter rank < 2 
# ----------------------------------------------------------------------------------
benchmark_data_prediction$RankEL <- as.numeric(benchmark_data_prediction$RankEL)
table(benchmark_data_prediction$response)

# Benchmark data 
rocobj_benchmark <- roc(benchmark_data_prediction$response, benchmark_data_prediction$prediction_rf, direction = "<")
#define object to plot and calculate AUC
auc_benchmark <- round(auc(benchmark_data_prediction$response, benchmark_data_prediction$prediction_rf),3)
# auc 0.1
pred <- prediction(benchmark_data_prediction$prediction_rf,benchmark_data_prediction$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_benchmark_0.1 <- round(auc_01@y.values[[1]],4)


# Benchmark data RankEL
rocobj_benchmark_el <- roc(benchmark_data_prediction$response, benchmark_data_prediction$RankEL, direction = ">")
#define object to plot and calculate AUC
auc_benchmark_el <- round(auc(benchmark_data_prediction$response, benchmark_data_prediction$RankEL),3)
# auc 0.1
pred <- prediction(benchmark_data_prediction$RankEL,benchmark_data_prediction$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_benchmark_el_0.1 <- round(auc_01@y.values[[1]],4)

table(benchmark_data_prediction$response)
### BENCHMARK DATA FIGURE 
# ---------------------------------------------------

benchmark_cols <- c("#bcbddc","#fec44f","#fc9272")   #43a2ca

Fig6A <- ggroc(list(rf_simple = rocobj_simple, 
                    rf_benchmark = rocobj_benchmark,
                    EL_bechmark  = rocobj_benchmark_el)) +
   ggtitle(paste0('(AUC IMPROVE simple = ', auc_simple, ')',
                  '(AUC IMPROVE bencmark = ', auc_benchmark , ')\n',
                  ' (AUC RankEL benchmark = ', auc_benchmark_el , ')\n',
                  '(AUC01 IMPROVE simple = ', auc_simple_0.1, ')',
                  '(AUC01 IMPROVE bencmark = ', auc_benchmark_0.1, ')\n',
                  ' (AUC01 RankEL bencmark = ', auc_benchmark_el_0.1,')')) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0),
               color="darkgrey", linetype="dashed") +
   scale_color_manual(values = benchmark_cols ,
                      breaks = c("rf_simple","rf_benchmark","EL_bechmark"), 
                      labels = c("CV","Benchmark","RankEL Benchmark")) +
  theme_bw() +
  labs(color  = "Model", y="Sensitivity (TPR)", x="Specificity")+
  theme(plot.title = element_text(size = 14),legend.position = "bottom",
        legend.text = element_text(size =14 ),
        axis.title = element_text(size = 18 ),
        legend.title = element_text(size = 16 )) + 
  guides(colour = guide_legend(override.aes = list(size=8)))
ggsave(Fig6A, file = "results/PaperPlots/Fig6/Fig6A.pdf", width = 8 , height = 6 )






### compare old run 


benchmark_data_prediction_old <- read.table("/Users/annieborch/Documents/GitHub/Immunugenicity/results/RandomForrest/benchmark_data/benchmark_data_prediction_test.txt", sep = "\t", header = T)
benchmark_data_prediction_old <- benchmark_data_prediction_old %>% mutate(pMHC = paste(HLA_allele,Peptide, sep = "_"))

benchmark_data_prediction[!benchmark_data_prediction$Expression %in% benchmark_data_prediction_old$Expression_Level,]



table(benchmark_data_prediction$Expression %in% benchmark_data_prediction_old$Expression_Level)

