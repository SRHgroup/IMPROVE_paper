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
benchmark_data_prediction$RankEL
benchmark_data_prediction <- benchmark_data_prediction %>% 
  rename("prediction_rf" = "mean_prediction_rf") %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide , sep = "_")) %>% 
distinct(pMHC, .keep_all = T) %>% 
  mutate(RankEL_minus = -RankEL) 

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
rocobj_benchmark_el <- roc(benchmark_data_prediction$response, benchmark_data_prediction$RankEL_minus, direction = "<")
#define object to plot and calculate AUC
auc_benchmark_el <- round(auc(benchmark_data_prediction$response, benchmark_data_prediction$RankEL_minus),3)
# auc 0.1
pred <- prediction(benchmark_data_prediction$RankEL_minus,benchmark_data_prediction$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_benchmark_el_0.1 <- round(auc_01@y.values[[1]],4)

min(benchmark_data_prediction$Expression)

### BENCHMARK DATA FIGURE 
# ---------------------------------------------------

benchmark_cols <- c("#bcbddc","#fec44f","#fc9272")   #43a2ca

Fig6A <- ggroc(list(rf_simple = rocobj_simple, 
                    rf_benchmark = rocobj_benchmark,
                    EL_bechmark  = rocobj_benchmark_el)) +
   # ggtitle(paste0('(AUC IMPROVE simple CV = ', auc_simple, ')',
   #                '(AUC IMPROVE simple bencmark = ', auc_benchmark , ')\n',
   #                ' (AUC RankEL benchmark = ', auc_benchmark_el , ')\n',
   #                '(AUC01 IMPROVE simple = ', auc_simple_0.1, ')',
   #                '(AUC01 IMPROVE bencmark = ', auc_benchmark_0.1, ')\n',
   #                ' (AUC01 RankEL bencmark = ', auc_benchmark_el_0.1,')')) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0),
               color="darkgrey", linetype="dashed") +
   scale_color_manual(values = benchmark_cols ,
                      breaks = c("rf_simple","rf_benchmark","EL_bechmark"), 
                      labels = c("IMPROVE simple CV","IMPROVE simple bencmark","RankEL benchmark")) +
  theme_bw() +
  labs(color  = "Model", y="Sensitivity (TPR)", x="Specificity")+
  theme(plot.title = element_text(size = 14),legend.position = "right",
        legend.text = element_text(size =14 ),
        axis.title = element_text(size = 18 ),
        legend.title = element_text(size = 16 )) + 
  guides(colour = guide_legend(override.aes = list(size=8)))
#ggsave(Fig6A, file = "results/PaperPlots/Fig6/Fig6A.pdf", width = 8 , height = 6 )

roc.test(rocobj_benchmark_el ,rocobj_benchmark)

simple_dat <- get_pr_dat(model = "IMPROVE_simple", predictions = pred_df_simple$prediction_rf, labels  = pred_df_simple$response)
benchmark_dat <- get_pr_dat(model = "IMPROVE_Becnhmark", predictions = benchmark_data_prediction$prediction_rf, labels  = benchmark_data_prediction$response)
EL_benchmark_dat <- get_pr_dat(model = "RankEL_Benchmark", predictions = -benchmark_data_prediction$RankEL, labels  = benchmark_data_prediction$response)


auc(improve_dat$recall, improve_dat$precision)

pr_all_dat_benchmark <- bind_rows(simple_dat,benchmark_dat) %>% bind_rows(.,EL_benchmark_dat ) 
pr_all_dat_benchmark  <- na.omit(pr_all_dat_benchmark )

pr_curve_benchmark_fig_6  <- 
  pr_all_dat_benchmark %>% 
  ggplot(. , aes(x =recall, y = percision )) +
  geom_point(aes(color = model))  +
  geom_line(aes(color = model)) +
  scale_x_log10() + 
 # scale_y_continuous(limits = c(0,0.3)) + 
  theme_bw()
ggsave(pr_curve_benchmark_fig_6, file = "results/PaperPlots/Fig6/pr_curve_benchmark_fig_6.pdf", height = 5, width = 7 )

PRAUC(pred_df_simple$prediction_rf, pred_df_simple$response)
PRAUC(benchmark_data_prediction$prediction_rf, benchmark_data_prediction$response)
PRAUC(-benchmark_data_prediction$RankEL, benchmark_data_prediction$response)



### COMPARE TO OTHER TOOLS 


# result simple_model 
# ------------------------------
benchmark_compareison <- read.table(file = "data/benchmark_comparison/neoepitopes_predictions.tsv", sep = "\t", header = T)
AUC_comparison <- read.table(file = "data/benchmark_comparison/neoepitopes_predictions_performance.tsv", sep = "\t", header = T)
# -------------
length(na.omit(benchmark_compareison$HLAthena.MSiCE))

auc_df <- data.frame()
for (col in colnames(benchmark_compareison)[5:length(colnames(benchmark_compareison))]) {
  print(col)
  if (col != "response") {
    df = benchmark_compareison[,c(col,"response")]
    df =  na.omit(df)
    y = df["response"] %>% pull()
    x = df[col]  %>% pull()
    rocobj <- roc(y,x, direction = "<")
    #define object to plot and calculate AUC
    AUC <- round(auc(y, x),3)
    # auc 0.1
    pred <- prediction(x,y)
    auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
    auc_01 <- round(auc_01@y.values[[1]],4)
    auc_df <- rbind(auc_df,c(col,AUC,auc_01))
  }
}
colnames(auc_df) <- c("model","AUC","AUC01")

auc_df <- rbind(auc_df,c("IMPROVE TME", 0.649,0.0146))
auc_df <- rbind(auc_df,c("IMPROVE", 0.631,0.0142))
auc_df <- rbind(auc_df,c("IMPROVE Simple", 0.643,0.0134))


auc_df$AUC <- as.numeric(auc_df$AUC)
auc_df$AUC01 <- as.numeric(auc_df$AUC01)


Fig6B <- auc_df %>%
  filter(!model %in% c("PRIME.score","PRIME.rank")) %>% 
  filter(model!="id") %>% 
  ggplot(., aes(x = reorder(model, -AUC01), y = AUC01)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 12))+
  labs(x = "", y = "AUC01")

Fig6C <- auc_df %>% 
  filter(!model %in% c("PRIME.score","PRIME.rank")) %>% 
  filter(model!="id") %>% 
  ggplot(., aes(x = reorder(model, -AUC), y = AUC)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 12))+
  labs(x = "", y = "AUC") +
#ggsave(Fig6B , file = "results/PaperPlots/Fig6/Fig6B.pdf", width = 12 , height = 6 )

pdf(file = 'results/PaperPlots/Fig6/Figure6_ALL.pdf', width = 8, height = 16)
ggdraw() +
  draw_plot(Fig6A, .0, .8, 1, 0.2) +
  draw_plot(Fig6B, .0, .4, 1, 0.4) +
  draw_plot(Fig6C, .0, .0, 1, 0.4) +
  draw_plot_label(c("A","B","C"), x = c(0,0,0), y = c(1.005,0.8,0.4), size = 20)  

dev.off()
# 
# 

### compare old run 


benchmark_data_prediction_old <- read.table("/Users/annieborch/Documents/GitHub/Immunugenicity/results/RandomForrest/benchmark_data/benchmark_data_prediction_test.txt", sep = "\t", header = T)
benchmark_data_prediction_old <- benchmark_data_prediction_old %>% mutate(pMHC = paste(HLA_allele,Peptide, sep = "_"))

benchmark_data_prediction[!benchmark_data_prediction$Expression %in% benchmark_data_prediction_old$Expression_Level,]



table(benchmark_data_prediction$Expression %in% benchmark_data_prediction_old$Expression_Level)

