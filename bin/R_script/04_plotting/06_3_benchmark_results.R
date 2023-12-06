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


### COMPARE TO OTHER TOOLS with In-house data 


# result simple_model 
# ------------------------------
benchmark_compareison <- read.table(file = "data/benchmark_comparison/neoepitopes_predictions.tsv", sep = "\t", header = T)
AUC_comparison <- read.table(file = "data/benchmark_comparison/neoepitopes_predictions_performance.tsv", sep = "\t", header = T)
# -------------
length(na.omit(benchmark_compareison$HLAthena.MSiCE))
nrow(benchmark_compareison)

Pred_Modelling_to_join <- Pred_Modelling %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide,sep = "_"))%>% 
  select(pMHC, Prime, Foreigness, RankEL_minus) 


benchmark_compareison <- benchmark_compareison %>%  
  mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_"), .after = "Mut_peptide") %>% 
  distinct(pMHC , .keep_all = T)

table(benchmark_compareison$pMHC %in% Pred_Modelling_to_join$pMHC)
benchmark_compareison_test <- Pred_Modelling_to_join  %>%
  left_join(.,benchmark_compareison, by = "pMHC")  %>% 
  select(-c(HLA_allele,Norm_peptide,Mut_peptide))

table(benchmark_data_prediction_filter$response)

auc_df <- data.frame()
for (col in colnames(benchmark_compareison_test)[2:length(colnames(benchmark_compareison_test))]) {
  print(col)
  if (col != "response") {
    df = benchmark_compareison_test[,c(col,"response")]
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
#auc_df <- rbind(auc_df,c("RankEL", 0.583,0.0109))


auc_df$AUC <- as.numeric(auc_df$AUC)
auc_df$AUC01 <- as.numeric(auc_df$AUC01)
auc_df$model

auc_df <- auc_df %>% 
  mutate(Model_name = 
           case_when(model == "DeepNetBim.immunogenicity.probability" ~ "NetDeppBim immunogenicty probalility",
                   #  model == "IMPROVE Simple" ~ "IMPROVE simple",
                     model == "IMPROVE TME" ~ "IMPROVE TME",
                     model == "IMPROVE" ~ "IMPROVE",
                     model == "Foreigness" ~ "Foreigness",
                     model == "PRIME.score" ~ "Prime",
                     model == "iTTCA.RF.probability" ~ "iTTCA RF probability",
                     model == "IEDB.immunogenicity" ~ "IEDB immunogenicity",
                     model == " MHCflurry.processing " ~ "MHCflurry prediction",
                     model == "MHCflurry.presentation" ~ "MHCflurry prediction percentile",
                     model == "MixMHCpred" ~ "MixMHCpred score",
                     #   model == "MixMHCpred" ~ "MixMHCpred rank",
                     model == "RankEL_minus" ~ "RankEL",
                     TRUE ~ "other"
           ))




Fig6A <- auc_df %>%
  filter(Model_name!="other") %>% 
  filter(model!="id") %>% 
  ggplot(., aes(x = reorder(Model_name, -AUC01), y = AUC01)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  labs(x="",y="AUC 10% score")+
  geom_hline(yintercept = 0.005, linetype = "dashed") + 
  theme(axis.text.x = element_text(angle = 90, size = 12,hjust = 0.9,vjust=0.5))+
  labs(x = "", y = "AUC01")

Fig6B <- auc_df %>% 
  filter(Model_name!="other") %>% 
  filter(model!="id") %>% 
  ggplot(., aes(x = reorder(Model_name, -AUC), y = AUC)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  coord_cartesian(ylim=c(0.4, 0.7))+
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  # scale_y_continuous(limits = c(0.4,1)) +
  theme(axis.text.x = element_text(angle = 90, size = 12,hjust = 0.9,vjust=0.5))+
  labs(x = "", y = "AUC") 
#ggsave(Fig6B , file = "results/PaperPlots/Fig6/Fig6B.pdf", width = 12



# result simple_model 
# ------------------------------
benchmark_data_prediction <- read.table(file = "results/benchmark_data/04_Benchmark_neoepitopes_prediction.tsv", sep = "\t", header = T)

benchmark_data_prediction <- benchmark_data_prediction %>% 
  rename("prediction_rf" = "mean_prediction_rf") %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide , sep = "_")) %>% 
distinct(pMHC, .keep_all = T) %>% 
  mutate(RankEL_minus = -RankEL)  

colnames(benchmark_data_prediction)

table(benchmark_data_prediction$response)
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

Fig6C <- ggroc(list(rf_simple = rocobj_simple, 
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
  theme(plot.title = element_text(size = 14),legend.position = "bottom",
        legend.text = element_text(size =14 ),
        axis.title = element_text(size = 18 ),
        legend.title = element_text(size = 16 )) + 
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(color=guide_legend(ncol=1, title.position = "top"))
#ggsave(Fig6A, file = "results/PaperPlots/Fig6/Fig6A.pdf", width = 8 , height = 6 )

roc.test(rocobj_benchmark_el ,rocobj_benchmark)

### COMPARE TO OTHER TOOLS with CEDAR data
# -----------------------------------------
load("results/benchmark_data/plotting_data/CEDAR_data_peformance.Rdata")
auc_df_bench_V2$model
auc_df_bench_V2 <- auc_df_bench_V2 %>% 
  mutate(Model_name = 
case_when(model == "immuno_probability_netdeepbim" ~ "NetDeppBim immunogenicty probalility",
         model == "IMPROVE" ~ "IMPROVE simple",
         model == "Foreigness" ~ "Foreigness",
         model == "Prime_new" ~ "Prime (without training pMHC)",
         model == "Prime" ~ "Prime (All data)",
         model == "iTTCA_score" ~ "iTTCA RF probability",
         model == "IEDB_score" ~ "IEDB immunogenicity",
         model == "mhcflurry_prediction" ~ "MHCflurry prediction",
         model == "mhcflurry_prediction_percentile" ~ "MHCflurry prediction percentile",
         model == "mixMHCpred_Score" ~ "MixMHCpred score",
       #  model == "mixMHCpred_Rank" ~ "MixMHCpred rank",
         model == "RankEL" ~ "RankEL",
         TRUE ~ "other"
         )) %>% 
  mutate(Model_group = 
  case_when(model == "immuno_probability_netdeepbim" ~ "deepbim",
            model == "IMPROVE" ~ "IMPROVE",
            model == "Foreigness" ~ "Foreigness",
            model == "Prime_new" ~ "Prime",
            model == "Prime" ~ "Prime",
            model == "iTTCA_score" ~ "iTTCA ",
            model == "IEDB_score" ~ "IEDB",
            model == "mhcflurry_prediction" ~ "MHCflurry",
            model == "mhcflurry_prediction_percentile" ~ "MHCflurry",
            model == "mixMHCpred_Score" ~ "MixMHCpred",
            model == "mixMHCpred_Rank" ~ "MixMHCpred",
            model == "RankEL" ~ "RankEL",
            TRUE ~ "other"
  ))

#colors_to_benchmark <- c("#3288bd",)

Fig6D <- auc_df_bench_V2 %>%
  filter(Model_name!="other") %>% 
  #filter(!model %in% c("PRIME.score","PRIME.rank")) %>% 
  #  filter(model!="id") %>% 
  ggplot(., aes(x = reorder(Model_name, -AUC01), y = AUC01)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  labs(x="",y="AUC 10% score")+
  geom_hline(yintercept = 0.005, linetype = "dashed") + 
  theme(axis.text.x = element_text(angle = 90, size = 12,hjust = 0.9,vjust=0.5))+
  labs(x = "", y = "AUC01")

Fig6E <- auc_df_bench_V2 %>% 
  filter(Model_name!="other") %>% 
  # filter(!model %in% c("PRIME.score","PRIME.rank")) %>% 
  ggplot(., aes(x = reorder(Model_name, -AUC), y = AUC)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  coord_cartesian(ylim=c(0.4, 0.7))+
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  # scale_y_continuous(limits = c(0.4,1)) +
  theme(axis.text.x = element_text(angle = 90, size = 12,hjust = 0.9,vjust=0.5))+
  labs(x = "", y = "AUC") 
#ggsave(Fig6B , file = "results/PaperPlots/Fig6/Fig6B.pdf", width = 12 , height = 6 )







# --------------------------------------------------
### simple model without prime 
# --------------------------------------------------
simple_wo_prime <- read.table("results/5_fold_CV/Simple/pred_df_wo_primeSimple.txt", header = T)

CEDAR_prediction_wo_prime <- read.table("results/benchmark_data/04_Benchmark_neoepitopes_prediction_wo_prime.tsv", header = T)


# In-house data 
rocobj_simple_wo_prime <- roc(simple_wo_prime$response, simple_wo_prime$prediction_rf, direction = "<")
#define object to plot and calculate AUC
auc_simple_wo_prime <- round(auc(simple_wo_prime$response, simple_wo_prime$prediction_rf),3)
# auc 0.1
pred <- prediction(simple_wo_prime$prediction_rf,simple_wo_prime$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_simple_wo_prime_0.1 <- round(auc_01@y.values[[1]],4)

# cedar data 
rocobj_simple_wo_prime_CEDAR <- roc(CEDAR_prediction_wo_prime$response, CEDAR_prediction_wo_prime$mean_prediction_rf, direction = "<")
#define object to plot and calculate AUC
auc_simple_wo_prime_CEDAR <- round(auc(CEDAR_prediction_wo_prime$response, CEDAR_prediction_wo_prime$mean_prediction_rf),3)
# auc 0.1
pred <- prediction(CEDAR_prediction_wo_prime$mean_prediction_rf,CEDAR_prediction_wo_prime$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_simple_wo_prime_CEDAR_0.1 <- round(auc_01@y.values[[1]],4)



Fig6F <- ggroc(list(rf_simple_wo_prime = rocobj_simple_wo_prime, 
                    rf_simple_wo_prime_benchmark = rocobj_simple_wo_prime_CEDAR)) +
  # ggtitle(paste0('(AUC IMPROVE simple CV = ', auc_simple, ')',
  #                '(AUC IMPROVE simple bencmark = ', auc_benchmark , ')\n',
  #                ' (AUC RankEL benchmark = ', auc_benchmark_el , ')\n',
  #                '(AUC01 IMPROVE simple = ', auc_simple_0.1, ')',
  #                '(AUC01 IMPROVE bencmark = ', auc_benchmark_0.1, ')\n',
  #                ' (AUC01 RankEL bencmark = ', auc_benchmark_el_0.1,')')) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0),
               color="darkgrey", linetype="dashed") +
  scale_color_manual(values = benchmark_cols ,
                     breaks = c("rf_simple_wo_prime","rf_simple_wo_prime_benchmark"), 
                     labels = c("IMPROVE simple without Prime \n 5-fold CV","IMPROVE simple withot Prime \n bencmark prediction")) +
  theme_bw() +
  labs(color  = "Model", y="Sensitivity (TPR)", x="Specificity")+
  theme(plot.title = element_text(size = 14),legend.position = "bottom",
        legend.text = element_text(size =14 ),
        axis.title = element_text(size = 18 ),
        legend.title = element_text(size = 16 )) + 
  guides(color=guide_legend(ncol=1, title.position = "top"))

  
  
  
  
  

pdf(file = 'results/PaperPlots/Fig6/Figure6.pdf', width = 12, height = 11)
ggdraw() +
  draw_plot(Fig6A, .0, .5, 0.33, 0.5) +
  draw_plot(Fig6B, .33, .5, 0.33, 0.5) +
  draw_plot(Fig6C, .66, .61, 0.33, 0.39) +
  draw_plot(Fig6D, .0, .0, 0.33, 0.5) +
  draw_plot(Fig6E, .33, .0, 0.33, 0.5) +
  draw_plot(Fig6F, .66, .11, 0.33, 0.39) +
draw_plot_label(c("A","B","C"), x = c(-0.002,.325,.655), y = c(1.005,1.005,1.005), size = 20)  +
draw_plot_label(c("D","E","D"), x = c(-0.002,.325,.655), y = c(0.515,0.515,0.515), size = 20)  
dev.off()






