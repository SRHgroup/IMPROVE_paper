 

load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")


library(GGally)
library(corrplot)
library(openxlsx)

# fore roc curves 
# Loading package
library(caTools)
library(randomForest)
library(pROC)
library(ROCit)
library(caret)
library(ROCR)
library(verification)

# survival packages 
library(survminer)
library(survival)

library(rstatix)

library(data.table)

library(MLmetrics)

unique(all_peptides$Sample)
# --------------------------------------------
#               Figure 3 
# --------------------------------
Pred_Modelling <- Pred_Modelling %>% mutate(assembel_score_rf = (Prediction+prediction_rf)/2)
class(Pred_Modelling$Prediction)
## NNalign 
# ---------------------------------------------
Pred_Modelling$Prediction <- as.numeric(Pred_Modelling$Prediction)
rocobj_nnalign <- roc(Pred_Modelling$response,Pred_Modelling$Prediction, direction = "<")
auc_nnalign <- round(auc(Pred_Modelling$response,Pred_Modelling$Prediction),3)
# Auc 01 
pred <- prediction(Pred_Modelling$Prediction,Pred_Modelling$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_nnalign_0.1 <- round(auc_01@y.values[[1]],4)

## random forrest 
rocobj_rf <- roc(Pred_Modelling$response, Pred_Modelling$prediction_rf, direction = "<")
auc_rf <- round(auc(Pred_Modelling$response, Pred_Modelling$prediction_rf),3)
# Auc 01 
pred <- prediction(Pred_Modelling$prediction_rf,Pred_Modelling$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_rf_0.1 <- round(auc_01@y.values[[1]],4)

## random forrest assembling 
rocobj_rf_assembel <- roc(Pred_Modelling$response, Pred_Modelling$assembel_score_rf, direction = "<")
auc_rf_assembel <- round(auc(Pred_Modelling$response, Pred_Modelling$assembel_score_rf),3)
# Auc 01 
pred_assembel <- prediction(Pred_Modelling$assembel_score_rf,Pred_Modelling$response)
auc_01_assembel <- performance(pred_assembel,measure = "auc", fpr.stop=0.1)
auc_rf_0.1_assembel <- round(auc_01_assembel@y.values[[1]],4)


## NNalign
rocobj_RankEL <- roc(all_peptides$response, all_peptides$RankEL_minus, direction = "<")
auc_RankEL <- round(auc(all_peptides$response, all_peptides$RankEL_minus),3)
# Auc 01 
pred_RankEL <- prediction(all_peptides$RankEL_minus,all_peptides$response)
auc_01_RankEL <- performance(pred_RankEL,measure = "auc", fpr.stop=0.1)
auc_0.1_RankEL <- round(auc_01_RankEL@y.values[[1]],4)



Pred_Modelling <- Pred_Modelling %>% 
  mutate(pMHC = paste(HLA_allele,Mut_peptide, sep = "_"))



##Fig3B plot roc curve 
# -------------------------------------------------------------
Fig3A <- ggroc(list(nnalign = rocobj_nnalign,RankEL =rocobj_RankEL ,rf = rocobj_rf ,assembel = rocobj_rf_assembel )) +
  # ggtitle(paste0('(AUC IMPROVE = ', auc_rf, ')','(AUC NNalign = ', auc_nnalign, ')\n',
  #                '(AUC RankEL = ', auc_RankEL, ')','(AUC Ensemble = ', auc_rf_assembel, ')\n',
  #                '(AUC01 IMPROVE = ', auc_rf_0.1, ') ',
  #                ' (AUC01 NNalign = ', auc_nnalign_0.1, ')\n',
  #                '(AUC01 RankEL = ', auc_0.1_RankEL, ') ','(AUC01 Ensemble = ', auc_rf_0.1_assembel, ')')) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0),
               color="darkgrey", linetype="dashed") +
  scale_color_manual(values = model_col_rf_TME_NNalign_assembling_RankEL , 
                     breaks = c("rf","nnalign","assembel","RankEL"),
                     labels = c("IMPROVE","NNalign","Ensemble","RankEL")) +
  theme_bw() +
  # scale_x_continuous(breaks = c(1,0.75,0.5,0.25,0), 
  #                   labels = c("0","0.25","0.5","0.75","1")) +
  labs(color  = "Model", y="Sensitivity (TPR)", x="Specificity")+
  theme(plot.title = element_text(size = 12),
        axis.title = element_text(size = 16 ),
        legend.position = "bottom",
        legend.text = element_text(size =14 ),
        legend.title = element_text(size = 16 ))+ 
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2))

#ggsave(Fig3A, file = "results/PaperPlots/Fig3/Fig3A.pdf", width = 5 , height = 4 )




roc.test(rocobj_nnalign , rocobj_rf )
roc.test(rocobj_RankEL, rocobj_rf  )

#are.paired(rocobj_nnalign , rocobj_rf)
# figure 3b feature importance 
colnames(Feature_importance) <- c("Feature","importance","Partition")

table(Feature_importance$Feature)
nrow(Feature_importance)
# bar plot
#feature_sum <- as.data.frame(aggregate(Feature_importance$importance, by=list(feature = Feature_importance$Feature), FUN=sum))
#feature_sum <- feature_sum %>% mutate(mean_imp = x/250)

feature_sum  <- Feature_importance %>% group_by(Partition,Feature) %>% summarize(mean_imp_part = mean(importance, na.rm=TRUE))

feature_sum <- feature_sum %>% mutate(mean_imp =  mean_imp_part/5)

feature_sum <- feature_sum %>% 
  mutate(Feature_category = 
           case_when(Feature %in% c("Aro","CysRed","PropHydroAro","Inst","HydroCore",
                                    "mw","pI","pI","PropAcidic","PropAro","PropAromatic","PropBasic",
                                    "PropSmall","Prime","HydroAll") ~ "Physicochemical properties",
                     Feature %in% c("CelPrev","Expression","VarAlFreq") ~ "Mutation qualities",
                     Feature %in% c("RankEL","RankBA","Stability","Stability","NetMHCExp") ~ "Peptide-MHC",
                     Feature %in% c("SelfSim","Foreigness","DAI") ~ "Comparing normal",
                     Feature=="PrioScore" ~ "Other"))



Fig3B <-  feature_sum %>% 
  ggplot(., aes(x = reorder(Feature, -mean_imp), y = mean_imp)) +
  geom_bar(aes(fill =Feature_category), stat = "identity")+
  theme(axis.text.x = element_text(angle = 90))  +
  theme_bw() +
  scale_fill_manual(breaks = c("Physicochemical properties","Peptide-MHC","Mutation qualities","Comparing normal","Other"),
                    values = c("#abdda4","#d7191c","#fdae61","#2b83ba","#ffffbf"))  +
  scale_x_discrete(labels=Abbreviations) +
  labs(x = "", y = "Mean importance", fill = "Feature category") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size=14), 
        legend.position = "bottom", 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14), 
        axis.title = element_text(size = 16)) +
  guides(fill=guide_legend(ncol=2, title.position = "top"))
#ggsave(Fig3B, file = "results/PaperPlots/Fig3/Fig3B.png", width = 8 , height = 6 )



Fig3C <- Pred_Modelling %>% 
  # filter(cohort == "bladder" ) %>% 
  ggplot(., aes( x = response_lab , y = prediction_rf))+
  geom_boxplot(aes(fill = response_lab),alpha = 0.5)+
  theme_bw() +
  scale_fill_manual(breaks = c("no","yes"), values = response_col) +
  scale_y_continuous(limits = c(0.2,0.91))+
  geom_boxplot(aes(fill = response_lab),alpha = 0.5)+ 
  labs(x = "", y = "IMPROVE \n prediction score", fill = "Immunogenic") +
  geom_signif(comparisons = list(c("no","yes")),
              na.rm = T,
              textsize = 7,
              map_signif_level = T,
              data = Pred_Modelling ,
              test = wilcox.test)  +
  facet_grid(.~cohort) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        strip.text.x = element_text(size = 18))

Fig3C_2 <- Pred_Modelling %>% 
  # filter(cohort == "bladder" ) %>% 
  ggplot(., aes( x = response_lab , y = Prediction))+
  geom_boxplot(aes(fill = response_lab),alpha = 0.5)+
  scale_y_continuous(limits = c(0,0.22))+
  theme_bw() +
  scale_fill_manual(breaks = c("no","yes"), values = response_col) +
  # scale_y_continuous(breaks = c(0,0.02,0.04,0.06,0.08,0.1,0.125,0.15,0.175,0.2,0.225))+
  geom_boxplot(aes(fill = response_lab),alpha = 0.5)+ 
  labs(x = "", y = "NNalign \n prediction score", fill = "Immunogenic") +
  geom_signif(comparisons = list(c("no","yes")),
              na.rm = T,
              map_signif_level = T,
              textsize = 7,
              data = Pred_Modelling ,
              test = wilcox.test)  +
  facet_grid(.~cohort) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16), 
        strip.text.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 18), legend.title =element_text(size = 20) ) + 
  guides(fill = guide_legend(override.aes = list(size=1.5)))

Fig3C_legend <- get_legend(Fig3C_2)
Fig3C_2 <- Fig3C_2 + theme(legend.position = "none")


#              Sup Figure 3 
# --------------------------------
cols_selected <- c("mw" , "Aro" , "Inst"  ,"CysRed", "HydroAll",
                   "RankEL","RankBA", "Expression",
                   "PrioScore","VarAlFreq" ,"CelPrev" ,"SelfSim", 
                   "Prime", "PropHydroAro","HydroCore" ,             
                   "PropAro","PropBasic","PropAcidic" ,  
                   "pI","DAI",  "Stability","Foreigness","NetMHCExp")
Abbreviations_cor <- c("mw","Aro","Inst","CysRed","HydroAll","RankEL","RankBA","Expression","PrioScore","VarAlFrac",
                       "CelPrev","SelfSim","Prime","PropHydroAro","HydroCore","PropBasic","PropAcidic",
                       "PropAromatic",'PI',"DAI","Stability","Foreigness","NetMHCExp")
Correlation_all_data  <- all_peptides %>% 
  dplyr::select(cols_selected) 
Correlation_all_data[cols_selected] <- sapply(Correlation_all_data[cols_selected],as.numeric)

#colnames(Correlation_all_data) <- Abbreviations_cor
SupFig3A <- ggcorr(Correlation_all_data , method = c("pairwise", "spearman"), 
                   hjust = 0.8, size = 6, label = TRUE,label_size = 3) +
  labs(fill = "spearman", color ="spearman" )
#ggsave(SupFig3A, file = "results/PaperPlots/Fig3/SupFig3A.pdf", width = 10, height = 10)


# subFigB percision recall curves 
improve_dat <- get_pr_dat(model = "improve", predictions = Pred_Modelling$prediction_rf, labels  = Pred_Modelling$response)
NNalign_dat <- get_pr_dat(model = "NNalign", predictions = Pred_Modelling$Prediction, labels  = Pred_Modelling$response)
rankel_dat <- get_pr_dat(model = "Rankel", predictions = Pred_Modelling$RankEL_minus, labels  = Pred_Modelling$response)

pr_all_dat <- bind_rows(improve_dat,NNalign_dat) %>% bind_rows(.,rankel_dat )
pr_all_dat <- na.omit(pr_all_dat)

SupFig3B  <- 
  pr_all_dat %>%  
  ggplot(. , aes(x =recall, y = percision )) +
  geom_point(aes(color = model))  +
  geom_line(aes(color = model)) +
  scale_color_manual(values = model_col_rf_NNalign_RankEL ,
                     breaks = c("improve","NNalign","Rankel"),
                     labels = c("IMPROVE","NNalign","RankEL")) +
  scale_x_log10() + 
  scale_y_continuous(limits = c(0,0.3)) + 
  theme_bw()+
  theme(legend.position = "none")

PRAUC(Pred_Modelling$Prediction, Pred_Modelling$response)
PRAUC(Pred_Modelling$RankEL_minus, Pred_Modelling$response)
PRAUC(Pred_Modelling$prediction_rf, Pred_Modelling$response)
#ggsave(pr_curve_all_fig_3, file = "results/PaperPlots/Fig3/pr_curveall_fig_3.png")



# Subfig3C nn align logoplot 

# SupFig3 B,C and D
SupFig3F <- pred_per_patient(y_val = 'prediction_rf', c = "mUC" )
SupFig3E <-  pred_per_patient(y_val = 'prediction_rf',c = "Melanoma", legpos ="bottom",h = 5)
SupFig3Elegend <- get_legend(SupFig3E )
SupFig3E <- SupFig3E + theme(legend.position = "none")
SupFig3D <- pred_per_patient(y_val = 'prediction_rf',c = "Basket")

pred_per_patient(y_val = 'prediction_rf_tme', c = "mUC")
pred_per_patient(y_val = 'prediction_rf_tme',c = "Melanoma")
pred_per_patient(y_val = 'prediction_rf_tme',c = "Basket")


Pred_Modelling  <- Pred_Modelling %>% 
  mutate(pMHC = paste(Mut_peptide,HLA_allele, sep = "_")) %>% 
  group_by(pMHC) %>% 
  add_tally(name = "dup")  %>% 
  mutate(duplicates_group = case_when(dup>1 ~ "Shared",
                                      TRUE ~ "Unique")) 



Pred_Modelling  <- Pred_Modelling %>% 
  group_by(Mut_peptide,Patient) %>% 
  add_tally(name = "dup_pep_patient")  %>% 
  group_by(Mut_peptide) %>% 
  add_tally(name = "peptide_dup")  %>%
  mutate(shared_peptides =
           case_when(peptide_dup>dup_pep_patient ~ "shared",
                     TRUE ~ "not-shared"))

all_peptides  <- all_peptides %>% 
  group_by(pMHC) %>% 
  add_tally(name = "dup")  %>% 
  mutate(duplicates_group = case_when(dup>1 ~ "Shared",
                                      TRUE ~ "Unique")) 
table(all_peptides$response)

shared_pep <- all_peptides %>% filter(duplicates_group=="Shared")
shared_pep_response <- all_peptides %>% filter(duplicates_group=="Shared", response==1)

table(shared_pep$response)
nrow(shared_pep)
shared_pep_distinct <- shared_pep %>% distinct(pMHC, .keep_all=T)
Unique_pep <- all_peptides %>% filter(duplicates_group=="Unique")

# gene_with_response <- all_peptides %>% filter(response==1) %>% ungroup()  %>% 
#   distinct(Gene_Symbol, .keep_all=T) %>% select(Gene_Symbol)
# write.table(gene_with_response, file="tabels/gene_with_response.txt", quote = F, row.names = F, col.names = F)
# 


shared_pep <- shared_pep %>%
  mutate(dup_gene = paste(pMHC,Gene_Symbol)) %>% 
  group_by(dup_gene) %>% 
  add_tally(name = "dup_gene_number") %>% 
  mutate(not_dup_gene = case_when(dup_gene_number!=dup ~ "no_dup_gene",
                                  TRUE ~"dup_gene"))


no_dup_gene <- shared_pep %>% filter(not_dup_gene=="no_dup_gene")


Fig_3F_RankEL <- Pred_Modelling %>% filter(response==0) %>% 
  mutate(Immu_lab = "non-immunogenic") %>% 
  ggplot(. , aes(x = duplicates_group, y = RankEL)) +
  geom_quasirandom(aes(color = duplicates_group),size = 3) + 
  geom_boxplot(aes(fill = duplicates_group),alpha = 0.5)+
  # geom_line(aes(group = pMHC)) +
  scale_color_manual(breaks = c("Shared","Unique"), values = c("#998ec3","#f1a340")) +
  scale_fill_manual(breaks = c("Shared","Unique"), values = c("#998ec3","#f1a340") ) +
  # geom_text(aes(label = pMHC)) +
  theme_bw()+
  scale_y_log10(limits = c(0.001, 48), breaks = c(0.001,0.005,0.1,1,5,15),
                labels = c("0.001","0.005","0.1","1","5","15")) + 
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "pMHC", y = "RankEL prediction score") +
  facet_grid(. ~ Immu_lab) +
  geom_signif(comparisons = list(c("Shared","Unique")),
              na.rm = T,
              map_signif_level = T,
              textsize = 7,
              data = Pred_Modelling,
              test = wilcox.test)  
#ggsave(Fig_3F_RankEL, file = "results/PaperPlots/Fig3/Fig_3F_RankEL.pdf", width = 3, height = 3)


Fig_3F_prediction_rf <- Pred_Modelling %>% filter(response==0) %>% 
  mutate(Immu_lab = "non-immunogenic") %>% 
  ggplot(. , aes(x = duplicates_group, y = prediction_rf)) +
  geom_quasirandom(aes(color = duplicates_group),size = 3) + 
  geom_boxplot(aes(fill = duplicates_group),alpha = 0.5)+
  # geom_line(aes(group = pMHC)) +
  scale_y_continuous(limits = c(0,0.9))+
  scale_color_manual(breaks = c("Shared","Unique"), values = c("#998ec3","#f1a340") ) +
  scale_fill_manual(breaks = c("Shared","Unique"), values = c("#998ec3","#f1a340") ) +
  # geom_text(aes(label = pMHC)) +
  theme_bw()+
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "pMHC",y = "IMPROVE prediction score") +
 # scale_y_continuous(limits = c(0.2,0.9))+
  facet_grid(. ~ Immu_lab) +
  geom_signif(comparisons = list(c("Shared","Unique")),
              na.rm = T,
              textsize = 7,
              map_signif_level = F,
              data = Pred_Modelling,
              test = wilcox.test)  

#ggsave(Fig_3F_prediction_rf, file = "results/PaperPlots/Fig3/Fig_3F_prediction_rf.pdf", width = 3, height = 3)






##### model per patein 
p_vec <- Pred_Modelling %>% filter(response==1) %>% pull(Patient) 
AUC_score_pp <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(AUC_score_pp) <- c("Patient","AUC","AUC01","type")
Pred_Modelling$Prediction <- as.numeric(Pred_Modelling$Prediction)
Pred_Modelling$prediction_rf <- as.numeric(Pred_Modelling$prediction_rf)
Pred_Modelling$prediction_rf_tme <- as.numeric(Pred_Modelling$prediction_rf_tme)
Pred_Modelling$response<- as.numeric(Pred_Modelling$response, levels = c(1,0))


Auc_pp(name = "NNalign" , var = "Prediction")
Auc_pp(name = "rf" , var = "prediction_rf")
Auc_pp(name = "rf_tme" , var = "prediction_rf_tme")
Auc_pp(name = "RankEL" , var = "RankEL_minus")

colnames(AUC_score_pp) <- c("Patient","AUC","AUC01","type")
AUC_score_pp[c("AUC","AUC01")] <- sapply(AUC_score_pp[c("AUC","AUC01")],as.numeric)



AUC_score_pp_long <- AUC_score_pp %>% 
  gather(., key = "cal_type", value = "value", -c(type,Patient))

AUC_score_pp_long$value <- as.numeric(AUC_score_pp_long$value)
all_peptides$Patient
cohort <- all_peptides %>% ungroup() %>%  select(cohort,Patient) %>% distinct(Patient, .keep_all=T)
AUC_score_pp_long <- AUC_score_pp_long %>% left_join(.,cohort)

num_response <- Pred_Modelling %>% filter(response==1) %>% group_by(Patient) %>% tally()

AUC_score_pp_long$type <- factor(AUC_score_pp_long$type, levels = c("RankEL","NNalign","rf","rf_tme"))




SupFig3C <- AUC_score_pp_long %>%
  filter(cal_type=="AUC") %>% filter(type!="rf_tme") %>% 
  ggplot(., aes(x = type, y= value, group =type )) +
  geom_boxplot(aes(x = type, y= value, group = type , fill = type), alpha = 0.8) +
  geom_beeswarm(aes(x = type, y= value, group =type , color = type))+
  geom_line(aes(group = Patient),alpha = 0.1) + 
  facet_grid(.~cal_type)+
  theme_bw() +
  scale_x_discrete(breaks = c("rf","NNalign","RankEL"),
                   labels = c("IMPROVE","NNalign","RankEL")) +
  scale_y_continuous(limits = c(0.3,1.15)) +
  scale_color_manual(values = model_col_rf_NNalign_RankEL ,
                     breaks = c("rf","NNalign","RankEL"),
                     labels = c("IMPROVE","NNalign","RankEL")) +
  scale_fill_manual(values = model_col_rf_NNalign_RankEL ,
                    breaks = c("rf","NNalign","RankEL"),
                    labels = c("IMPROVE","NNalign","RankEL")) +
  labs(x="",y="AUC score", fill = "model",color = "model")+
  geom_signif(comparisons = list(c("NNalign","rf"),
                                 c("RankEL","rf")),
              na.rm = T,
              map_signif_level = F,
              data = AUC_score_pp_long %>% filter(cal_type=="AUC")%>% filter(type!="rf_tme"),
              test = wilcox.test,
              textsize = 5,
              step_increase = 0.15,
              test.args = list(paired = T, exact = F))  +
  theme(legend.position = "none",axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), strip.text.x =element_text(size = 14) )
#ggsave(Fig3new_auc, file="results/PaperPlots/Fig3/Fig3new_auc.pdf", height = 4, width = 5 )


Fig3D  <- AUC_score_pp_long %>% filter(cal_type=="AUC01") %>% filter(type!="rf_tme") %>% 
  #filter(type %in% c("NNalign","rf","rf_assemling")) %>% 
  ggplot(., aes(x = type, y= value, group =type )) +
  geom_boxplot(aes(x = type, y= value, group = type , fill = type), alpha = 0.8) +
  geom_beeswarm(aes(x = type, y= value, group =type , color = type))+
  geom_line(aes(group = Patient),alpha = 0.1) + 
  scale_y_continuous(limits = c(0,0.1)) +
  facet_grid(.~cal_type)+
  theme_bw() +
  scale_x_discrete(breaks = c("rf","NNalign","RankEL"),
                   labels = c("RF","NNalign","RankEL")) +
  scale_color_manual(values = model_col_rf_NNalign_RankEL ,
                     breaks = c("rf","NNalign","RankEL"),
                     labels = c("IMPROVE","NNalign","RankEL")) +
  scale_fill_manual(values = model_col_rf_NNalign_RankEL ,
                    breaks = c("rf","NNalign","RankEL"),
                    labels = c("IMPROVE","NNalign","RankEL")) +
  labs(x="",y="AUC 10% score", fill = "model",color = "model")+
  geom_signif(comparisons = list(c("NNalign","rf"),
                                 c("RankEL","rf")),
              na.rm = T,
              map_signif_level = T,
              data = AUC_score_pp_long %>% filter(cal_type=="AUC01")%>% filter(type!="rf_tme"),
              test = wilcox.test,
              textsize = 5,
              step_increase = 0.15,
              test.args = list(paired = T, exact = F))  +
  theme(legend.position = "none",axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 14), strip.text.x =element_text(size = 16) )
#ggsave(Fig3new_auc01, file="results/PaperPlots/Fig3/Fig3new_auc01.pdf", height = 4, width = 6 )




# -----------------------------------------------------------------
# fig 3 gather 


pdf(file = 'results/PaperPlots/Fig3/Figure3_ALL.pdf', width = 12, height = 16)
ggdraw() +
  draw_plot(Fig3A ,  .0, .45, .4, 0.35) +
  draw_plot(Fig3C,       .4, .61, 0.59, 0.18) +
  draw_plot(Fig3C_2 ,     .4, .45, 0.59, 0.18) +
  draw_plot(Fig3C_legend, .64, .4, 0.1,0.1)  +
  draw_plot(Fig3D, .0, .21, 0.4, 0.23) +
  draw_plot(Fig3B, .4, .09, 0.6, 0.35) +
  draw_plot(Fig_3F_prediction_rf, .25, .0, .25, 0.2) +
  draw_plot(Fig_3F_RankEL,        .0, .0, .25, 0.2) +
  draw_plot_label(c("A","B","C","D","E","F"), x = c(0,0,.4,.0,0.4,0), y = c(1,0.8,0.8,0.45,0.45,0.2), size = 22)  

dev.off()



pdf(file = 'results/PaperPlots/Fig3/SupplementrayFigure3_ALL.pdf', width = 12, height = 16)
ggdraw() +
  draw_plot(SupFig3A ,   .0, .61, 0.70, 0.39) +
  draw_plot(SupFig3B, .65, .88, 0.3, 0.12) +
  draw_plot(SupFig3C, .65, .74, 0.3, 0.152) +
  draw_plot(SupFig3D, .0, .405, 1,0.2)  +
  draw_plot(SupFig3E, .0, .216, 1, 0.2) +
  draw_plot(SupFig3F, .0, .03, 1, 0.2) +
  draw_plot(SupFig3Elegend, .4, -.07, .25, 0.2) +
  
  draw_plot_label(c("A","B","C","D","E","F","G"), x = c(0,0.64,.64,.64,.0,0,0), y = c(-0.1,-0.1,0.9,0.76,0.617,0.427,0.237), size = 22)  

dev.off()



# ------------------------------------------------------------------------------
#                 Figure 4 -- TME 
# ------------------------------------------------------------------------------

Fig4A <- plot_response(y_val ='HLAexp', y_label = "HLAexp",log = TRUE, limit = c(4,20000), map_level = T)

min(all_peptides)
Fig4C <- plot_response(y_val ='CYT', y_label = "CYT",log = TRUE, limit = c(0.1,300), map_level = T)
Fig4D <- plot_response(y_val ='MCPmean', y_label = "MCPmean",log = TRUE, limit = c(0.1,3500), map_level = F)

# 
# CYT_nubmer_response_cor <- cor(x = TME_response_sample$CYT,
#                                y = TME_response_sample$number_responses ,
#                                method = 'spearman')
# Fig4A <- TME_response_sample %>%
#   ggplot(. ,aes( x = CYT , y = number_responses)) + #,label = Patient
#   geom_point(aes(color = cohort)) +
#   #  geom_text(hjust=0, vjust=0) +
#   geom_smooth(method = "lm",  colour = '#525252') +
#   scale_color_manual(breaks = c("Basket","Melanoma","mUC"), values = cohort_col) +
#   scale_x_log10() +
#   theme_bw() +
#   #facet_grid(.~cohort)+
#   labs(x = "CYT", y = "Immunogenic \n neopeptides") +
#   annotate('text', label = paste('Spearman =', round(CYT_nubmer_response_cor, digits = 3)), y = 28, x = 0.5, size = 4) +
#   theme(legend.position = "none",axis.title = element_text(size = 14))
# #ggsave(Fig4A, file="results/PaperPlots/Fig4/Fig4A.pdf", height =4 , width = 5.5 )
# 
# 
# 
# HLA_nubmer_response_cor <- cor(x = TME_response_sample_hla$HLAexp,
#                                y = TME_response_sample_hla$number_responses_pr_HLA ,
#                                method = 'spearman')


# 
# Fig4B <- TME_response_sample_hla %>%
#   ggplot(. ,aes( x = HLAexp , y = number_responses_pr_HLA)) + #,label = Patient
#   geom_point(aes(color = cohort)) +
#   #  geom_text(hjust=0, vjust=0) +
#   geom_smooth(method = "lm",  colour = '#525252') +
#   scale_color_manual(breaks = c("Basket","Melanoma","mUC"), values = cohort_col) +
#   scale_x_log10() +
#   theme_bw() +
#   #facet_grid(.~cohort)+
#   labs(x = "HLA expression", y = "Immunogenic \n neopeptides", color = "Cohort") +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12 ),
#         legend.title = element_text(size = 12 ),
#         axis.title = element_text(size = 14)) +
#   guides(colour = guide_legend(override.aes = list(size=3))) + #title.position = "top"
#   annotate('text', label = paste('Spearman =', round(HLA_nubmer_response_cor, digits = 3)), y = 9, x = 30, size = 4)
# #ggsave(Fig4B, file="results/PaperPlots/Fig4/FigFig4B.pdf", height =4 , width = 5.5 )

# fig4A_B <- get_legend(Fig4B)
# Fig4B <- Fig4B + theme(legend.position = "none")
# 
# # pdf(file = 'results/PaperPlots/Fig4/Figure4B_B.pdf', width = 5, height = 7)
# # ggdraw() +
# #   draw_plot(Fig4A, .0, .06, 1, 0.47) +
# #   draw_plot(Fig4B, .0, .53, 1, 0.47) +
# #   draw_plot(fig4A_B , .0, -.46, 1, 1) 
# # 
# # dev.off()



####################### RF model TME ###########################################3

## random forrest TME

rocobj_rf_TME <- roc(Pred_Modelling$response, Pred_Modelling$prediction_rf_tme, direction = "<")
#define object to plot and calculate AUC
auc_rf_TME <- round(auc(Pred_Modelling$response, Pred_Modelling$prediction_rf_tme),3)
# auc 0.1
pred <- prediction(Pred_Modelling$prediction_rf_tme,Pred_Modelling$response)
auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
auc_rf_TME_0.1 <- round(auc_01@y.values[[1]],4)


Fig4E <- ggroc(list(rf = rocobj_rf, 
                    rf_TME = rocobj_rf_TME)) +
  # ggtitle(paste0('(AUC IMPROVE = ', auc_rf, ')',
  #                ' (AUC IMPROVE TME = ', auc_rf_TME , ')\n',
  #                '(AUC01 IMPROVE = ', auc_rf_0.1, ')',
  #                ' (AUC01 IMPROVE TME = ', auc_rf_TME_0.1,')')) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0),
               color="darkgrey", linetype="dashed") +
  scale_color_manual(values = model_col_rf_TME ,
                     breaks = c("rf_TME","rf"), 
                     labels = c("IMPROVE TME","IMPROVE")) +
  theme_bw() +
  labs(color  = "Model", y="Sensitivity (TPR)", x="Specificity")+
  theme(plot.title = element_text(size = 14),legend.position = "bottom",
        legend.text = element_text(size =14 ),
        axis.title = element_text(size = 18 ),
        legend.title = element_text(size = 16 )) + 
  guides(colour = guide_legend(override.aes = list(size=4), ncol=1, title.position = "top")) 
#ggsave(Fig4C, file = "results/PaperPlots/Fig4/Fig4C.pdf", width = 5 , height = 3.8 )
#ggsave(Fig4C, file = "results/PaperPlots/Fig4/Fig4C.png", width = 5 , height = 3.8 )

roc.test(rocobj_rf, rocobj_rf_TME)


Fig4_test <- Pred_Modelling %>%
  ggplot(., aes(x = response_lab, y = prediction_rf_tme)) +
  geom_quasirandom(aes(color = response_lab),size = 3) +
  geom_boxplot(aes(fill = response_lab), alpha = 0.5)+
  scale_color_manual(breaks = c("no","yes"), values = response_col) +
  geom_signif(comparisons = list(c("no","yes")),
              na.rm = T,
              map_signif_level = F,
              textsize = 5,
              data = Pred_Modelling,
              test = wilcox.test)  +
  scale_fill_manual(breaks = c("no","yes"), values = response_col ) +
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "", y = "IMPROVE TME prediction score", color = "Immunogenic", fill = "Immunogenic")


nrow(rocobj_rf)
Fig4F <- AUC_score_pp_long %>% filter(cal_type=="AUC01") %>% 
  filter(type %in% c("rf","rf_tme")) %>% 
  ggplot(., aes(x = type, y= value, group =type )) +
  geom_boxplot(aes(x = type, y= value, group = type , fill = type), alpha = 0.8) +
  geom_beeswarm(aes(x = type, y= value, group =type , color = type))+
  geom_line(aes(group = Patient),alpha = 0.1) + 
  #facet_grid(.~cal_type)+
  theme_bw() +
  scale_x_discrete(breaks = c("rf","rf_tme"),
                   labels = c("IMPROVE","IMPROVE TME")) +
  scale_y_continuous(limits = c(0,0.10)) +
  scale_color_manual(values = model_col_rf_TME ,
                     breaks = c("rf_tme","rf"),
                     labels = c("IMPROVE TME","IMPROVE")) +
  scale_fill_manual(values = model_col_rf_TME ,
                    breaks = c("rf_tme","rf"),
                    labels = c("IMPROVE TME","IMPROVE")) +
  labs(x="",y="AUC 10% score", fill = "model",color = "model")+
  geom_signif(comparisons = list(c("rf","rf_tme")),
              na.rm = T,
              map_signif_level = F,
              data = AUC_score_pp_long %>% filter(cal_type=="AUC01") %>% filter(type %in% c("rf","rf_tme")),
              test = wilcox.test,
              textsize = 5,
              step_increase = 0.1,
              test.args = list(paired = T, exact = F))  +
  theme(legend.position = "none",axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), strip.text.x =element_text(size = 14) )
#ggsave(Fig4F, file="results/PaperPlots/Fig4/Fig4F.pdf", height = 4, width = 5 )




colnames(Feature_importance_TME) <- c("Feature","importance","Partition")
# bar plot
Feature_importance_TME %>% group_by(Feature) %>% summarise()

feature_sum_tme  <- Feature_importance_TME %>% group_by(Partition,Feature) %>% summarize(mean_imp_part = mean(importance, na.rm=TRUE))

feature_sum_tme<- feature_sum_tme %>% mutate(mean_imp =  mean_imp_part/5)



feature_sum_tme <- feature_sum_tme %>% 
  mutate(Feature_category = 
                    case_when(Feature %in% c("Aro","CysRed","PropHydroAro","Inst","HydroCore",
                                             "mw","pI","pI","PropAcidic","PropAro","PropAromatic","PropBasic",
                                             "PropSmall","Prime","HydroAll") ~ "Physicochemical properties",
                              Feature %in% c("CelPrev","Expression","VarAlFreq") ~ "Mutation qualities",
                              Feature %in% c("RankEL","RankBA","Stability","Stability","NetMHCExp") ~ "Peptide-MHC",
                              Feature %in% c("SelfSim","Foreigness","DAI") ~ "Comparing normal",
                              Feature=="PrioScore" ~ "Other",
                              TRUE ~ "Patient TME"))



Fig4G <-  feature_sum_tme %>% 
  ggplot(., aes(x = reorder(Feature, -mean_imp), y = mean_imp)) +
  geom_bar(aes(fill = Feature_category), stat = "identity")+
  theme(axis.text.x = element_text(angle = 90))  +
  theme_bw() +
  scale_fill_manual(breaks = c("Physicochemical properties","Peptide-MHC","Mutation qualities","Comparing normal","Other","Patient TME"),
                    values = c("#abdda4","#d7191c","#fdae61","#2b83ba","#ffffbf","#df65b0"))  +
  
  labs(x = "", y = "Mean importance", fill = "Feature category") +
  scale_x_discrete(labels=Abbreviations) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size=14), 
        legend.position = "bottom",
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14), 
        axis.title = element_text(size = 14)) +
  guides(fill=guide_legend(ncol=3, title.position = "top"))
#ggsave(Fig4D, file = "results/PaperPlots/Fig4/Fig4D.png", width = 6 , height = 5)



pdf(file = 'results/PaperPlots/Fig4/Figure4_ALL.pdf', width = 11, height = 8)
ggdraw() +
  draw_plot(Fig4A, .0, .7, 0.2, 0.3) +
 # draw_plot(Fig4B, .2, .76, 0.4, 0.24) +
  draw_plot(Fig4C, .2, .7, 0.2, 0.3) +
  draw_plot(Fig4D, .4, .7, 0.2, 0.3) +
#  draw_plot(Fig4B_legend, .3, .7, 0.2, 0.24) +
  
  draw_plot(Fig4E, .0, .1, 0.38, 0.6) +
  draw_plot(Fig4F, .6, .7, .4, 0.3) +
  draw_plot(Fig4G ,.38, .0, .6, 0.72) +
  
  draw_plot_label(c("A","B","C","D","E","F"), x = c(0,0.2,0.4,0,0.59,0.38), y = c(1.005,1.005,1.005,0.74,1.005,0.74), size = 20)  

dev.off()



# -----------------------------------------------------------------

# data prep TME
all_peptides_temt <- all_peptides %>% 
  mutate(sample_hla = paste(Sample,HLA_allele)) %>% 
  group_by(Sample) %>% 
  add_tally(name = "number_screened") %>% 
  filter(response=="1") %>% 
  group_by(Sample,response) %>% 
  add_tally(name = "number_responses") %>% 
  group_by(Sample,response,HLA_allele) %>% 
  add_tally(name = "number_responses_pr_HLA") %>% 
  mutate(fraction_response = number_responses/number_screened) %>% 
  distinct(sample_hla, .keep_all = T) %>% 
  ungroup() %>% 
  select(sample_hla,Sample,HLA_allele,Patient,number_screened,number_responses,number_responses_pr_HLA,fraction_response, cohort)

TME_cols <- c('CYT','Tcells','TcellsCD8', 'CytoxLympho','Blinage','NKcells',
              'MyeloidDC','Monocytes','Neutrophils','Endothelial' ,'Fibroblasts',
              'HLAexp','MCPmean')

TME_response_sample_hla <- all_peptides %>% ungroup() %>% 
  mutate(sample_hla = paste(Sample,HLA_allele)) %>% 
  distinct(sample_hla, .keep_all = T) %>% 
  left_join(., all_peptides_temt) %>% 
  dplyr::select(sample_hla,number_responses,number_screened,number_responses_pr_HLA,TME_cols,Sample, Patient, cohort) 

table(TME_response_sample_hla$Sample)


HLA_summarize <- TME_response_sample_hla %>%  
  group_by(Sample) %>% 
  summarize(HLA_exp_mean = mean(HLAexp, na.rm=TRUE)) %>% 
  distinct(Sample, .keep_all = T) 

TME_response_sample <- TME_response_sample_hla %>% 
  distinct(Sample, .keep_all=T) %>% 
  left_join(., HLA_summarize) 

length(unique(TME_response_sample$HLA_exp_mean))


TME_response_sample$number_responses[is.na(TME_response_sample$number_responses)==T] <- 0
TME_response_sample$number_responses_pr_HLA[is.na(TME_response_sample$number_responses_pr_HLA)==T] <- 0
TME_response_sample_hla$number_responses[is.na(TME_response_sample_hla$number_responses)==T] <- 0
TME_response_sample_hla$number_responses_pr_HLA[is.na(TME_response_sample_hla$number_responses_pr_HLA)==T] <- 0

#######sup figure 4 

Correlation_TME  <-TME_response_sample_hla %>% ungroup() %>% select(TME_cols) %>% select(!c("CYT","HLAexp","MCPmean")) #,HLA_exp_sum) %>% select(-HLAexp)
#Correlation_TME  <-TME_response_sample %>% ungroup() %>% select(TME_cols,HLA_exp_mean) %>% select(-HLAexp)
#Correlation_TME  <- all_peptides %>% ungroup() %>% select(TME_cols) #%>% select(!c("CYT","HLAexp"))

Correlation_TME <- sapply(Correlation_TME,as.numeric)
table(all_peptides$Sample,all_peptides$response)

SupFig4A <- ggcorr(Correlation_TME , method = c("pairwise", "spearman"), 
                   hjust = 0.8, size = 5, label = TRUE, label_size = 3,digits = 2,legend.position = "bottom") #label_round = 2 
SupFig4A



SupFig4B_1 <- plot_response(y_val ='Monocytes', y_label = "Monocytes",log = TRUE, limit = c(0.1,130))
SupFig4B_2 <- plot_response(y_val ='Blinage', y_label = "Blinage",log = TRUE, limit = c(0.1,70000))
SupFig4B_3<- plot_response(y_val ='TcellsCD8', y_label = "TcellsCD8", log = F,  limit = c(0,55))
SupFig4B_4 <- plot_response(y_val ='MyeloidDC', y_label = "MyeloidDC", log = F,  limit = c(0,23))
SupFig4B_5 <- plot_response(y_val ='Tcells', y_label = "Tcells", log = F,  limit = c(0.,40))
SupFig4B_6 <- plot_response(y_val ='NKcells', y_label = "NKcells", log = F,  limit = c(0.,18))
SupFig4B_7 <- plot_response(y_val ='Fibroblasts', y_label = "Fibroblasts", log = T,  limit = c(0.1,4500))
SupFig4B_8 <- plot_response(y_val ='Endothelial', y_label = "Endothelial", log = F,  limit = c(0,30))
SupFig4B_9 <- plot_response(y_val ='Neutrophils', y_label = "Neutrophils", log = F,  limit = c(0,70))
SupFig4B_10 <- plot_response(y_val ='CytoxLympho', y_label = "CytoxLympho", log = F,  limit = c(0,40))




SupFig4C <- all_peptides %>%
  ggplot(., aes(x = response_lab, y = HLAexp)) +
  geom_quasirandom(aes(color = response_lab),size = 3) +
  geom_boxplot(aes(fill = response_lab), alpha = 0.5)+
  scale_y_log10( limits = c(4,18000)) +
  scale_color_manual(breaks = c("no","yes"), values = response_col) +
  facet_grid(.~HLA_type) +
  geom_signif(comparisons = list(c("no","yes")),
              na.rm = T,
              map_signif_level = F,
              textsize = 5,
              data = all_peptides,
              test = wilcox.test)  +
  scale_fill_manual(breaks = c("no","yes"), values = response_col ) +
  theme_bw()+
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "", y = "HLAexp", color = "Immunogenic", fill = "Immunogenic")


Correlation_TME  <- TME_response_sample_hla %>% ungroup() %>% select(CYT,MCPmean,HLAexp)  
Correlation_TME <- sapply(Correlation_TME,as.numeric)
table(all_peptides$Sample,all_peptides$response)


SupFig4D <- ggcorr(Correlation_TME , method = c("pairwise", "spearman"), 
                   hjust = 0.8, size = 5, label = TRUE, label_size = 4, digits = 2, legend.position = "bottom") #label_round = 2 
SupFig4D


# ------------------------------------

improve_tme_dat <- get_pr_dat(model = "improve_tme", predictions = Pred_Modelling$prediction_rf_tme, labels  = Pred_Modelling$response)
pr_all_dat <- bind_rows(pr_all_dat, improve_tme_dat) 
pr_all_dat <- na.omit(pr_all_dat)

PRAUC(Pred_Modelling$prediction_rf_tme, Pred_Modelling$response)
PRAUC(Pred_Modelling$prediction_rf, Pred_Modelling$response)
SupFig4E <- 
  pr_all_dat %>%   filter(model %in% c("improve","improve_tme")) %>% 
  ggplot(. , aes(x =recall, y = percision )) +
  geom_point(aes(color = model))  +
  geom_line(aes(color = model)) +
  scale_x_log10() + 
  scale_y_continuous(limits = c(0,0.3)) + 
  theme_bw() +
  labs(color = "Model")+
  scale_color_manual(values = model_col_rf_TME ,
                     breaks = c("improve","improve_tme"),
                     labels = c("IMPROVE","IMPROVE TME")) +
  theme(legend.position = "bottom")

SupFig4F <- AUC_score_pp_long %>%
  filter(cal_type=="AUC") %>% filter(type %in% c("rf","rf_tme")) %>% 
  ggplot(., aes(x = type, y= value, group =type )) +
  geom_boxplot(aes(x = type, y= value, group = type , fill = type), alpha = 0.8) +
  geom_beeswarm(aes(x = type, y= value, group =type , color = type))+
  geom_line(aes(group = Patient),alpha = 0.1) + 
  facet_grid(.~cal_type)+
  theme_bw() +
  scale_x_discrete(breaks = c("rf","rf_tme"),
                   labels = c("IMPROVE","IMPROVE TME")) +
  scale_y_continuous(limits = c(0.4,1.1)) +
  scale_color_manual(values = model_col_rf_TME ,
                     breaks = c("rf_tme","rf"),
                     labels = c("RF TME","RF")) +
  scale_fill_manual(values = model_col_rf_TME ,
                    breaks = c("rf_tme","rf"),
                    labels = c("IMPROVE TME","IMPROVE")) +
  labs(x="",y="AUC score", fill = "model",color = "model")+
  geom_signif(comparisons = list(c("rf","rf_tme")),
              na.rm = T,
              map_signif_level = F,
              data = AUC_score_pp_long %>% filter(cal_type=="AUC") %>% filter(type %in% c("rf","rf_tme")),
              test = wilcox.test,
              step_increase = 0.1,
              textsize = 5,
              test.args = list(paired = T, exact = F))  +
  theme(legend.position = "none", axis.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 14), strip.text.x =element_text(size = 16)  )


# sub 4 H 
Pred_Modelling <-Pred_Modelling  %>% 
  mutate(delta_rf = prediction_rf_tme-prediction_rf) %>% 
  mutate(improved_tme_pred = case_when(delta_rf>0  ~ "improved", 
                                       TRUE ~ "not-improved"))


Pred_Modelling_immunogenic <- Pred_Modelling %>% filter(response==1)
Pred_Modelling_nonimmunogenic <- Pred_Modelling %>% filter(response==0)
HLA_exp_cor_1 <- cor(x = Pred_Modelling_immunogenic$CYT, 
                     y = Pred_Modelling_immunogenic$delta_rf ,
                     method = 'spearman')

HLA_exp_cor_0 <- cor(x = Pred_Modelling_nonimmunogenic$CYT, 
                     y = Pred_Modelling_nonimmunogenic$delta_rf ,
                     method = 'spearman')


SupFig4G_0 <- Pred_Modelling_nonimmunogenic  %>% 
  # filter(response==1) %>% 
  ggplot(., aes(y = CYT, x = delta_rf )) + 
  annotate('text', label = paste('Spearman =', round(HLA_exp_cor_0, digits = 2)), y = 80, x = -0.12,size = 5) +
  geom_point(aes(color = improved_tme_pred), alpha = 0.2) + 
  geom_smooth(method = "lm",  colour = '#525252') +
  scale_color_manual(breaks = c("improved","not-improved"),
                     values = c("#5ab4ac","#d8b365") , labels = c("higher in IMPROVE TME","lower in IMPROVE TME")) +
  labs(color = "Delta score", y = "CYT", x = "Delta prediction score (RF TME - RF)") +
  #  geom_point(aes(y=prediction_rf_TME), color = "blue", alpha = 0.2) + 
  scale_y_log10() + 
  facet_grid(. ~response_lab,labeller = as_labeller(c("yes" = "Immunogenic",
                                                      "no" = "Non-immunogenic"))) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

SupFig4G_1 <- Pred_Modelling_immunogenic  %>% 
  # filter(response==1) %>% 
  ggplot(., aes(y = CYT, x = delta_rf )) + 
  annotate('text', label = paste('Spearman =', round(HLA_exp_cor_1, digits = 2)), y = 80, x = -0.08, size = 5) +
  geom_point(aes(color = improved_tme_pred), alpha = 0.2) + 
  geom_smooth(method = "lm",  colour = '#525252') +
  scale_color_manual(breaks = c("improved","not-improved"),
                     values = c("#5ab4ac","#d8b365") , labels = c("Higher in RF TME","Lower in RF TME")) +
  labs(color = "Delta score", y = " ", x = "Delta prediction score (RF TME - RF)") +
  #  geom_point(aes(y=prediction_rf_TME), color = "blue", alpha = 0.2) + 
  scale_y_log10() + 
  facet_grid(. ~response_lab,labeller = as_labeller(c("yes" = "Immunogenic",
                                                      "no" = "Non-immunogenic"))) +
  theme_bw() +
  theme(legend.position = "bottom", axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
#ggsave(SupFig4G, file = "results/PaperPlots/Fig4/SupFig4G.png", width = 6 , height = 4 )

SupFig4G_1_legend <- get_legend(SupFig4G_1)
SupFig4G_1 <- SupFig4G_1 + theme(legend.position = "none")


pdf(file = 'results/PaperPlots/Fig4/SupplementrayFigure4_ALL.pdf', width = 14, height = 12)
ggdraw() +
  draw_plot(SupFig4B_1, 0.0 ,0.82, 0.12, 0.18) +
  draw_plot(SupFig4B_2, 0.12 ,0.82, 0.12, 0.18) +
  draw_plot(SupFig4B_3, 0.24 ,0.82, 0.12, 0.18) +
  draw_plot(SupFig4B_4, 0.36 ,0.82, 0.12, 0.18) +
  draw_plot(SupFig4B_5, 0.48 ,0.82, 0.12, 0.18) +
  
  draw_plot(SupFig4B_6, 0.0 ,0.64, 0.12, 0.18) +
  draw_plot(SupFig4B_7, 0.12 ,0.64, 0.12, 0.18) +
  draw_plot(SupFig4B_8, 0.24 ,0.64, 0.12, 0.18) +
  draw_plot(SupFig4B_9, 0.36 ,0.64, 0.12, 0.18) +
  draw_plot(SupFig4B_10, 0.48 ,0.64, 0.12, 0.18) +
  
  draw_plot(SupFig4A, .63 ,0.58, 0.4, 0.44) +
  
  
  draw_plot(SupFig4C , .0, .3, 0.5,0.32)  +
  draw_plot(SupFig4D , .5, 0.33, 0.2, 0.3) +
  draw_plot(SupFig4E, .7, .3, .3, 0.28) +
  
  draw_plot(SupFig4F, .0, .06, .3, 0.24) +
  draw_plot(SupFig4G_0, .3, .06, .337, 0.24) +
  draw_plot(SupFig4G_1, .64, .06, .337, 0.24) +
  draw_plot(SupFig4G_1_legend , .6, -.02, 0.1, 0.1)  +

  draw_plot_label(c("A","B","C","D","E","F","G"), x = c(-.01,0.61,.0,.5,.7,0,0.31), y = c(1,1,0.62,0.62,0.62,.31,0.31), size = 22)

dev.off()


