
# Useful functions when working with logistic regression

load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")

library(ROCR)
library(grid)
library(broom)
library(caret)
library(tidyr)
library(dplyr)
library(scales)
library(ggplot2)
library(ggthemr)
library(ggthemes)
library(gridExtra)
library(data.table)
library(openxlsx)

# calculate cut of for rf TME
# ----------------------------------------
min(pred_df_TME_include$prediction_rf_tme)
results_best_sensitivity_rf_tme  <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(results_best_sensitivity_rf_tme) <- c("cutoff","sensitivity","speceficity")
for (cut in seq(from = 0.2, to = 0.7, by = 0.005)) {
# print(i)
  cm_info <- ConfusionMatrixInfo( data = Pred_Modelling, predict = "prediction_rf_tme", 
                                  actual = "response", cutoff = cut )

  CONFMAT <- as.data.frame(table(cm_info$data$type))
  
  sensitivity = CONFMAT$Freq[CONFMAT$Var1=="TP"]/(CONFMAT$Freq[CONFMAT$Var1=="TP"]+CONFMAT$Freq[CONFMAT$Var1=="FN"])
  speceficity = CONFMAT$Freq[CONFMAT$Var1=="TN"]/(CONFMAT$Freq[CONFMAT$Var1=="TN"]+CONFMAT$Freq[CONFMAT$Var1=="FP"])
  
  vec <-  c(cut,sensitivity,speceficity)
  print(vec)
  results_best_sensitivity_rf_tme <-  rbind(results_best_sensitivity_rf_tme,vec)
}





colnames(results_best_sensitivity_rf_tme) <- c("cutoff","sensitivity","speceficity")


Fig5C_cutoff_figure_TME <-
  results_best_sensitivity_rf_tme %>% gather(., key = "type",value="value", -cutoff) %>% 
  ggplot(.,aes(x = cutoff, y = value)) +
  geom_line(aes(color = type), size = 1)+
  theme_bw() +
  scale_y_continuous(breaks = seq(0.1, 1,0.1))+
  scale_x_continuous(breaks = seq(0.2, 0.7, 0.1), limits = c(0.1,0.82))+
  labs(x="IMPROVE TME cutoff", y = "Performance values", color = "Type") +
  scale_color_manual(breaks = c("sensitivity","speceficity"),
                     values = c("#8dd3c7","#bebada") ) +
  theme(legend.position="bottom",axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16) ,
        legend.text = element_text(size = 16)  ) +
   guides(colour = guide_legend(override.aes = list(size=5)))

#ggsave(Fig5C_cutoff_figure_TME , file = "results/PaperPlots/Fig5/Fig5C_cutoff_figure_TME.pdf", width = 10 , height = 3 )
#write.xlsx(results_best_sensitivity, file = "results/PaperPlots/Fig5/tabel_results_best_sensitivity_rf.xlsx", sep = "\t")

# read on the graph specific cutoffs 
cut_off_tme_cross = 0.415
cut_off_tme_90 = 0.335
  
Fig5C_cutoff_figure_TME <- Fig5C_cutoff_figure_TME +   
  geom_vline(xintercept = cut_off_tme_cross) 

Fig5C_cutoff_figure_TME
Fig5C_legend <- get_legend(Fig5C_cutoff_figure_TME)
Fig5C_cutoff_figure_TME <- Fig5C_cutoff_figure_TME + theme(legend.position="none")
# calculate cut of for rf without TME
# ----------------------------------------

results_best_sensitivity_rf <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(results_best_sensitivity_rf) <- c("cutoff","sensitivity","speceficity")
for (cut in seq(from = 0.217, to = 0.73, by = 0.005)) {
  # print(i)
  cm_info <- ConfusionMatrixInfo( data = Pred_Modelling, predict = "prediction_rf", 
                                  actual = "response", cutoff = cut )
  
  CONFMAT <- as.data.frame(table(cm_info$data$type))
  
  sensitivity = CONFMAT$Freq[CONFMAT$Var1=="TP"]/(CONFMAT$Freq[CONFMAT$Var1=="TP"]+CONFMAT$Freq[CONFMAT$Var1=="FN"])
  speceficity = CONFMAT$Freq[CONFMAT$Var1=="TN"]/(CONFMAT$Freq[CONFMAT$Var1=="TN"]+CONFMAT$Freq[CONFMAT$Var1=="FP"])
  
  vec <-  c(cut,sensitivity,speceficity)
  print(vec)
  results_best_sensitivity_rf <-  rbind(results_best_sensitivity_rf,vec)
}
colnames(results_best_sensitivity_rf) <- c("cutoff","sensitivity","speceficity")



Fig5C_cutoff_figure <-
  results_best_sensitivity_rf %>% gather(., key = "type",value="value", -cutoff) %>% 
  ggplot(.,aes(x = cutoff, y = value)) +
  geom_line(aes(color = type), size = 1)+
  theme_bw() +
  scale_y_continuous(breaks = seq(0.1, 1,0.1))+
  scale_x_continuous(breaks = seq(0.1, 0.8, 0.1), limits = c(0.1,0.82))+
  labs(x="IMPROVE cutoff", y = "Performance values", color = "Type") +
  scale_color_manual(breaks = c("sensitivity","speceficity"),
                     values = c("#8dd3c7","#bebada") ) +
  theme(legend.position="bottom",axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16) ,
        legend.text = element_text(size = 16)  ) +
  guides(colour = guide_legend(override.aes = list(size=5)))

#ggsave(p_rf, file = "results/PaperPlots/Fig5/Fig5B_cutoff_rf_figure.pdf", width = 5 , height = 3 )
cut_off_cross = 0.412

Fig5C_cutoff_figure_legend <- get_legend(Fig5C_cutoff_figure)
Fig5C_cutoff_figure <- Fig5C_cutoff_figure + 
  theme(legend.position =  "none") +
  geom_vline(xintercept = cut_off_cross) 
Fig5C_cutoff_figure

Pred_Modelling$RankEL_minus
results_best_sensitivity  <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(results_best_sensitivity) <- c("cutoff","sensitivity","speceficity")
#for (cut in seq(from = 0.01, to = 0.99, by = 0.005)) {
  for (cut in seq(from = -3, to = -0.0001, by = 0.05)) {
  # print(i)
  cm_info <- ConfusionMatrixInfo( data = Pred_Modelling, predict = "RankEL_minus", 
                                  actual = "response", cutoff = cut )
  
  CONFMAT <- as.data.frame(table(cm_info$data$type))
  
  sensitivity = CONFMAT$Freq[CONFMAT$Var1=="TP"]/(CONFMAT$Freq[CONFMAT$Var1=="TP"]+CONFMAT$Freq[CONFMAT$Var1=="FN"])
  speceficity = CONFMAT$Freq[CONFMAT$Var1=="TN"]/(CONFMAT$Freq[CONFMAT$Var1=="TN"]+CONFMAT$Freq[CONFMAT$Var1=="FP"])
  
  vec <-  c(cut,sensitivity,speceficity)
  print(vec)
  results_best_sensitivity <-  rbind(results_best_sensitivity,vec)
}
colnames(results_best_sensitivity) <- c("cutoff","sensitivity","speceficity")

p_rankel <-
  results_best_sensitivity %>% gather(., key = "type",value="value", -cutoff) %>% 
  ggplot(.,aes(x = cutoff, y = value)) +
  geom_line(aes(color = type), size = 1)+
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1,0.1))+
 # scale_x_continuous(breaks = seq(0.1, 0.8, 0.1), limits = c(0.01,0.99))+
  labs(x="-RankEL cutoff", y = "Performance values", color = "Type") +
  scale_color_manual(breaks = c("sensitivity","speceficity"),
                     values = c("#8dd3c7","#bebada") ) +
  theme(legend.position="none",axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16) ,
        legend.text = element_text(size = 16)  ) +
  guides(colour = guide_legend(override.aes = list(size=5))) 

rankel_cutoff_cross = -0.45
rankel_cutoff_90_speceficity = -1.37

p_rankel <- p_rankel +
geom_vline(xintercept = rankel_cutoff_cross) +
#  geom_vline(xintercept = rankel_cutoff_90_speceficity) +
  geom_vline(xintercept = rankel_cutoff_cross) 

  # save data for figure 
save(cut_off_90,cut_off_cross,cut_off_tme_90,cut_off_tme_cross,
     Fig5C_cutoff_figure_TME,Fig5C_cutoff_figure,
     Fig5C_legend,
     p_rankel,rankel_cutoff_90_speceficity,rankel_cutoff_cross, file="data/04_plotting/Survival_data/05_cut_offs.Rdata")

