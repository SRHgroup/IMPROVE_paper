

load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")
library(grid)
# ------------------------------------------------------------------------------
#                 Figure 2
# ------------------------------------------------------------------------------
colnames(all_peptides)

table(all_peptides$response)
#fig2B
Fig2A_tab <- all_peptides %>% 
  arrange(response) %>% 
  # distinct(Mut_identity, .keep_all = T) %>% 
  mutate(response = as.factor(response)) %>% 
  group_by(Mutation_Consequence,response)  %>% 
  tally()  %>% 
  spread(key = "response", value = "n")

setDT(Fig2A_tab)[, frac_neg := (`0` / sum(`0`))*100]
setDT(Fig2A_tab)[, frac_pos := (`1` / sum(`1`))*100]

Fig2A_tab <- Fig2A_tab %>% as.tibble() %>% 
  dplyr::select(Mutation_Consequence,frac_neg,frac_pos) %>% 
  gather(., key = "type", value = "fraction" , -Mutation_Consequence) %>% 
  mutate(p_val = case_when(Mutation_Consequence == "F" & type == "frac_neg" ~ "p = 0.10",
                           Mutation_Consequence == "M" & type == "frac_neg"~ "p = 0.21",
                           Mutation_Consequence == "I" & type == "frac_neg"~ "p = 0.43",
                           Mutation_Consequence == "D" & type == "frac_neg"~ "p = 0.96"))


Fig2A <- Fig2A_tab  %>% 
  ggplot(aes(x=Mutation_Consequence, y = fraction, fill = type))+
  geom_bar(aes(fill = type), stat = "identity",, position = "dodge")+
  scale_fill_manual(breaks = c("frac_neg","frac_pos"), 
                    values = c("#91bfdb", "#ef8a62") , 
                    labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("M","F","I","D"),
                   labels = c("Missense","Frame-\n shift","Inframe \n insertion","Inframe \n deletion"))+
  theme_bw()+
  scale_y_continuous(limits = c(0,85))+
  # geom_text_repel(aes(label = round(fraction,2)), position = position_stack(vjust = 0.9),direction = "y", 
  #           box.padding = unit(0.01, "lines"))+
  geom_text(aes(label = p_val), vjust = -1) +
  labs(x = "", y = "Percent neopeptides", fill = "Immunogenic") +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 14),
        axis.title.y= element_text(size = 18),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=9))
#ggsave(Fig2A, file="results/PaperPlots/Fig2/Fig2A.pdf", height = 4, width = 4 )

Fig2_legend <- get_legend(Fig2A)
Fig2A <- Fig2A + theme(legend.position = "none")

# stats
# mutation F 
all_peptides$good <- ifelse(all_peptides$Mutation_Consequence == "F",  "good", "bad")
prop_test(all_peptides$good,all_peptides$response,"good","response")
table(all_peptides$response,all_peptides$good)
# mutation M 
all_peptides$good <- ifelse(all_peptides$Mutation_Consequence == "M",  "good", "bad")
prop_test(all_peptides$good,all_peptides$response,"good","response")
table(all_peptides$response,all_peptides$good)
# mutation I
all_peptides$good <- ifelse(all_peptides$Mutation_Consequence == "I",  "good", "bad")
prop_test(all_peptides$good,all_peptides$response,"good","response")
table(all_peptides$response,all_peptides$good)
# mutation D
all_peptides$good <- ifelse(all_peptides$Mutation_Consequence == "D",  "good", "bad")
prop_test(all_peptides$good,all_peptides$response,"good","response")
table(all_peptides$response,all_peptides$good)


######### mutation possition figures 
#source("bin/R_script/04_plotting/source/04_2_Mutation_possition_figure.R")
load("data/04_plotting/Mutation_possistion/mut_figures.Rdata")

######### RNA confirmation 
#source("bin/R_script/04_plotting/source/04_2_validation_mut_rna.R")
load("data/04_plotting/RNA_validation/RNA_val_figuress.Rdata")




# fig 2 B and C in 04_2_mutations_possition 
Fig2D_reponses_CB_IB <- all_peptides %>% # filter(response==1) %>% 
  ggplot(., aes( x = IB_CB_cat, y = SelfSim)) +
  geom_quasirandom(aes(color = IB_CB_cat), alpha = 0.7) +
  geom_boxplot(aes(fill = IB_CB_cat), alpha = 0.7) +
  scale_fill_manual(breaks = c("CB","IB"), labels = c("Conserved","Improved"), values = c("#2ca25f","#88419d"))+
  scale_color_manual(breaks = c("CB","IB"),labels = c("Conserved","Improved"), values = c("#2ca25f","#88419d"))+
  labs(y = "Self-similarity", color = "MHC binding",fill = "MHC binding", x = "MHC binding") +
  geom_signif(comparisons = list(c("CB","IB")),
              na.rm = T,
              map_signif_level = T,
              data = all_peptides, # %>% filter(response==1) ,
              textsize = 6,
              test = wilcox.test,
              test.args=list( p.adjust.method = "bonferroni"))  +
  scale_y_log10(limits = c(0.78,1.03)) + 
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14), 
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "none")  +
  facet_grid(.~response_lab,labeller = as_labeller(c("yes" = "Immunogenic",
                                                     "no" = "Non-immunogenic"))) 
#ggsave(Fig2D_reponses_CB_IB, file="results/PaperPlots/Fig2/Fig2D_reponses_CB_IB.pdf", height = 4, width = 4 )

# figure 2 D in 01_5_validation_mut 

 # figrue E 

cols_selected_AUC <- c('Aro', 'Inst', 'CysRed','RankEL_minus','RankBA_minus','NetMHCExp_minus',
                       'Expression','SelfSim','Prime','PropHydroAro','HydroCore','HydroAll','pI',
                       'minus_PropSmall','PropAro','PropBasic','minus_PropAcidic','DAI','Stability_minus','Foreigness',
                       'CelPrev','PrioScore','VarAlFreq')



roclist <- list()
AUC_score <- data.frame(matrix(ncol = 2, nrow = 0))
AUC_score_01 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(all_peptides)
df <- all_peptides %>% 
  dplyr::select(c(cols_selected_AUC,"response") )
#df <- sapply(df , as.numeric)
df <- as.data.frame(df)

# 
for (x in colnames(df)) {
  print(x)
  roclist[[length(roclist) + 1]] <- roc(df$response, df[,x])
  AUC_score  <- rbind(AUC_score,c(x, auc(df$response, df[,x])))
  names(roclist)[length(roclist)] <- x
  pred <- prediction(df[,x],df$response)
  auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
  AUC_score_01  <- rbind(AUC_score_01,c(x, round(auc_01@y.values[[1]],4)))
  
  
}



colnames(AUC_score_01) <- c("Feature","AUC")
AUC_score_01$AUC <- as.numeric(AUC_score_01$AUC)

colnames(AUC_score) <- c("Feature","AUC")
AUC_score$AUC <- as.numeric(AUC_score$AUC)

AUC_score <- AUC_score  %>% arrange(desc(AUC))
roclist <- roclist[order(match(names(roclist),AUC_score$Feature))]


AUC_score_01  <- AUC_score_01  %>% mutate(Feature_category = 
                                            case_when(Feature %in% c("Aro","CysRed","PropHydroAro","Inst","minus_inst","HydroCore",
                                                                     "mw","minus_pI","pI","minus_PropAcidic","PropAcidic","PropAro","PropAromatic","PropBasic",
                                                                     "PropSmall","minus_PropSmall","Prime","HydroAll") ~ "Physicochemical properties",
                                                      Feature %in% c("CelPrev","Expression","VarAlFreq") ~ "Mutation qualities",
                                                      Feature %in% c("RankEL_minus","RankBA_minus","Stability","Stability_minus","NetMHCExp_minus","NetMHCExp") ~ "Peptide-MHC",
                                                      Feature %in% c("SelfSim","Foreigness","DAI") ~ "Comparing normal",
                                                      Feature=="PrioScore" ~ "Other"))

# Figure 2 I
Fig2J_AUC <- AUC_score_01 %>% filter(!Feature %in% c("response")) %>% 
  ggplot(., aes(x = reorder(Feature,-AUC), y = AUC )) + 
  geom_bar(aes(fill = Feature_category), stat = "identity")+
  scale_fill_manual(breaks = c("Physicochemical properties","Peptide-MHC","Mutation qualities","Comparing normal","Other"),
                    values = c("#abdda4","#d7191c","#fdae61","#2b83ba","#ffffbf"))  +
  scale_y_continuous(breaks = seq(0, 0.02, by = 0.001)) +
  theme_bw()+
  scale_x_discrete(labels=Abbreviations) +
  geom_hline(yintercept = 0.005, linetype = "dashed") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.9,vjust=0.4, size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) +
  labs(x = "", fill = "Feature category", y = "AUC01")
#ggsave(Fig2J_AUC, file = "results/PaperPlots/Fig2/Fig2E_AUC.png", width = 10 , height = 5 )

# Figure 2 F,G,H,I
Fig2E <- plot_response(y_val ='Expression', y_label = "Expression")
Fig2F <- plot_response(y_val ='RankBA', y_label = "RankBA",log = TRUE, limit = c(0.001,150))
Fig2G <- plot_response(y_val ='HydroCore', y_label = "HydroCore",log = FALSE, limit = c(-5,5.4))
Fig2H <- plot_response(y_val ='PropHydroAro', y_label = "PropHydroAro",log = FALSE, limit = c(0,1.1))

all_peptides <- all_peptides %>% mutate(Mut_Exp = 
                                          case_when(RankEL >= 0.5 & Expression >= 2 ~ "MH_EH", 
                                                    RankEL >= 0.5 & Expression < 2 ~ "MH_EL",
                                                    RankEL < 0.5 & Expression < 2 ~ "ML_EL",
                                                    RankEL < 0.5 & Expression >= 2 ~ "ML_EH",
                                                    TRUE ~"What") )


fraction_Mut_Exp_tab <- all_peptides %>% 
  group_by(Mut_Exp,response_lab) %>% 
  tally() %>% 
  spread(., key = "response_lab", value = "n") %>% 
  mutate(fration_immunugenic = yes/(yes+no)) 


# Figure 2 J
Fig2J <- all_peptides %>% 
  ggplot(., aes(x = RankEL, y = Expression )) + 
  geom_point(aes(color = response_lab, alpha = response_lab )) + 
  scale_x_log10(breaks = c( 0.001,.01,.1, .5, 2,8,20), 
                minor_breaks = NULL, 
                labels = c( '0.001','0.01','0.1', '0.5', '2','8','20')) +
  scale_y_log10(breaks = c( .01,.1, .5, 2,10, 40,100,300,1000), 
                minor_breaks = NULL, 
                labels = c( '0.01','0.1', '0.5','2' ,'10', '40','100','300','1000'), limits = c(0.01,3600)) + 
  # scale_size_continuous() +
  annotate("text", x=0.002, y=3000, label= "3.24%", size = 5) + 
  annotate("text", x=0.002, y=0.01, label= "2.44%", size = 5) + 
  annotate("text", x=11, y=3000, label= "2.37%", size = 5) + 
  annotate("text", x=11, y=0.01, label= "2.28%", size = 5) + 
  scale_color_manual(breaks = c("no","yes"), values = response_col) +
  geom_vline(xintercept = c(0.5), colour = '#737373', linetype = 2) +
  geom_hline(yintercept = c(2), colour = '#737373', linetype = 2) +
  theme_bw(base_size = 22) +
  # facet_wrap(.~Patient) +
  # stat_ellipse() +
  labs(y = 'Expression', 
       x = 'RankEL',
       #  size = 'Estimated Frequency %',
       color = "Immunogenic",
       alpha = "Immunogenic") + 
  theme(legend.position = "none", axis.title = element_text(size = 20))
  

#ggsave(Fig2J, file="results/PaperPlots/Fig2/Fig2J.pdf", height = 5, width = 5 )
#ggsave(Fig2J, file="results/PaperPlots/Fig2/Fig2J.png", height = 8, width = 8 )

# Z-test 
# expression
all_peptides$good_or_bad_exp_EL <- ifelse(all_peptides$Expression > 2,  "good", "bad")
prop_test(all_peptides$response,all_peptides$good_or_bad_exp_EL, "reponse","good_bad")
table(all_peptides$response,all_peptides$good_or_bad_exp_EL)
# mut rank 
all_peptides$good_or_bad_exp_EL <- ifelse(all_peptides$RankEL < 0.5, "good", "bad")
prop_test(all_peptides$response,all_peptides$good_or_bad_exp_EL, "reponse","good_bad")
table(all_peptides$response,all_peptides$good_or_bad_exp_EL)
# combination
all_peptides$good_or_bad_exp_EL <- ifelse(all_peptides$Expression > 2 & all_peptides$RankEL < 0.5, "good", "bad")
prop_test(all_peptides$response,all_peptides$good_or_bad_exp_EL, "reponse","good_bad")
table(all_peptides$response,all_peptides$good_or_bad_exp_EL)



all_peptides <- all_peptides %>% mutate(Helix_Hydro  = 
                                          case_when(PropHydroAro >= 0.4 & HydroCore >= 0.3 ~ "HeP_HyP", 
                                                    PropHydroAro >= 0.4 & HydroCore < 0.3 ~ "HeP_HyN",
                                                    PropHydroAro < 0.4 & HydroCore < 0.3 ~ "HeN_HyN",
                                                    PropHydroAro < 0.4 & HydroCore >= 0.3 ~ "HeN_HyP",
                                                    TRUE ~"WTF") )


fraction_He_Hp_tab <- all_peptides %>% 
  group_by(Helix_Hydro,response_lab) %>% 
  tally() %>% 
  spread(., key = "response_lab", value = "n") %>% 
  mutate(fration_immunugenic = yes/(yes+no)) 


Fig2K <- all_peptides %>% 
  ggplot(., aes(x = HydroCore, y = PropHydroAro )) + 
  geom_point(aes(color = response_lab, alpha = response_lab)) + 
  # scale_x_log10(breaks = c( 0.001,.01,.1, .5, 2,8,20), 
  #               minor_breaks = NULL, 
  #               labels = c( '0.001','0.01','0.1', '0.5', '2','8','20')) +
  # scale_x_log10(breaks = c( .01,.1, .5, 2,10, 40,100,200,1000),
  #                minor_breaks = NULL,
  #                labels = c( '0.01','0.1', '0.5','2' ,'10', '40','100','200','1000')) +
  #scale_size_continuous() +
  scale_color_manual(breaks = c("no","yes"), values = response_col, labels = c("No","Yes")) +
  scale_alpha_manual(breaks = c("no","yes"), values = c(0.1,0.8), labels = c("No","Yes")) +
  geom_vline(xintercept = c(0.4), colour = '#737373', linetype = 2) +
  geom_hline(yintercept = c(0.3), colour = '#737373', linetype = 2) +
  annotate("text", x=3.8, y=1, label= "4.68%", size = 5) + 
  annotate("text", x=3.8, y=0, label= "3.07%", size = 5) + 
  annotate("text", x=-3.8, y=1, label= "3.00%", size = 5) + 
  annotate("text", x=-3.8, y=0, label= "1.73%", size = 5) + 
  theme_bw(base_size = 22) +
  # facet_wrap(.~Patient) +
  # stat_ellipse() +
  labs(x = 'HydroCore', 
       y = 'PropHydroAro',
       #  size = 'Estimated Frequency %',
       color = "Immunogenic",
       alpha = "Immunogenic")+
  theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size=3)))
#ggsave(Fig2K, file="results/PaperPlots/Fig2/Fig2K.pdf", height = 7, width = 10 )

legend_Fig2K <- get_legend(Fig2K)
Fig2K <- Fig2K + theme(legend.position = "none",axis.title = element_text(size = 20))

# Z-test 
# expression
all_peptides$good_or_bad_exp_EL <- ifelse(all_peptides$PropHydroAro  > 0.4,  "good", "bad")
prop_test(all_peptides$response,all_peptides$good_or_bad_exp_EL, "reponse","good_bad")
table(all_peptides$response,all_peptides$good_or_bad_exp_EL)
# mut rank 
all_peptides$good_or_bad_exp_EL <- ifelse(all_peptides$HydroCore > 0.3, "good", "bad")
prop_test(all_peptides$response,all_peptides$good_or_bad_exp_EL, "reponse","good_bad")
table(all_peptides$response,all_peptides$good_or_bad_exp_EL)
# combination
all_peptides$good_or_bad_exp_EL <- ifelse(all_peptides$PropHydroAro  > 0.4 & all_peptides$HydroCore > 0.3, "good", "bad")
prop_test(all_peptides$response,all_peptides$good_or_bad_exp_EL, "reponse","good_bad")
table(all_peptides$response,all_peptides$good_or_bad_exp_EL)

# -----------------------
# supplementary fig 2 

# supfig2A in 01_5_validation_mut


    # supFigure C


unique(AUC_score$Feature)
AUC_score <- AUC_score %>%   mutate(Feature_category = 
                                      case_when(Feature %in% c("Aro","CysRed","PropHydroAro","Inst","minus_inst","HydroCore",
                                                               "mw","minus_pI","pI","minus_PropAcidic","PropAcidic","PropAro","PropAromatic","PropBasic",
                                                               "PropSmall","minus_PropSmall","Prime","HydroAll") ~ "Physicochemical properties",
                                                Feature %in% c("CelPrev","Expression","VarAlFreq") ~ "Mutation qualities",
                                                Feature %in% c("RankEL_minus","RankBA_minus","Stability","Stability_minus","NetMHCExp_minus","NetMHCExp") ~ "Peptide-MHC",
                                                Feature %in% c("SelfSim","Foreigness","DAI") ~ "Comparing normal",
                                                Feature=="PrioScore" ~ "Other"))




# --------------------------------------------
#             gather fig 2 

pdf(file = 'results/PaperPlots/Fig2/Figure2_ALL.pdf', width = 17, height = 18)
ggdraw() +
  draw_plot(Fig2A,                 .0, .73, .21, 0.22) +
  draw_plot(Fig2B_mut_pos,        .22, .7, .27, 0.27) +
  draw_plot(Fig2C_mut_type,       .5, .73, .24, 0.24) +
  draw_plot(Fig2D_reponses_CB_IB, .74, .73, .26, 0.24) +
#  draw_plot(CB_IB_legend, .67, .56, .3, 0.3) +
  
  draw_plot(Fig2E_RNA, .0, .5, 0.24, 0.19) +
  draw_plot(Fig2E,     .25, .5, 0.17, 0.19) +
  draw_plot(Fig2F ,    .43, .5, 0.17, 0.19) +
  draw_plot(Fig2G,     .62, .5, 0.17, 0.19) +
  draw_plot(Fig2H,     .81, .5, 0.17, 0.19) +
  
  draw_plot(Fig2J_AUC, .0, .2, 1,0.3) + 

  draw_plot(Fig2J, .0, .0, .34, 0.21) +
  draw_plot(Fig2K, .37, .0, .34, 0.21) +
  draw_plot(Fig2_legend, .62, .001, .3, .3) +
  draw_plot_label(c("A","B","C","D"), x = c(0,.22,.5,.74), y = c(0.98,0.98,0.98,0.98), size = 22) + 
  draw_plot_label(c("E","F","G","H","I"), x = c(0,.26,.43,.63,.82), y = c(0.703,0.703,0.703,0.703,0.703), size = 22)  +
  draw_plot_label(c("J","K","L"), x = c(0,.0,.37), y = c(0.51,0.21,0.21), size = 22)  

dev.off()

# --------------------------------

# supplementray figure 2 
# ------------------------------------


# ------------------------
all_peptides <- all_peptides %>% 
  mutate(Self_Similarity_missense = case_when(Mutation_Consequence=="M" ~ SelfSim,
                                              TRUE ~ "NA"))

all_peptides$Self_Similarity_missense <- ifelse(all_peptides$Mutation_Consequence=="M",all_peptides$SelfSim, NA )                       
all_peptides$self_sim_IB <- ifelse(all_peptides$IB_CB_cat=="IB",all_peptides$SelfSim, NA )
all_peptides$self_sim_CB <- ifelse(all_peptides$IB_CB_cat=="CB",all_peptides$SelfSim, NA )


plot1 <- plot_response_sup_fig2(y_val ='Aro', y_label = "Aro",log = FALSE, limit = c(0,0.88))
plot2 <- plot_response_sup_fig2(y_val ='CelPrev', y_label = "CelPrev",log = FALSE, limit = c(0,1.112))
plot3 <- plot_response_sup_fig2(y_val ='CysRed', y_label = "CysRed",log = FALSE, limit = c(0,24500))
plot4 <- plot_response_sup_fig2(y_val ='Inst', y_label = "Inst",log = FALSE, limit = c(-65,380))
plot6 <- plot_response_sup_fig2(y_val ='HydroAll', y_label = "HydroAll",log = FALSE, limit = c(-3.7,4.9))
all_peptides$mw <- as.numeric(all_peptides$mw)
plot8 <- plot_response_sup_fig2(y_val ='mw', y_label = "mw",log = FALSE, limit = c(650,1750))
plot9 <- plot_response_sup_fig2(y_val ='NetMHCExp', y_label = "NetMHCExp",log = TRUE, limit = c(0.001,350))
plot10 <- plot_response_sup_fig2(y_val ='PropAcidic', y_label = "PropAcidic",log = FALSE, limit = c(0,0.90))
plot11 <- plot_response_sup_fig2(y_val ='PropAro', y_label = "PropAro",log = FALSE, limit = c(0,1.2))
plot24 <- plot_response_sup_fig2(y_val ='PropSmall', y_label = "PropSmall",log = FALSE, limit = c(0,1.2))
max(all_peptides$Prime)
plot7 <- plot_response_sup_fig2(y_val ='RankEL', y_label = "RankEL",log = TRUE, limit = c(0.0005,50))
plot12 <- plot_response_sup_fig2(y_val ='Prime', y_label = "Prime",log = FALSE, limit = c(0,0.4))
plot13 <- plot_response_sup_fig2(y_val ='Stability', y_label = "Stability",log = TRUE, limit = c(0.005,270))
plot14 <- plot_response_sup_fig2(y_val ='VarAlFreq', y_label = "VarAlFreq",log = FALSE, limit = c(0,1.2))
min(all_peptides$DAI)
plot15 <- plot_response_sup_fig2(y_val ='DAI', y_label = "DAI",log = TRUE, limit = c(0.0001,2400))
plot15_1 <- plot_response_sup_fig2(y_val ='DAI_IB', y_label = "DAI CB",log = TRUE, limit = c(0.0001,900))
plot15_2 <- plot_response_sup_fig2(y_val ='DAI_CB', y_label = "DAI IB",log = TRUE, limit = c(0.0001,900))

plot16 <- plot_response_sup_fig2(y_val ='Foreigness', y_label = "Foreigness",log = FALSE, limit = c(0.0000000000000000001,1.2))
max(all_peptides$foreignness_score)

plot17 <- plot_response_sup_fig2(y_val ='pI', y_label = "pI",log = FALSE, limit = c(3,15))
plot18 <- plot_response_sup_fig2(y_val ='PrioScore', y_label = "PrioScore",log = FALSE, limit = c(0,120))
plot19 <- plot_response_sup_fig2(y_val ='PropBasic', y_label = "PropBasic",log = FALSE, limit = c(0,1.2))
plot20 <- plot_response_sup_fig2(y_val ='ValMutRNACoef', y_label = "ValMutRNACoef",log = FALSE, limit = c(0,1.1))
plot21 <- plot_response_sup_fig2(y_val ='SelfSim', y_label = "SelfSim",log = FALSE, limit = c(0.75,1.05))
plot22 <- plot_response_sup_fig2(y_val ='self_sim_CB', y_label = "SelfSim CB",log = FALSE, limit = c(0.75,1.05))
plot23 <- plot_response_sup_fig2(y_val ='self_sim_IB', y_label = "SelfSim IB",log = FALSE, limit = c(0.75,1.05))

SubFig2C <- AUC_score %>% filter(!Feature %in% c("response")) %>% 
  ggplot(., aes(x = reorder(Feature,-AUC), y = AUC )) + 
  geom_bar(aes(fill = Feature_category), stat = "identity")+
  scale_fill_manual(breaks = c("Physicochemical properties","Peptide-MHC","Mutation qualities","Comparing normal","Other"),
                    values = c("#abdda4","#d7191c","#fdae61","#2b83ba","#ffffbf"))  +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.05)) +
  theme_bw()+
  scale_x_discrete(labels=Abbreviations) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.9,vjust = 0.4, size = 16),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = "bottom") +
  labs(x = "", fill = "Feature category") +
  guides(fill=guide_legend(ncol=3,title.position ="top"))
#ggsave(SubFig2C , file = "results/PaperPlots/Fig2/SubFig2C.png", width = 8 , height = 5 )


pdf(file = 'results/PaperPlots/Fig2/Sup2.pdf', height = 22, width = 18)
ggdraw() +
  draw_plot(plot15, .0, .85, 0.165, 0.15) +
  draw_plot(plot16 , .165, .85, 0.165, 0.15) +
  draw_plot(plot9 , .33, .85, 0.165, 0.15) +
  draw_plot(plot17, .495, .85, 0.165, 0.15) +
  draw_plot(plot18, .66, .85,0.165, 0.15) +
  draw_plot(plot19, .825, .85,0.165, 0.15) +
  
  draw_plot(plot21, .0, .7, 0.165, 0.15) +
  draw_plot(plot22 , .165, .7, 0.165, 0.15) +
  draw_plot(plot23, .33, .7, 0.165, 0.15) +
  draw_plot(plot20, .495, .7, 0.165, 0.15) +
  draw_plot(plot14, .66, .7, 0.165, 0.15) +
  
  draw_plot(plot1, .0, .55, 0.165, 0.15) +
  draw_plot(plot2 , .165, .55, 0.165, 0.15) +
  draw_plot(plot3, .33, .55, 0.165, 0.15) +
  draw_plot(plot6, .495, .55, 0.165, 0.15) +
  draw_plot(plot4, .66, .55, 0.165, 0.15) +
  draw_plot(plot8, .825, .55, 0.165, 0.15) +
  
  draw_plot(plot10, .0, .4, 0.165, 0.15) +
  draw_plot(plot11, .165, .4, 0.165, 0.15) +
  draw_plot(plot24, .33, .4, 0.165, 0.15) +
  draw_plot(plot12, .495, .4, 0.165, 0.15) +
  draw_plot(plot7, .66, .4, 0.165, 0.15) +
  draw_plot(plot13 , .825, .4, 0.165, 0.15) +
  
  draw_plot(SupFig2C_RNA_cohort, .0, .12, 0.5, 0.25) +
  draw_plot(SubFig2C, .5, .1, 0.5, 0.28) +
  
  draw_plot_label(c("A","B","C","D"), x = c(0,0,0,.495), y = c(1,0.72,.4,.4), size = 22) 


dev.off()




