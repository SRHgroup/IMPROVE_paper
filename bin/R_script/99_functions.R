## 
# packages 
# ----------

library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggbeeswarm)
library(ggrepel)
library(cowplot)
library(ggpubr)

# functions 

# -----------------------------------------------------------------
##        Z-test 
# -----------------------------------------------------------------

# define color 
response_col <- c("#91bfdb", "#ef8a62")
el_col <- c("#fdae61") 
cohort_col <- c("#f4a582","#c51b7d","#a1d76a")
HLA_col <- c("#5ab4ac","#f6e8c3","#8c510a")

# define color 
model_col_rf_NNalign <- c("#762a83","#7fbf7b") # rf , nnalign
model_col_rf_NNalign_assembling <- c("#762a83","#7fbf7b","#7fcdbb") # rf , nnalign, rf assembling
model_col_rf_TME_NNalign <- c("#e7d4e8","#762a83","#7fbf7b") # rf_tme,rf ,  nnalign
model_col_rf_TME_NNalign_assembling_RankEL <- c("#762a83","#7fbf7b","#7fcdbb","#d7191c") #  rf , nnalign, rf assembling, rankel
model_col_rf_TME_NNalign_RankEL <- c("#e7d4e8","#762a83","#7fbf7b","#d7191c") #  rf_tme,rf ,  nnalign, rankel
model_col_rf_NNalign_RankEL <- c("#762a83","#7fbf7b","#d7191c") #  rf ,  nnalign, rankel
model_col_rf_TME <- c("#e7d4e8","#762a83") # rf_tme,rf ,  nnalign
model_col_f5 <- c("#fdae61","#762a83","#e7d4e8") 

model_col_f5_random <- c("grey","#d7191c","#762a83","#e7d4e8") 


## abrreviations 
# -------------------------------------------
Abbreviations <- c(
  'mw'="mw",
  'Aro'="Aro",
  'Inst'="Inst",
  'minus_inst'="Inst",
  'CysRed'="CysRed",
  'RankEL' = "RankEL",
  'RankBA' = "RankBA",
  'RankBA_minus' = "RankBA",
  'RankEL_minus' = "RankEL",
  'Expression' = "Expression",
  'PrioScore' = "PrioScore",
  'VarAlFreq' = "VarAlFreq",
  'CelPrev' = "CelPrev",
  'Self_Similarity' = "SelfSim",
  'self_sim_IB'= "SelfSimIB",
  'self_sim_CB'= "SelfSimCB",
  'Prime' = "Prime",
  'PropHydroAro' = "PropHydroAro",
  'HydroCore' = "HydroCore",
  'PropSmall' = 'PropSmall',
  'minus_PropSmall' = 'PropSmall',
  'PropBasic' = "PropBasic",
  'PropAcidic' = "PropAcidic",
  'minus_PropAcidic' = "PropAcidic",
  'PropAromatic' = "PropAromatic",
  'pI' = 'pI',
  'minus_pI' = 'pI',
  'DAI' = "DAI",
  'ValMutRNACoef' = "ValMutRNACoef",
  'Stability' = "Stability",
  'Stability_minus' = "Stability",
  'NetMHCExp' = "NetMHCExp",
  'NetMHCExp_minus' = "NetMHCExp",
  'Foreigness' = "Foreigness" ,
  'HydroAll' = "HydroAll",
  "CYT" = "CYT",
  "HLAexp" = "HLAexp",
  "Tcells" = "Tcells",
  "TcellsCD8" = "TcellsCD8" ,
  "CytoxLympho" = "CytoxLympho" ,
  "Blineage" = "Blineage",
  "NKcells" = "NKcells",
  "Monocytes" = "Monocytes",
  "MyeloidDC" = "MyeloidDC",
  "Neutrophils" = "Neutrophils",
  "Endothelial" = "Endothelial",
  "Fibroblasts" = "Fibroblast")


# ---- 
# Theme 
# --------

theme(legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 0.8, size = 10),
      strip.background = element_blank(),
      strip.text = element_text(face="bold", size=9))



#########
prop_test <- function(x,y, x_name, y_name) {
  tab <- as.data.frame(table(x,y))
  x_name_c <-  as.character(x_name)
  y_name_c <- as.character(y_name)
  colnames(tab) <- c(x_name_c,y_name_c,"count")
  tab$count <- as.numeric(tab$count)
  tab <- spread(key = x_name,value = count, data = tab)
  tab[,1] <- NULL
  colnames(tab) <- c(x_name_c,y_name_c)
  rownames(tab) <- c("no","yes")
  tab
  prop.test(tab[,2],rowSums(tab))
  
}




# -----------------------------------------------------------------
##        plot_response 
# -----------------------------------------------------------------

plot_response <- function(data = all_peptides,
                          y_val,
                          log = TRUE,
                          y_label = "",
                          valid = FALSE,
                          y_position = NULL,
                          facet = FALSE,
                          limit = c(0,100),
                          map_level = T) {
 # mini <- data %>% select(y_val) %>% min() %>% as.numeric()
#  maxi <- data %>% select(y_val) %>% max() %>% as.numeric()

  
  p <- data %>% 
    ggplot(., aes_string('response_lab', y_val)) +
    geom_quasirandom(aes(color = response_lab),size = 3) + 
    geom_boxplot(aes(fill = response_lab),alpha = 0.5)+
    scale_color_manual(breaks = c("no","yes"), values = response_col) +
    scale_fill_manual(breaks = c("no","yes"), values = response_col ) +
    theme_bw()+
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x = "", y = y_label, color = "Immunogenic", fill = "Immunogenic") +
    # facet_wrap(Sample ~ .) +
    geom_signif(comparisons = list(c("no","yes")),
                na.rm = T,
                map_signif_level = map_level,
                data = data,
                textsize = 5,
                y_position = y_position,
                test = wilcox.test,
                test.args=list( p.adjust.method = "bonferroni"))  
  if ( log==TRUE) {
    p <- p + scale_y_log10( breaks = c(0.1,0.5,2,10,100,1000,10000), 
                            labels = c("0.1","0.5","2","10","100","1000","10000"), 
                            limits = limit)
  }
  
  if ( log==FALSE) {
    p <- p + scale_y_continuous(limits = limit)
  }
  
  if ( y_val == "Expression") {
    p <- p + scale_y_log10( breaks = c(0.01,0.1,1,10,50,3000), 
                            labels = c("0.01","0.1","1","10","50","3000"),
                            limits = c(0.01,25000)) #
  }
  
  # if (y_val == "Self_Similarity") {
  #   p <- p +   facet_grid(.~ IB_CB_cat) 
  # }
  if (facet == TRUE) {
    p <- p +   facet_grid(.~ HLA_type) 
  }
  

  # if (valid == T) {
  #   ggsave(p, file = paste0("results/PaperPlots/Fig2/", y_label, ".pdf"), width = 5 , height = 4)
  # }
  # ggsave(p, file = paste0("results/PaperPlots/Fig2/", y_label, ".pdf"), width = 5 , height = 4)
  # 
return(p)
}



## sup fig 2 
plot_response_sup_fig2 <- function(data = all_peptides,
                          y_val,
                          log = TRUE,
                          y_label = "",
                          valid = FALSE,
                          y_position = NULL,
                          facet = FALSE,
                          limit = c(0,100)) {
  # mini <- data %>% select(y_val) %>% min() %>% as.numeric()
  #  maxi <- data %>% select(y_val) %>% max() %>% as.numeric()
  
  
  p <- data %>% 
    ggplot(., aes_string('response_lab', y_val)) +
    geom_quasirandom(aes(color = response_lab),size = 3) + 
    geom_boxplot(aes(fill = response_lab),alpha = 0.5)+
    scale_color_manual(breaks = c("no","yes"), values = response_col) +
    scale_fill_manual(breaks = c("no","yes"), values = response_col ) +
    theme_bw()+
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 14)) +
    labs(x = "", y = y_label, color = "Immunogenic", fill = "Immunogenic") +
    # facet_wrap(Sample ~ .) +
    geom_signif(comparisons = list(c("no","yes")),
                na.rm = T,
                map_signif_level = F,
                data = data,
                textsize = 7,
                y_position = y_position,
                test = wilcox.test,
                test.args=list( p.adjust.method = "bonferroni"))  
  if ( log==TRUE) {
    p <- p + scale_y_log10( breaks = c(0.01,0.1,0.5,2,10,50), 
                            labels = c("0.01","0.1","0.5","2","10","50"), 
                            limits = limit)
  }
  
  if ( log==FALSE) {
    p <- p + scale_y_continuous(limits = limit)
  }
  
  if ( y_val == "Expression") {
    p <- p + scale_y_log10( breaks = c(0.01,1,2,5,10,50,3000), 
                            labels = c("0.01","1","2","5","10","50","3000"),
                            limits = c(0.01,NA)) #
  }
  
  # if (y_val == "Self_Similarity") {
  #   p <- p +   facet_grid(.~ IB_CB_cat) 
  # }
  if (facet == TRUE) {
    p <- p +   facet_grid(.~ HLA_type) 
  }
  
  
  # if (valid == T) {
  #   ggsave(p, file = paste0("results/PaperPlots/Fig2/boxplots/", y_label, ".pdf"), width = 5 , height = 4)
  # }
  # ggsave(p, file = paste0("results/PaperPlots/Fig2/boxplots/", y_label, ".pdf"), width = 5 , height = 4)
  # 
  return(p)
}
# -------------------------------------------------


# find pval and hz for survival 
library(broom)

HR_and_pval <- function(
  col = All_samples_clinical_subset$pred_neo_model_cut_TME_high_low,
  surv_object = os) {
  
  fit <- coxph(surv_object ~ col, data = All_samples_clinical_subset)
  df <-  fit |> 
    tidy(conf.int = TRUE, exponentiate = TRUE) |> 
    dplyr::select(term, estimate, starts_with("conf"),p.value)
  HR_survival = paste("HR = ",round(df$estimate,2)," (",round(df$conf.low,2),"-",round(df$conf.high,2),")",sep="" )
  pval_survival = paste("p = ",round(df$p.value,3), sep = "")
  surv_list <- list(pval_survival = pval_survival,
                    HR_survival = HR_survival)
return(surv_list)
  
}




# Survival plot
# ------------------------------------------------------------------
surv_plot <- function(y_fit,                     # survfit object with calculated statistics.
                      inputdata = All_samples_clinical_subset,  
                      legend_lab_names = c("high","low"), 
                      ylab_name = "OS",
                      legend_position = c(.92, .95),
                      xlimits = c(0,100),
                      palette = c('#1f78b4', '#33a02c'),
                      pcoord = c(70, .7),
                      plot_title="",
                      HR = 'HR',
                      hrx = 5,
                      p_value = "NS")
{
  survplot <- ggsurvplot(y_fit,                     # survfit object with calculated statistics.
                         data = inputdata,          # data used to fit survival curves.
                         risk.table = TRUE,         # show risk table.
                         pval = p_value,               # show p-value of log-rank test.
                         pval.coord = pcoord ,
                         xlim = xlimits,         
                         xlab = "Months",   
                         title = plot_title,
                         ylab = ylab_name,
                         break.time.by = 10,         # break X axis in time intervals by 500.
                         ggtheme = theme_classic(base_size = 16),   # customize plot and risk table with a theme.
                         risk.table.y.text.col = T, # colour risk table text annotations.
                         risk.table.y.text = FALSE,  
                         risk.table.fontsize = 6,
                         tables.theme = theme_cleantable(),
                         surv.scale = 'percent',
                         legend = legend_position,
                         legend.title="",
                        legend.labs = legend_lab_names,
                         legend.text = element_text(size = 16),
                         palette = palette)
  # ggtheme(theme_survminer(font.x = c(14, "plain", "black"),
  #                         font.y = c(14, "plain", "black"),
  #                         font.caption = c(15, "plain", "black")))
  # change margins 
  survplot$plot <- survplot$plot + theme(plot.margin = unit(c(3,3,0,0), "lines"), 
                                         axis.title.x = element_text(size = 21),
                                         axis.title.y = element_text(size = 21),
                                         legend.text = element_text(size = 19))
  survplot$table <- survplot$table + theme(plot.title = element_text(size = 18))
  survplot$plot <- survplot$plot+ annotate("text", x = pcoord[1] + hrx, y = pcoord[2]-0.1, label = HR, size = 5)
  return(survplot)
}


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


change_format_hla_allel_BARCODE_DAT <- function(data) {
  data$HLA_allel <- 
    str_c(
      str_sub(data$HLA.x, start = 1L, end = 3L), 
      ":", 
      str_sub(data$HLA.x, start = 4L, end = 5L)
    )
  
  data$HLA_allel <- paste("HLA",data$HLA_allel, sep='-')
  return(data)
}



####### function for patients
pred_per_patient <- function(y_val = predition_rf, c = "mUC",
                             legpos = "none",
                             h = 4) {
  
  p <- Pred_Modelling  %>% 
    filter(cohort==c) %>% 
    ggplot(., aes_string( x = 'response_lab' , y = y_val, fill = 'response_lab'))+
    geom_boxplot(aes(fill = response_lab),alpha = 0.5)+ 
    scale_fill_manual(breaks = c("no","yes"), values = response_col ) +
    labs(y = "Prediction score", x = "", color = "Immunugenic", fill = "Immunogenic") +
    scale_y_continuous(limits = c(0.2,0.9)) +
    geom_signif(comparisons = list(c("yes","no")),
                na.rm = T,
                textsize = 4,
                map_signif_level = F,
                data = Pred_Modelling  %>% 
                  filter(cohort==c),
                test = wilcox.test) +
    theme_bw() +
   # facet_grid(.~Patient) +
    facet_wrap( ~ Patient, nrow = 2)+
    theme(legend.position = legpos, 
          strip.text.x = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 14))
  
  
  return(p)
 # ggsave(p , file=paste0("results/PaperPlots/Fig3/SupFig3D-E_feature_pred_",c,".pdf"), height = h, width = 22)
}

####### Top figures #########3


top_figure_percent <- function(top_number = 20) {
  set.seed(10)
  top_pred <- Pred_Modelling %>% 
    group_by(Patient) %>% 
    top_n(top_number,prediction_rf) %>% 
    group_by(Patient,response) %>% 
    tally() %>% 
    spread(., key = response, value = n ) %>% 
    mutate("number"  = `1`) %>% 
    mutate("cat" = "top_50") %>% 
    mutate("type" = "prediction_rf")
  
  top_pred_TME <- Pred_Modelling %>% 
    group_by(Patient) %>% 
    top_n(top_number,prediction_rf_tme) %>% 
    group_by(Patient,response) %>% 
    tally() %>% 
    spread(., key = response, value = n ) %>% 
    mutate("number"  = `1`) %>% 
    mutate("cat" = "top_50") %>% 
    mutate("type" = "prediction_rf_TME")
  
  top_EL <- Pred_Modelling %>% 
    group_by(Patient) %>% 
    top_n(top_number,-RankEL) %>% 
    group_by(Patient,response) %>% 
    tally() %>% 
    spread(., key = response, value = n ) %>% 
    mutate("number"  = `1`) %>% 
    mutate("cat" = "top_50") %>% 
    mutate("type" = "rank_EL")
  
  Random_top <- Pred_Modelling %>% 
    group_by(Patient) %>% 
    sample_n(., top_number)%>% 
    group_by(Patient,response) %>% 
    tally() %>% 
    spread(., key = response, value = n ) %>% 
    mutate("number"  = `1`) %>% 
    mutate("cat" = "top_50") %>% 
    mutate("type" = "Random")
  
  top_tab <- full_join(top_pred,top_EL) %>% full_join(.,top_pred_TME) %>% full_join(.,Random_top)
  top_tab$number[is.na(top_tab$number)==T] <- 0
  top_tab$type <- factor(top_tab$type, levels=c("Random","rank_EL","prediction_rf","prediction_rf_TME"))
  
  # add number of screened and number of responses 
  
  top_tab <- Pred_Modelling %>% 
    filter(response==1) %>% 
    group_by(Patient) %>% 
    tally(n = "number_response") %>% 
    right_join(.,top_tab ) %>% 
    mutate(fraction_of_response = number/number_response)
  
  
  p1 <- top_tab %>% 
    ggplot(., aes(x = type, y = fraction_of_response)) +
    geom_line(aes(group = Patient),alpha = 0.2) +
    geom_boxplot(aes(fill = type),alpha = 0.5) +
    geom_quasirandom(aes(color = type),size = 3, alpha = 0.6) + 
    scale_x_discrete(labels = c("Random","RankEL","IMPROVE","IMPROVE TME") ,breaks = c("Random","rank_EL","prediction_rf","prediction_rf_TME"))+
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1")) +
    scale_fill_manual(breaks = c("Random","rank_EL","prediction_rf","prediction_rf_TME"), values = model_col_f5_random  ,
                      labels = c("Random","RankEL","IMPROVE","IMPROVE TME") ) +
    scale_color_manual(breaks = c("Random","rank_EL","prediction_rf","prediction_rf_TME"), values = model_col_f5_random  ,
                       labels = c("Random","RankEL","IMPROVE","IMPROVE TME") ) +
    #geom_line(aes(group = Sample)) +
    geom_signif(comparisons = list(c("prediction_rf","rank_EL"),
                                   c("prediction_rf_TME","prediction_rf"),
                                   c("prediction_rf","Random")),
                na.rm = T,
                map_signif_level = F,
                data = top_tab,
                test = wilcox.test,
                step_increase = 0.1,
                test.args = list(paired = T,exact = F)
    ) +
    theme_bw() +
    theme(legend.position="bottom",axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16), axis.title = element_text(size = 18),
          legend.title = element_text(size = 16) ,
          legend.text = element_text(size = 16)  ) +
    labs(x = "", fill = "Model",color = "Model", y = "Fraction of immunugenic neoepitope") + 
    guides(fill=guide_legend(ncol=4,title.position ="left")) +
    ggtitle(paste("Top hit in top:",top_number))
  ggsave(p1, file=paste0("results/PaperPlots/Fig5/top_percent_pval_test",top_number,".pdf"), height = 5, width = 10 )
  
  p2 <- top_tab %>% 
    ggplot(., aes(x = type, y = fraction_of_response)) +
    geom_line(aes(group = Patient),alpha = 0.2) +
    geom_boxplot(aes(fill = type),alpha = 0.5) +
    geom_quasirandom(aes(color = type),size = 3, alpha = 0.6) + 
    scale_x_discrete(labels = c("Random","RankEL","IMPROVE","IMPROVE TME") ,breaks = c("Random","rank_EL","prediction_rf","prediction_rf_TME"))+
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1")) +
    scale_fill_manual(breaks = c("Random","rank_EL","prediction_rf","prediction_rf_TME"), values = model_col_f5_random  ,
                      labels = c("Random","RankEL","IMPROVE","IMPROVE TME") ) +
    scale_color_manual(breaks = c("Random","rank_EL","prediction_rf","prediction_rf_TME"), values = model_col_f5_random  ,
                       labels = c("Random","RankEL","IMPROVE","IMPROVE TME") ) +
    #geom_line(aes(group = Sample)) +
    geom_signif(comparisons = list(c("prediction_rf","rank_EL"),
                                   c("prediction_rf_TME","prediction_rf"),
                                   c("prediction_rf","Random")),
                na.rm = T,
                map_signif_level = T,
                data = top_tab,
                test = wilcox.test,
                step_increase = 0.1,
                test.args = list(paired = T,exact = F)
    ) +
    theme_bw() +
    theme(legend.position="bottom",axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16), axis.title = element_text(size = 18),
          legend.title = element_text(size = 16) , title = element_text(size = 16),
          legend.text = element_text(size = 16)  ) +
    labs(x = "", fill = "Model",color = "Model", y = "Fraction of immunugenic neoepitope") + 
    guides(fill=guide_legend(ncol=4,title.position ="left")) +
    ggtitle(paste("Top hit in top:",top_number))
  ggsave(p2, file=paste0("results/PaperPlots/Fig5/top_percent_no_test_pval",top_number,".pdf"), height = 5, width = 10 )
  
}




## Abbrevi

# AUC pr patient 
#patien="MM-03"
Auc_pp <- function(name = "rf" , var = "prediction_rf") {
  for (patien in unique(p_vec)) {
    df_p <- Pred_Modelling %>% filter(Patient==patien)
    name <- name
    var_pred <- df_p[,var] %>% pull()
    rocobj<- roc(df_p$response, var_pred ,direction = "<")
    #define object to plot and calculate AUC
    auc <- round(auc(df_p$response, var_pred),4)
    # auc 0.1
    pred <- prediction(var_pred,df_p$response)
    auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
    auc_rf_0.1 <- round(auc_01@y.values[[1]],4)
    
    vec <- c(patien,auc,auc_rf_0.1,name)
    AUC_score_pp  <- rbind(AUC_score_pp,vec)
  }
  AUC_score_pp <<- AUC_score_pp
  return(AUC_score_pp)
}



#### for cut-offs  



# Useful functions when working with logistic regression
library(ROCR)
library(grid)
library(caret)
library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(data.table)


# ------------------------------------------------------------------------------------------
# [AccuracyCutoffInfo] : 
# Obtain the accuracy on the trainining and testing dataset.
# for cutoff value ranging from .4 to .8 ( with a .05 increase )
# @train   : your data.table or data.frame type training data ( assumes you have the predicted score in it ).
# @test    : your data.table or data.frame type testing data
# @predict : prediction's column name (assumes the same for training and testing set)
# @actual  : actual results' column name
# returns  : 1. data : a data.table with three columns.
#            		   each row indicates the cutoff value and the accuracy for the 
#            		   train and test set respectively.
# 			 2. plot : plot that visualizes the data.table

AccuracyCutoffInfo <- function( train, test, predict, actual )
{
  # change the cutoff value's range as you please 
  cutoff <- seq( .4, .8, by = .05 )
  
  accuracy <- lapply( cutoff, function(c)
  {
    # use the confusionMatrix from the caret package
    data_train <- as.factor( as.numeric( train[[predict]] > c ) )
    cm_train <- confusionMatrix(data_train, as.factor(train[[actual]]) )
    data_test <- as.factor( as.numeric( test[[predict]] > c ) )
    cm_test  <- confusionMatrix( data_test, as.factor(test[[actual]]) )
    
    dt <- data.table( cutoff = c,
                      train  = cm_train$overall[["Accuracy"]],
                      test   = cm_test$overall[["Accuracy"]] )
    return(dt)
  }) %>% rbindlist()
  
  # visualize the accuracy of the train and test set for different cutoff value 
  # accuracy in percentage.
  accuracy_long <- gather( accuracy, "data", "accuracy", -1 )
  
  plot <- ggplot( accuracy_long, aes( cutoff, accuracy, group = data, color = data ) ) + 
    geom_line( size = 1 ) + geom_point( size = 3 ) +
    scale_y_continuous( label = percent ) +
    ggtitle( "Train/Test Accuracy for Different Cutoff" )
  
  return( list( data = accuracy, plot = plot ) )
}


# ------------------------------------------------------------------------------------------
# [ConfusionMatrixInfo] : 
# Obtain the confusion matrix plot and data.table for a given
# dataset that already consists the predicted score and actual outcome.
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome 
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cutoff  : cutoff value for the prediction score 
# return   : 1. data : a data.table consisting of three column
#            		   the first two stores the original value of the prediction and actual outcome from
#			 		   the passed in data frame, the third indicates the type, which is after choosing the 
#			 		   cutoff value, will this row be a true/false positive/ negative 
#            2. plot : plot that visualizes the data.table 

ConfusionMatrixInfo <- function( data, predict, actual, cutoff )
{	
  # extract the column ;
  # relevel making 1 appears on the more commonly seen position in 
  # a two by two confusion matrix	
  predict <- data[[predict]]
  actual  <- relevel( as.factor( data[[actual]] ), "1" )
  
  result <- data.table( actual = actual, predict = predict )
  
  # calculating each pred falls into which category for the confusion matrix
  result[ , type := ifelse( predict >= cutoff & actual == 1, "TP",
                            ifelse( predict >= cutoff & actual == 0, "FP", 
                                    ifelse( predict <  cutoff & actual == 1, "FN", "TN" ) ) ) %>% as.factor() ]
  
  # jittering : can spread the points along the x axis 
  plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
    geom_violin( fill = "white", color = NA ) +
    geom_jitter( shape = 1 ) + 
    geom_hline( yintercept = cutoff, color = "blue", alpha = 0.6 ) + 
    scale_y_continuous( limits = c( 0, 1 ) ) + 
    scale_color_discrete( breaks = c( "TP", "FN", "FP", "TN" ) ) + # ordering of the legend 
    guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
    ggtitle( sprintf( "Confusion Matrix with Cutoff at %.2f", cutoff ) )
  
  return( list( data = result, plot = plot ) )
}


# ------------------------------------------------------------------------------------------
# [ROCInfo] : 
# Pass in the data that already consists the predicted score and actual outcome.
# to obtain the ROC curve 
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cost.fp : associated cost for a false positive 
# @cost.fn : associated cost for a false negative 
# return   : a list containing  
#			 1. plot        : a side by side roc and cost plot, title showing optimal cutoff value
# 				 	   		  title showing optimal cutoff, total cost, and area under the curve (auc)
# 		     2. cutoff      : optimal cutoff value according to the specified fp/fn cost 
#		     3. totalcost   : total cost according to the specified fp/fn cost
#			 4. auc 		: area under the curve
#		     5. sensitivity : TP / (TP + FN)
#		     6. specificity : TN / (FP + TN)

ROCInfo <- function( data, predict, actual, cost.fp, cost.fn )
{
  # calculate the values using the ROCR library
  # true positive, false postive 
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # cost with the specified false positive and false negative cost 
  # false postive rate * number of negative instances * false positive cost + 
  # false negative rate * number of positive instances * false negative cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 0 ) + 
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 1 )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # optimal cutoff value, and the corresponding true positive and false positive rate
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # normalize the cost to assign colors to 1
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)   
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) + 
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )				
  
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    scale_y_continuous( labels = comma ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )	
  
  # the main title for the two arranged plot
  sub_title <- sprintf( "Cutoff at %.2f - Total Cost = %d, AUC = %.3f", 
                        best_cutoff, best_cost, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2, 
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  return( list( plot 		  = plot, 
                cutoff 	  = best_cutoff, 
                totalcost   = best_cost, 
                auc         = auc,
                sensitivity = best_tpr, 
                specificity = 1 - best_fpr ) )
}



## percision recall 
# ----------------------
get_pr_dat <- function(model = "improve", predictions = Pred_Modelling$prediction_rf, labels  = Pred_Modelling$response) {
  pred <- prediction(predictions, labels)
  perf <- performance(pred, "prec", "rec")
  recall <- as.data.frame(perf@x.values)
  percision <- as.data.frame(perf@y.values)
  pr_dat <- cbind(recall,percision)
  colnames(pr_dat) <- c("recall","percision")
  pr_dat$model <- model
  return(pr_dat)
}


