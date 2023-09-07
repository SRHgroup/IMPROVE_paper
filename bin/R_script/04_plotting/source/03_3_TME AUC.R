# performance pr TME features 
load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")

TME_cols <- c('CYT','Monocytes',
              'Tcells','TcellsCD8', 'CytoxLympho','Blinage','NKcells',
              'MyeloidDC','Neutrophils','Endothelial' ,'Fibroblasts','HLAexp')


# all_peptides_temt <- all_peptides %>% 
#   group_by(Sample) %>% 
#   add_tally(name = "number_screened") %>% 
#   filter(response=="1") %>% 
#   group_by(Sample,response) %>% 
#   add_tally(name = "number_responses") %>% 
#   group_by(Sample,response,HLA_allele) %>% 
#   add_tally(name = "number_responses_pr_HLA") %>% 
#   mutate(fraction_response = number_responses/number_screened)


TME_response <- 
all_peptides %>% 
  mutate(sampe_hla = paste(Sample,HLA_allele)) %>% 
  distinct(sampe_hla, .keep_all = T)   


TME_response_cor <- all_peptides %>%
  group_by(Sample) %>% 
  summarize(HLA_exp = mean(HLAexp, na.rm=TRUE)) %>% 
  left_join(.,TME_response) %>% 
  distinct(Sample, .keep_all = T) 

# Correlation_TME  <-TME_response_cor %>% ungroup() %>%   select(TME_cols) 
# Correlation_TME <- sapply(Correlation_TME,as.numeric)


# SupFig4C <- ggcorr(Correlation_TME , method = c("pairwise", "spearman"), 
#                    hjust = 0.8, size = 5, label = TRUE, label_size = 2) 
# ggsave(SupFig4C, file = "results/PaperPlots/Fig4/SupFig4C.pdf", width = 8, height = 8)


# all_peptides
# table(df$response)
# df <- all_peptides %>% 
#   dplyr::select(c(TME_cols,response,Partition))
# 
# AUC_df <- data.frame(matrix(ncol = 4, nrow = 0))
# 
# # 
# 
# df$response <- factor(df$response, levels = c(1,0))
# for (x in TME_cols) {
#   print(x)
#   for(i in unique(df$Partition)) {
#   #  print(i)
#     df_p <- df %>% filter(Partition!=i)
#     print(unique(df_p$Partition))
#     #table(df_p$response)
#    # auc(df_p$response, df_p[,x])
#     AUC_score  <- auc(df_p$response, df_p[,x])
#     pred <- prediction(df_p[,x],df_p$response)
#     auc_01 <- performance(pred,measure = "auc", fpr.stop=0.1)
#     AUC_df  <- rbind(AUC_df,c(x,i, round(AUC_score,4),round(auc_01@y.values[[1]],4)))
#   }
#   
# }
# 
# colnames(AUC_df) <- c("Feature","Partition","AUC","AUC_01")
# 
# AUC_df$AUC <- as.numeric(AUC_df$AUC)
# AUC_df$AUC_01 <- as.numeric(AUC_df$AUC_01)
# 
# 
# AUC_mean <- AUC_df %>% 
#   group_by(Feature) %>%
#   summarize(AUC = mean(AUC, na.rm=TRUE)) 
# 
# AUC_mean_01 <- AUC_df %>% 
#   group_by(Feature) %>%
#   summarize(AUC_01 = mean(AUC_01, na.rm=TRUE)) %>% 
#   mutate(Partition = "mean" ) %>% 
#   right_join(AUC_mean) %>% select(Feature,Partition,AUC,AUC_01)
# 
# 
# 
# AUC_df <- bind_rows(AUC_df,AUC_mean_01)
# 
# max(AUC_df$AUC)
#   AUC_heatmap_TME <- AUC_df %>% 
#     ggplot(., aes(x = Partition, y = reorder(Feature, -AUC))) +
#     geom_tile(aes(fill = AUC), colour = "white", size=0.25) +
#     scale_fill_distiller(name = "AUC", palette = "Spectral", na.value=c(aplha=0.1, aplha=0.1), 
#                          breaks = c(0.45,0.5,0.55,0.6), 
#                          labels =c("0.45","0.50","0.55","0.60")) +  
#     theme_bw() +
#     labs(x = "", y = "", fill = "AUC") + 
#     scale_x_discrete( expand = c(0, 0), drop=FALSE, 
#                       labels = c("Train 1","Train 2","Train 3","Train 4","Train 5","Mean"))+
#     scale_y_discrete(expand = c(0, 0),labels=Abbreviations)+ 
#     theme(axis.text.x = element_text(size=14),
#     #  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size=14),
#           axis.text.y = element_text(hjust = 0, size=14),
#           legend.title = element_text(size=16), 
#           axis.ticks = element_blank(),
#           panel.spacing = unit(0.1, "lines"),
#           axis.title.y = element_text(size=20, vjust = 1.5),  #, vjust = 1.03, hjust = -5, angle = 0
#           panel.grid.minor = element_line(colour="grey", size=12, linetype=1),
#           strip.text.x = element_text(size = 16), 
#           legend.spacing = grid::unit(0,"cm"),
#           legend.text = element_text(size = 12),
#           legend.position = "right") 
#   ggsave(AUC_heatmap_TME, file = paste("results/PaperPlots/Fig4/AUC_heatmap_TME","all",".pdf", sep = ""), width = 7.5 , height = 5 )
#   

# for every train set 
  # -------------------------------

for (i in c(0,1,2,3,4)) {
AUC_heatmap_TME <- AUC_df %>% filter(Partition==i) %>% 
  ggplot(., aes(x = Partition, y = reorder(Feature, -AUC))) +
  geom_tile(aes(fill = AUC), colour = "white", size=0.25) +
  scale_fill_distiller(name = "AUC", palette = "Spectral", na.value=c(aplha=0.1, aplha=0.1), 
                       breaks = c(0.001,0.004,0.008,0.012,0.016), 
                       labels =c("0.001","0.004","0.008","0.012","0.016")) +  
  theme_bw() +
  labs(x = "", y = "", fill = "AUC") + 
  scale_x_discrete( expand = c(0, 0), drop=FALSE)+
  scale_y_discrete(expand = c(0, 0),labels=Abbreviations)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size=14),
        axis.text.y = element_text(hjust = 0, size=12),
        legend.title = element_text(size=14), 
        axis.ticks = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        axis.title.y = element_text(size=20, vjust = 1.5),  #, vjust = 1.03, hjust = -5, angle = 0
        panel.grid.minor = element_line(colour="grey", size=12, linetype=1),
        strip.text.x = element_text(size = 16), 
        legend.spacing = grid::unit(0,"cm"),
        legend.text = element_text(size = 12),
        legend.position = "right") 
ggsave(AUC_heatmap_TME, file = paste("results/PaperPlots/Fig4/AUC_heatmap_TME/AUC_heatmap_TME_",i,".pdf", sep = ""), width = 7 , height = 5 )

}



