##### mutation possistions 

library(janitor)
# after running the script 
# -------------------------------------------------------------------------------
mutation_neg <- read.table(file = "data/04_plotting/Mutation_possistion/dfMaster_1mut_core_neg_frac.csv", sep = "\t", header = TRUE)

mutation_pos <- read.table(file = "data/04_plotting/Mutation_possistion/dfMaster_1mut_core_pos_frac.csv", sep = "\t", header = TRUE)

mutation_type <- read.table(file = "data/04_plotting/Mutation_possistion/dfMaster_1mut_core_MutType_frac.csv", sep = "\t", header = TRUE)


mutation_raw <- read.table(file = "data/04_plotting/Mutation_possistion/dfMaster_1mut_core.csv", sep = "\t", header = TRUE)

###### make prop test mutpos per pept len  #####
table_pos <- mutation_raw %>% filter(Target==1) %>% group_by(PeptLen, CoreMutPos) %>% tally(n = "count_pos")
total_pos <- mutation_raw %>% filter(Target==1) %>% group_by(PeptLen) %>% tally(n = "total_pos")
table_pos <- table_pos %>% left_join(., total_pos)

table_neg <- mutation_raw %>% filter(Target==0) %>% group_by(PeptLen, CoreMutPos) %>% tally(n = "count_neg")
total_neg <- mutation_raw %>% filter(Target==0) %>% group_by(PeptLen) %>% tally(n = "total_neg")
table_neg <- table_neg %>% left_join(., total_neg)

table_all <- table_neg %>% left_join(.,table_pos)
table_all[is.na(table_all)] = 0

table_all_with_zero <- table_all %>% filter(count_pos==0)
table_all <- table_all %>% filter(count_pos>0)
Mut_pos_stats <- data.frame(matrix(ncol = 3, nrow = 0))

for (r in 1:nrow(table_all)){
  print(r)
  prop1 = table_all$count_neg[r]
  prop2 = table_all$count_pos[r]
  
  total1 = table_all$total_neg[r]
  total2 = table_all$total_pos[r]
  
  prop_test <- print(prop.test(c(prop1,prop2),c(total1,total2 )))
  p_val <- prop_test$p.value
  vec <- c(table_all$PeptLen[r], table_all$CoreMutPos[r],p_val)
  Mut_pos_stats<- rbind(Mut_pos_stats,vec)
} 
colnames(Mut_pos_stats) <- c("PeptLen","CoreMutPos","p-val")

# ----------------------------------------
########## mutpos all peptlen 
# ----------------------------------------------

Mut_pos_core <- mutation_raw %>% #filter(pep_length==10) %>% 
  group_by(CoreMutPos,Target) %>% filter(pep_length==10) %>% 
  tally() %>% 
  spread(., key = "Target", value = "n") %>% 
  adorn_totals("row")
colnames(Mut_pos_core) <- c("CoreMutPos", "no","yes")
Mut_pos_core$yes[is.na(Mut_pos_core$yes)==T] <- 0
Mut_pos_core_stats <- data.frame(matrix(ncol = 4, nrow = 0))

for (r in 1:nrow(Mut_pos_core)){
  print(Mut_pos_core$CoreMutPos[r])
  prop1 = Mut_pos_core$yes[r]
  prop2 = Mut_pos_core$no[r]
  
  total1 = Mut_pos_core$yes[11] - prop1
  total2 = Mut_pos_core$no[11] - prop2
  
  prop_test <- print(prop.test(c(prop1,prop2),c(total1,total2 )))
  p_val <- prop_test$p.value
  vec <- c(Mut_pos_core$CoreMutPos[r], Mut_pos_core$no[r],Mut_pos_core$yes[r], p_val)
  print(vec)
  Mut_pos_core_stats <- rbind(Mut_pos_core_stats,vec)
  
  
} 
colnames(Mut_pos_core_stats) <- c("CoreMutPos","no","yes","p-val")

prop.test(c(11,312),c(386,13646))

1662/(14953-1662)
# -----------------------------------------------


colnames(mutation_raw$Target)
mutation_pos$response <- "Immunogenic"
mutation_neg$response <- "Non-immunogenic"
mutation <- bind_rows(mutation_pos,mutation_neg)

mutation$PeptLen <- as.factor(mutation$PeptLen)


mutation <- mutation %>% 
  mutate(Core_type = case_when(CoreMutPos==0 ~ "OC" ,
                               TRUE ~ "Core")) 
mutation$Core_type <- factor(mutation$Core_type, levels = c("OC","Core"))
#mutation$CoreMutPos[mutation$CoreMutPos==0] <- "gap"

mut_pos <- ggplot(mutation, aes(x =CoreMutPos, y=percent )) +
  geom_line(aes(color = PeptLen)) +
  geom_point(aes(color = PeptLen), size = 3) +
  facet_grid(response ~ Core_type , scales = "free_x", drop = T, space = "free_x") +
  theme_bw() +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9), 
                     labels = c("gap","1","2","3","4","5","6","7","8","9") ) + 
  theme_bw() + 
  labs(x = "Binding core mutation position", color = "Peptide length", y = "Percent neopeptides" ) +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14), 
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "bottom")  +
  guides(colour = guide_legend(override.aes = list(size=2)))
#ggsave(mut_pos, file="results/PaperPlots/Fig2/Fig2_mut_pos.pdf", height = 4, width = 6 )


gt = ggplot_gtable(ggplot_build(mut_pos))
gt$widths[5] = 8*gt$widths[5]
Fig2B_mut_pos <- as_ggplot(gt)
# pdf(file="results/PaperPlots/Fig2/Fig2_mut_pos.pdf", height = 5, width = 5)
# grid.draw(gt)
# dev.off()






mutation_type$PeptLen <- as.factor(mutation_type$PeptLen)

mutation_type <- mutation_type %>% mutate(length_mut_type = paste(PeptLen,MutType, sep = "_"))
unique(mutation_type$length_mut_type)


mutation_type <- mutation_type %>% 
  mutate(Core_type = case_when(CoreMutPos==0 ~ "OC" ,
                               TRUE ~ "Core")) 
mutation_type$Core_type <- factor(mutation_type$Core_type, levels = c("OC","Core"))

mut_type <- ggplot(mutation_type, aes(x =CoreMutPos, y=percent, 
                                      group = interaction(PeptLen,MutType),
                                      color = MutType)) +
  geom_line() +
  geom_point(aes(color= MutType), size = 3) +
  scale_color_manual(breaks = c("Conserved","Improved"), values = c("#2ca25f","#88419d"))+
  theme_bw() +
  scale_y_continuous(breaks = c(0,5,10,15,20,25))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9), 
                     labels = c("gap","1","2","3","4","5","6","7","8","9") ) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        strip.text.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14), 
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "bottom")  +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  facet_grid(.~ Core_type, scales = "free_x", drop = T, space = "free_x") +
  labs(x = "Binding core mutation position", color = "MHC binding", y = "Percent neopeptide" )
#ggsave(mut_type, file="results/PaperPlots/Fig2/Fig2_mut_type.pdf", height = 4, width = 6 )

CB_IB_legend <- get_legend(mut_type)
mut_type <- mut_type + theme(legend.position = "none")

gt = ggplot_gtable(ggplot_build(mut_type))
gt$widths[5] = 8*gt$widths[5]
Fig2C_mut_type <- as_ggplot(gt)
# pdf(file = 'results/PaperPlots/Fig2/Fig2_mut_type.pdf', height = 5, width = 5)
# grid.draw(gt)
# dev.off()

save(Fig2C_mut_type ,Fig2B_mut_pos, file = "data/04_plotting/Mutation_possistion/mut_figures.Rdata")



