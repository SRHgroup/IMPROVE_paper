#source("bin/R_script/04_plotting/04_0_Paper_plots_prep.R")
# special library for Fig1
library(janitor)
# --
load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")
# ------------------------------------------------------------------------------
#                 Figure 1 - main 
# ------------------------------------------------------------------------------
nrow(all_peptides)
table(all_peptides$Patient)
all_peptides$response <- as.factor(all_peptides$response)

### investigations .. 
screened <- all_peptides %>% group_by(Patient) %>% tally(n = "screened") 
responses_tab <- all_peptides %>% filter(response==1) %>% group_by(Patient) %>% tally(n = "responses") 
screened <- screened %>% left_join(., responses_tab)


responses_tab <- all_peptides %>% filter(response==1) %>% group_by(Patient) %>% tally() 
tab <- all_peptides %>% group_by(Patient) %>% summarise(count = n_distinct(HLA_allele))

#                 Figure 1 - main 
# ------------------------------------------------------------------------------

Fig1C <- all_peptides %>% 
  group_by(cohort,Patient,response) %>% 
  tally() %>% 
  spread(., key = response, value = n) %>% 
  replace(is.na(.),0) %>% 
  mutate(total= `0`+`1`) %>% 
  mutate(fraction= `1`/total) %>%
  mutate(fraction_to_order= fraction) %>%
  mutate(Immunugenic = `1`) %>%
  gather(., key = type, value = number, -c(cohort,Patient,fraction_to_order)) %>% 
  filter(!type %in% c("0","1")) %>% 
  ggplot(. ) +
  geom_point(aes(x = reorder(Patient,-fraction_to_order) , y = number,color = type)) +
  #facet_grid(. ~ cohort, drop = T, space = "free", scales = "free") + 
  labs(y = "Neoepitopes", x = "", color = "Type") + 
  scale_y_log10(breaks = c(0.01,0.1,1,10,100,500), label = c("0.01","0.1","1","10","100","500"))+
  theme_bw() +
  scale_color_manual(breaks = c("total", "Immunugenic","fraction"), 
                     label = c("Total screened", "Total immunogenic","Fraction immunogenic") , 
                     values = c("black", "grey","red")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.title =element_text(size = 14) ) +
  guides(colour = guide_legend(override.aes = list(size=4)))


# Figure 1 all 
# ---------------------
pdf(file = 'results/PaperPlots/Fig1/Fig1_ALL.pdf', width = 14, height = 8)
ggdraw() +
  draw_plot(Fig1C, .45, .3, 0.55, 0.7) +
  draw_plot_label(c("A","B","C"), x = c(0,0,0.45), y = c(1.01,0.35,1.01), size = 22)

dev.off()



#                Supplementary Figure 1 - main 
# ------------------------------------------------------------------------------
# Supplementary Figure 1A
# ------------------------
SupFig1A_tab <- all_peptides %>% 
  group_by(cohort,response) %>% 
  tally() %>% 
  spread(., key = response, value = n) %>% 
  replace(is.na(.),0) %>% 
  mutate(total= `0`+`1`) %>% 
  mutate(yes= (`1`/total)*100) %>%
  mutate(no= (`0`/total)*100) %>%
  select(cohort,yes,no) %>% 
 # mutate(no= no-yes) %>%
  gather(., key = type, value = number, -cohort)

SupFig1A_tab$cohort <- factor(SupFig1A_tab$cohort, levels = c("Basket","Melanoma","mUC"))
SupFig1A_tab$type <- factor(SupFig1A_tab$type, levels = c("no","yes"))
SupFig1A <- SupFig1A_tab  %>% 
  ggplot(., aes(x = cohort, y = number,fill = type, alpha = type)) +
#  geom_bar(aes(fill = cohort, alpha = type ), stat = "identity") + 
  geom_col(position = "identity") +
  scale_alpha_manual(breaks = c("no","yes"), values = c(0.3,0.9), labels = c("No","Yes")) +
  scale_fill_manual(breaks = c("no","yes"), values = response_col, labels = c("No","Yes")  ) +
  theme_bw() + 
 # scale_fill_manual(breaks = c("Basket","Melanoma","mUC"), values = cohort_col) +
  #scale_y_log10()+
     scale_y_log10(breaks = c(0,1,2,3,5,10,20,50,100),
                   labels = c("0","1","2","3","5","10","20","50","100"),
                   limits = c(1,100)) +
  labs(x = "", y = "Percent \n neopetides", fill = "Immunogenic") + 
  theme(legend.position = "bottom",axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12) 
          )+
  guides(nrow = 2,alpha = "none") 

# Supplementary Figure 1B
# ------------------------------------------------------
SubFig1B <- all_peptides %>% 
  group_by(cohort,Patient) %>% 
  ggplot(., aes(x =RankEL)) +
  geom_density() + 
  theme_bw()+
  scale_x_log10(breaks = c( .01,.1, .5, 2,8,20), 
                minor_breaks = NULL, 
                labels = c( '0.01','0.1', '0.5', '2','8','20')) +
  labs(x ="RankEL", y = "Density" ) +
  theme(legend.position = "none",
        axis.title = element_text(size = 16))
# Supplementary Figure 1C
# --------------------------------------------------------------------
SubFig1C <- all_peptides %>% 
  group_by(cohort,Patient) %>% 
  ggplot(., aes(x =Expression)) +
  geom_density() + 
  theme_bw()+
  scale_x_log10(breaks = c( .01,.1, .5, 2,10, 100,1000), 
                minor_breaks = NULL, 
                labels = c( '0.01','0.1', '0.5','2' ,'10', '100','1000'),
                limit = c(0.01,3600)) +
  labs(x ="Expression", y = "Density" ) +
  theme(legend.position = "none",
        axis.title = element_text(size = 16))


# Supplementary Figure 1D
# ------------------------------------------------------
SupFig1D_test_type <- all_peptides %>%
  group_by(response,HLA_type) %>%
  tally() %>%
  spread(., key = response, value = n) 
  

SupFig1D_test_type[is.na(SupFig1D_test_type)] = 0
SupFig1D_test_type <- SupFig1D_test_type %>% mutate(total= `0`+`1`)
SupFig1D_test_type$yes <- as.numeric(SupFig1D_test_type$`1`)
SupFig1D_test_type[2:5] <- sapply(SupFig1D_test_type[2:5], as.numeric)


SupFig1D_test_type <- SupFig1D_test_type %>%
  adorn_totals("row")


HLA_test_type <- data.frame(matrix(ncol = 4, nrow = 0))
row_numb <- nrow(SupFig1D_test_type)-1
for (r in 1:row_numb){
#  p_val = 0
  print(r)
  pro1 = SupFig1D_test_type$yes[r]
  prop2 = SupFig1D_test_type$yes[4] - SupFig1D_test_type$yes[r]
  total1 = SupFig1D_test_type$total[r]
  total2 = SupFig1D_test_type$total[4] - SupFig1D_test_type$total[r]
  prop_test <- prop.test(c(pro1,prop2),c(total1,total2 ))
  p_val  = ifelse(prop_test$estimate[1]>prop_test$estimate[2], prop_test$p.value, "1")
  vec <- c(SupFig1D_test_type$HLA_type[r],p_val,pro1,total1)
  HLA_test_type <- rbind(HLA_test_type,vec)
  
} 
colnames(HLA_test_type) <- c("HLA_type","p-val","yes","total")


################## ALL HLA ALLELE ############################3

table(all_peptides$response)
SupFig1D_test <- all_peptides %>%
  group_by(HLA_allele,response,HLA_type) %>%
  tally() %>%
  spread(., key = response, value = n) 
nrow(SupFig1D_test )

SupFig1D_test[is.na(SupFig1D_test)] = 0
SupFig1D_test <- SupFig1D_test %>% mutate(total= `0`+`1`)
SupFig1D_test$yes <- as.numeric(SupFig1D_test$`1`)
SupFig1D_test[3:6] <- sapply(SupFig1D_test[3:6], as.numeric)

library(janitor)
SupFig1D_test <- SupFig1D_test %>%
  adorn_totals("row")

SupFig1D_test$yes <- as.numeric(SupFig1D_test$`1`)
SupFig1D_test$total <- as.numeric(SupFig1D_test$total)
prop.test(SupFig1D_test$yes,SupFig1D_test$total)

HLA_test <- data.frame(matrix(ncol = 5, nrow = 0))

nrow(SupFig1D_test)
for (r in 1:nrow(SupFig1D_test)){
  pro1 = SupFig1D_test$yes[r]
  prop2 = SupFig1D_test$yes[37] - SupFig1D_test$yes[r]
  
  total1 = SupFig1D_test$total[r]
  total2 = SupFig1D_test$total[37] - SupFig1D_test$total[r]
  prop_test <- print(prop.test(c(pro1,prop2),c(total1,total2)))
  p_val  = ifelse(prop_test$estimate[1]>prop_test$estimate[2], prop_test$p.value, "1")
  vec <- c(SupFig1D_test$HLA_allele[r],SupFig1D_test$HLA_type[r],p_val,pro1,total1)
  HLA_test <- rbind(HLA_test,vec)
  
} 
colnames(HLA_test) <- c("HLA_allele","HLA_type","p-val","yes","total")



HLA_test$yes <- as.numeric(HLA_test$yes)
HLA_test$total <- as.numeric(HLA_test$total)
HLA_test <- HLA_test %>% 
  mutate(fraction = yes/total) %>%
  mutate(percent = fraction*100) 
HLA_test$`p-val` <- as.numeric(HLA_test$`p-val`)
HLA_test$P_val_round <- round(HLA_test$`p-val`,2)
HLA_test <- HLA_test %>% mutate(p_val_sign = case_when(`p-val` > 0.05 ~ "NS.",
                                                       `p-val` < 0.05 & `p-val` > 0.01 ~ "*",
                                                       `p-val` < 0.01 & `p-val` > 0.001 ~ "**",
                                                       `p-val` < 0.001 ~ "***"))


HLA_test <- HLA_test %>% mutate(HLA_type_p = case_when(HLA_type=="HLA-C" ~ paste("*",HLA_test$HLA_type, sep = "\n"),
                                           TRUE ~ paste("NS.",HLA_test$HLA_type, sep = "\n")))

HLA_test$HLA_type_p <- factor(HLA_test$HLA_type_p, levels = c("NS.\nHLA-A","NS.\nHLA-B","*\nHLA-C"))
SupFig1D <- HLA_test %>%
  ggplot(. , aes(x = reorder(HLA_allele, -percent) , y = percent)) +
  geom_bar(aes(color = HLA_type_p, fill = HLA_type_p ), stat = "identity",colour="black")  +
  #scale_color_manual(values = HLA_col) +
  scale_fill_manual(values = HLA_col)  +
                    # breaks = c("NS.\nHLA-A","NS.\nHLA-B.","NS.\nHLA-C"), 
                    # labels = c("NS.\nHLA-A","NS.\nHLA-B","NS.\nHLA-C")) +
  labs(y = "Percent immunogenic", fill = "HLA type", color = "HLA type", x = "") +
  theme_bw() +
  geom_text(aes(label = p_val_sign), vjust = -0.2) +
 # annotate("text", x=5, y=0, label= "NS.", vjust = 2) + 
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=9))

all_peptides %>% filter(HLA_allele=="HLA-C06:02") %>% group_by(Patient,response) %>% tally()



### figure 1 gather 
pdf(file = 'results/PaperPlots/Fig1/SupFigure1_all.pdf', width = 10, height = 8)
ggdraw() +
  draw_plot(SupFig1A , .0, .68, .4, 0.32) +
  draw_plot(SubFig1B , .4, .73, .3, 0.27) +
  draw_plot(SubFig1C, .7, .73, .3, 0.27) +
  draw_plot(SupFig1D, .0, .0, 1, 0.7) +
  draw_plot_label(c("A","B","C","D"), x = c(-.01,.4,.7,-.01), y = c(1.01,1.01,1.01,0.72), size = 22)

dev.off()

