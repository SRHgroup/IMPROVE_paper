# find if mu tations are validated in RNA 

library(openxlsx)
library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
mutaion_val <- read.csv(file = "data/02_feature_data/Validated_RNA_seq/confirm_rnaseq_all_cohorts.csv", sep = " ")
load("data/04_plotting/Rdata/Final_prep_plotting.Rdata")
source("bin/R_script/99_functions.R")


mutaion_val$Sample <- sprintf("%02d", as.numeric(mutaion_val$Sample))
mutaion_val <- mutaion_val %>% 
  mutate(Patient = 
           case_when(Cohort=="Bladder" ~ paste("mUC",Sample, sep = "-"),
                     Cohort=="Melanoma" ~ paste("MM",Sample, sep = "-"),
                     Cohort=="Basket" ~ paste("RH",Sample, sep = "-"))) %>% 
  mutate(cohort = case_when(Cohort=="Melanoma" ~ "Melanoma", 
                            Cohort=="Bladder" ~ "mUC",
                            Cohort=="Basket" ~ "Basket")) 

mutaion_val <- mutaion_val %>% 
  mutate(identi_pep_patient = paste(Mut_peptide,HLA_allele,Patient, sep = "_")) %>% 
  distinct(identi_pep_patient, .keep_all = T)

colnames(all_peptides)
unique(all_peptides$cohort)
all_peptides_with_rnamut_orginal <- all_peptides %>% 
  mutate(identi_pep_patient = paste(Mut_peptide,HLA_allele,Patient, sep = "_")) %>% 
  distinct(identi_pep_patient, .keep_all = T) %>% 
  select(identi_pep_patient, Patient,response) 

all_peptides_with_rnamut_orginal <- all_peptides_with_rnamut_orginal %>% 
  left_join(.,mutaion_val ) %>% 
  mutate(response_lab = case_when(response=="1" ~ "yes",
                                  TRUE ~ "no"))
table(all_peptides_with_rnamut_orginal$cohort)


all_peptides_with_rnamut_orginal$rna_bin[all_peptides_with_rnamut_orginal$rna_confirm=="0/0;"] <- "no_coverage"


all_peptides <- all_peptides %>% left_join(.,mutaion_val)


all_peptides_with_rnamut_orginal$rna_bin[is.na(all_peptides_with_rnamut_orginal$rna_bin==T)] <- "Not validated"

#fig2B
Fig2C_tab <- all_peptides_with_rnamut_orginal %>% 
  arrange(response) %>% 
  # distinct(Mut_identity, .keep_all = T) %>% 
  mutate(response = as.factor(response)) %>% 
  group_by(rna_bin,response)  %>% 
  tally()  %>% 
  spread(key = "response", value = "n")

setDT(Fig2C_tab)[, frac_neg := (`0` / sum(`0`))*100]
setDT(Fig2C_tab)[, frac_pos := (`1` / sum(`1`))*100]

Fig2C_tab <- Fig2C_tab %>% as.tibble() %>% 
  dplyr::select(rna_bin,frac_neg,frac_pos) %>% 
  gather(., key = "type", value = "fraction" , -rna_bin) %>% 
  mutate(p_val = case_when(rna_bin == "0" & type == "frac_neg" ~ "p = 0.47",
                           rna_bin == "1" & type == "frac_neg" ~ "p = 0.47",
                           rna_bin == "no_coverage" & type == "frac_neg" ~ "p = 0.22",
                           rna_bin == "Not validated" & type == "frac_neg" ~ "p = 0.50"))

all_peptides_with_rnamut_orginal_stats <- 
  all_peptides_with_rnamut_orginal  %>% 
  filter(rna_bin != "Not validated")
all_peptides_with_rnamut_orginal_stats$good <- ifelse(all_peptides_with_rnamut_orginal_stats$rna_bin == "1",  "good", "bad")
prop_test(all_peptides_with_rnamut_orginal_stats$good,all_peptides_with_rnamut_orginal_stats$response,"good","response")
table(all_peptides_with_rnamut_orginal_stats$response,all_peptides_with_rnamut_orginal_stats$good)

all_peptides_with_rnamut_orginal_stats$good <- ifelse(all_peptides_with_rnamut_orginal_stats$rna_bin == "0",  "good", "bad")
prop_test(all_peptides_with_rnamut_orginal_stats$good,all_peptides_with_rnamut_orginal_stats$response,"good","response")
table(all_peptides_with_rnamut_orginal_stats$response,all_peptides_with_rnamut_orginal_stats$good)

all_peptides_with_rnamut_orginal_stats$good <- ifelse(all_peptides_with_rnamut_orginal_stats$rna_bin == "Not validated",  "good", "bad")
prop_test(all_peptides_with_rnamut_orginal_stats$good,all_peptides_with_rnamut_orginal_stats$response,"good","response")
table(all_peptides_with_rnamut_orginal_stats$response,all_peptides_with_rnamut_orginal_stats$good)

all_peptides_with_rnamut_orginal_stats$good <- ifelse(all_peptides_with_rnamut_orginal_stats$rna_bin == "no_coverage",  "good", "bad")
prop_test(all_peptides_with_rnamut_orginal_stats$good,all_peptides_with_rnamut_orginal_stats$response,"good","response")
table(all_peptides_with_rnamut_orginal_stats$response,all_peptides_with_rnamut_orginal_stats$good)



Fig2E_RNA <- Fig2C_tab  %>% filter(rna_bin %in% c(1,0,"no_coverage")) %>% 
  ggplot(aes(x=rna_bin, y = fraction, fill = type))+
  geom_bar(aes(fill = type), stat = "identity",, position = "dodge")+
  scale_fill_manual(breaks = c("frac_neg","frac_pos"), 
                    values = c("#91bfdb", "#ef8a62") , 
                    labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("1","0","no_coverage"),
                   labels = c("Found in \n RNA","Not found\n in RNA","No sufficent\n coverage"))+
  theme_bw()+
  scale_y_continuous(limits = c(0,60))+
  # geom_text_repel(aes(label = round(fraction,2)), position = position_stack(vjust = 0.9),direction = "y", 
  #           box.padding = unit(0.01, "lines"))+
  geom_text(aes(label = p_val), vjust = -1.8, size = 3.5) +
  labs(x = "", y = "Percent neopeptides", fill = "Immunogenic") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, size = 14),
        strip.background = element_blank(),
        axis.title = element_text(size = 18),
        strip.text = element_text(size=9))
#ggsave(Fig2E_RNA, file="results/PaperPlots/Fig2/Fig2E_RNA.pdf", height = 4, width = 4 )

table(all_peptides$response_lab,all_peptides$rna_bin,all_peptides$cohort)
table(all_peptides_with_rnamut_orginal$response)
table(all_peptides_with_rnamut_orginal$rna_bin,all_peptides_with_rnamut_orginal$response)


Fig2C_tab_basket <- all_peptides_with_rnamut_orginal %>% 
  filter(cohort=="Basket") 

Fig2C_tab_basket$good <- ifelse(Fig2C_tab_basket$rna_bin == "1",  "good", "bad")
prop_test(Fig2C_tab_basket$good,Fig2C_tab_basket$response,"good","response")
table(Fig2C_tab_basket$response,Fig2C_tab_basket$good)

Fig2C_tab_basket$good <- ifelse(Fig2C_tab_basket$rna_bin == "0",  "good", "bad")
prop_test(Fig2C_tab_basket$good,Fig2C_tab_basket$response,"good","response")
table(Fig2C_tab_basket$response,Fig2C_tab_basket$good)

Fig2C_tab_basket$good <- ifelse(Fig2C_tab_basket$rna_bin == "Not validated",  "good", "bad")
prop_test(Fig2C_tab_basket$good,Fig2C_tab_basket$response,"good","response")
table(Fig2C_tab_basket$response,Fig2C_tab_basket$good)

Fig2C_tab_basket$good <- ifelse(Fig2C_tab_basket$rna_bin == "no_coverage",  "good", "bad")
prop_test(Fig2C_tab_basket$good,Fig2C_tab_basket$response,"good","response")
table(Fig2C_tab_basket$response,Fig2C_tab_basket$good)



#fig2B
Fig2C_tab_basket <- Fig2C_tab_basket %>% 
  arrange(response) %>% 
  # distinct(Mut_identity, .keep_all = T) %>% 
  mutate(response = as.factor(response)) %>% 
  group_by(rna_bin,response)  %>% 
  tally()  %>% 
  spread(key = "response", value = "n")

setDT(Fig2C_tab_basket)[, frac_neg := (`0` / sum(`0`))*100]
setDT(Fig2C_tab_basket)[, frac_pos := (`1` / sum(`1`))*100]

Fig2C_tab_basket <- Fig2C_tab_basket %>% as.tibble() %>% 
  dplyr::select(rna_bin,frac_neg,frac_pos) %>% 
  gather(., key = "type", value = "fraction" , -rna_bin) %>% 
  mutate(p_val = case_when(rna_bin == "0" & type == "frac_neg" ~ "p = 0.83",
                           rna_bin == "1" & type == "frac_neg" ~ "p = 0.15",
                           rna_bin == "no_coverage" & type == "frac_neg" ~ "p = 1",
                           rna_bin == "Not validated" & type == "frac_neg" ~ "p = 0.13"))
Fig2C_tab_basket$cohort <- "Basket"

Fig2C_tab_Bladder <- all_peptides_with_rnamut_orginal %>% 
  filter(cohort=="mUC") 
  
  
Fig2C_tab_Bladder$good <- ifelse(Fig2C_tab_Bladder$rna_bin == "1",  "good", "bad")
prop_test(Fig2C_tab_Bladder$good,Fig2C_tab_Bladder$response,"good","response")
table(Fig2C_tab_Bladder$response,Fig2C_tab_Bladder$good)

Fig2C_tab_Bladder$good <- ifelse(Fig2C_tab_Bladder$rna_bin == "0",  "good", "bad")
prop_test(Fig2C_tab_Bladder$good,Fig2C_tab_Bladder$response,"good","response")
table(Fig2C_tab_Bladder$response,Fig2C_tab_Bladder$good)

Fig2C_tab_Bladder$good <- ifelse(Fig2C_tab_Bladder$rna_bin == "Not validated",  "good", "bad")
prop_test(Fig2C_tab_Bladder$good,Fig2C_tab_Bladder$response,"good","response")
table(Fig2C_tab_Bladder$response,Fig2C_tab_Bladder$good)

Fig2C_tab_Bladder$good <- ifelse(Fig2C_tab_Bladder$rna_bin == "no_coverage",  "good", "bad")
prop_test(Fig2C_tab_Bladder$good,Fig2C_tab_Bladder$response,"good","response")
table(Fig2C_tab_Bladder$response,Fig2C_tab_Bladder$good)
  
Fig2C_tab_Bladder <- Fig2C_tab_Bladder %>% 
  arrange(response) %>% 
  # distinct(Mut_identity, .keep_all = T) %>% 
  mutate(response = as.factor(response)) %>% 
  group_by(rna_bin,response)  %>% 
  tally()  %>% 
  spread(key = "response", value = "n")

setDT(Fig2C_tab_Bladder)[, frac_neg := (`0` / sum(`0`))*100]
setDT(Fig2C_tab_Bladder)[, frac_pos := (`1` / sum(`1`))*100]

Fig2C_tab_Bladder <- Fig2C_tab_Bladder %>% as.tibble() %>% 
  dplyr::select(rna_bin,frac_neg,frac_pos) %>% 
  gather(., key = "type", value = "fraction" , -rna_bin) %>% 
  mutate(p_val = case_when(rna_bin == "0" & type == "frac_neg" ~ "p = 0.80",
                           rna_bin == "1" & type == "frac_neg" ~ "p = 0.81",
                           rna_bin == "no_coverage" & type == "frac_neg" ~ "p = 0.55",
                           rna_bin == "Not validated" & type == "frac_neg" ~ "p = 0.64"))
Fig2C_tab_Bladder$cohort <- "mUC"

Fig2C_tab_melanoma <- all_peptides_with_rnamut_orginal %>% 
  filter(cohort=="Melanoma") 
  
Fig2C_tab_melanoma$good <- ifelse(Fig2C_tab_melanoma$rna_bin == "1",  "good", "bad")
prop_test(Fig2C_tab_melanoma$good,Fig2C_tab_melanoma$response,"good","response")
table(Fig2C_tab_melanoma$response,Fig2C_tab_melanoma$good)

Fig2C_tab_melanoma$good <- ifelse(Fig2C_tab_melanoma$rna_bin == "0",  "good", "bad")
prop_test(Fig2C_tab_melanoma$good,Fig2C_tab_melanoma$response,"good","response")
table(Fig2C_tab_melanoma$response,Fig2C_tab_melanoma$good)

Fig2C_tab_melanoma$good <- ifelse(Fig2C_tab_melanoma$rna_bin == "no_coverage",  "good", "bad")
prop_test(Fig2C_tab_melanoma$good,Fig2C_tab_melanoma$response,"good","response")
table(Fig2C_tab_melanoma$response,Fig2C_tab_melanoma$good)  

Fig2C_tab_melanoma$good <- ifelse(Fig2C_tab_melanoma$rna_bin == "Not validated",  "good", "bad")
prop_test(Fig2C_tab_melanoma$good,Fig2C_tab_melanoma$response,"good","response")
table(Fig2C_tab_melanoma$response,Fig2C_tab_melanoma$good)  

Fig2C_tab_melanoma <- Fig2C_tab_melanoma %>% 
  arrange(response) %>% 
  # distinct(Mut_identity, .keep_all = T) %>% 
  mutate(response = as.factor(response)) %>% 
  group_by(rna_bin,response)  %>% 
  tally()  %>% 
  spread(key = "response", value = "n")

setDT(Fig2C_tab_melanoma)[, frac_neg := (`0` / sum(`0`))*100]
setDT(Fig2C_tab_melanoma)[, frac_pos := (`1` / sum(`1`))*100]

Fig2C_tab_Melanoma <- Fig2C_tab_melanoma %>% as.tibble() %>% 
  dplyr::select(rna_bin,frac_neg,frac_pos) %>% 
  gather(., key = "type", value = "fraction" , -rna_bin) %>% 
  mutate(p_val = case_when(rna_bin == "0" & type == "frac_neg" ~ "p = 0.31",
                           rna_bin == "1" & type == "frac_neg" ~ "p = 0.84",
                           rna_bin == "no_coverage" & type == "frac_neg" ~ "p = 0.24",
                           rna_bin == "Not validated" & type == "frac_neg" ~ "p = 0.33"))

Fig2C_tab_Melanoma$cohort <- "Melanoma"


Fig2C_tab_per_cohort <- bind_rows(Fig2C_tab_Melanoma,Fig2C_tab_Bladder,Fig2C_tab_basket)


SupFig2C_RNA_cohort <- Fig2C_tab_per_cohort  %>% filter(rna_bin %in% c(1,0,"no_coverage")) %>% 
  ggplot(aes(x=rna_bin, y = fraction, fill = type))+
  geom_bar(aes(fill = type), stat = "identity",, position = "dodge")+
  scale_fill_manual(breaks = c("frac_neg","frac_pos"), 
                    values = c("#91bfdb", "#ef8a62") , 
                    labels = c("No","Yes")) +
  scale_x_discrete(breaks = c("1","0","no_coverage"),
                   labels = c("Found\n in RNA","Not found \n in RNA","No sufficent\n coverage"))+
  theme_bw()+
  facet_grid(.~cohort) +
  scale_y_continuous(limits = c(0,60))+
  # geom_text_repel(aes(label = round(fraction,2)), position = position_stack(vjust = 0.9),direction = "y", 
  #           box.padding = unit(0.01, "lines"))+
  geom_text(aes(label = p_val), vjust = -3, size = 3.5) +
  labs(x = "", y = "Percent neopeptides", fill = "Immunogenic") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 11),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
      #  strip.background = element_blank(),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size=16))
#ggsave(SupFig2C_RNA_cohort, file="results/PaperPlots/Fig2/SupFig2C_RNA_cohort.pdf", height = 5, width = 8 )


save(SupFig2C_RNA_cohort ,Fig2E_RNA, file = "data/04_plotting/RNA_validation/RNA_val_figuress.Rdata")

