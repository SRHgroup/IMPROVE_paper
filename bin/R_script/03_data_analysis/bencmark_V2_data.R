## thrid benchmark data 

library(tidyverse)
library(openxlsx)

bench_data <- read.xlsx("data/benchmark_comparison/100391_sup_1.xlsx")

bench_data <- bench_data[,-1]
colnames(bench_data) <- bench_data[1,]
bench_data <- bench_data[-1,]
bench_data <- bench_data %>% 
  rename("Mut_peptide" = `Mutant Peptide` ) %>% 
  rename("Norm_peptide" = `Wild Type Peptide` ) %>% 
  rename("HLA_allele" = MHC ) %>% 
  rename( "NetMHCExp" = `NetMHCpanExp rank` ) %>% 
  rename( "Expression" = `HPA expression`)# %>% 
 # rename( "RankEL" = `NetMHCpan 4.0`)
  
  
bench_data$Patient <- paste("Patient", bench_data$Patient)

# add forginess score 
# ---------------------

load("/Users/annieborch/Documents/IMPROVE/IMPROVE_paper/data/02_feature_data/foreignness_score/foreignness_score_output_becnhmark_data_V2.Rdata")

bench_data <- bench_data %>% left_join(., fs_EXTRA_data %>% select(-nmer))

bench_data <- bench_data %>% rename("Foreigness" = foreignness_score)

write.table(bench_data, file = "data/01_data/01_bench_data_V2_for_feature_cal.tsv", sep = "\t", quote = F)


# add Expression fro NetMHCexppan 
# ---------------------




# feature calclulations 
# -----------------------------------------------
# Run python feature calculation script first
# python3 bin/python_script/feature_calculations.py --file data/01_data/01_bench_data_V2_for_feature_cal.tsv --dataset "becnh_data_neoepitopes" --ProgramDir "/Users/annieborch/Documents/programs/" --outfile "data/02_feature_data/02_Benchmark_data_V2_neoepitopes_calculated_features.tsv" --TmpDir "/Users/annieborch/Documents/programs/"
# -----------------------------------------------

beanch_feature <- read.table("data/02_feature_data/02_Benchmark_data_V2_neoepitopes_calculated_features.tsv", header = T, sep = "\t")
beanch_feature <- beanch_feature %>% 
  rename("RankEL" = "RankEL_4.1")  %>% 
  rename("RankBA" = "RankBA_4.1")  %>% 
  rename("DAI" = "DAI_4.1")  

write.table(beanch_feature, file = "data/02_feature_data/02_2_Benchmark_data_V2_neoepitopes_calculated_features.tsv", sep = "\t", quote = F )
# -----------------------------------------------
# Run python prediction script
# python3 bin/python_script/Predict_immunogenicity.py --file "data/02_feature_data/02_2_Benchmark_data_V2_neoepitopes_calculated_features.tsv" --model Simple  --outfile "benchmark_data/04_Benchmark_data_V2_neoepitopes_prediction.tsv" 

# -----------------------------------------------
