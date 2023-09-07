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

# feature calclulations 
# -----------------------------------------------
# Run python feature calculation script first
# python3 bin/python_script/feature_calculations.py --file data/01_data/01_cedar_benchmark_neoepitopes.tsv --dataset "Validated_neoepitopes" --ProgramDir "/Users/annieborch/Documents/programs/" --outfile "data/02_feature_data/02_Benchmark_neoepitopes_calculated_features.tsv" --TmpDir "/Users/annieborch/Documents/programs/"
# -----------------------------------------------



benchmark_data <- read.csv("data/02_feature_data/02_Benchmark_neoepitopes_calculated_features.tsv", header = TRUE, sep = '\t')

colnames(benchmark_data)
benchmark_data
benchmark_data <- benchmark_data %>% 
  rename( "DAI" = "DAI_4.1") %>% 
  rename("RankEL" = "RankEL_4.1") %>% 
  rename("RankBA" = "RankBA_4.1")  %>% 
  rename("Foreigness" = "Foreignness")  %>% 
  mutate(Patient = "Unknown") %>% 
  filter(RankEL<2)


write.table(benchmark_data, file = "data/02_feature_data/02_2_Benchmark_neoepitopes_calculated_features.tsv", sep = "\t", quote = F, row.names=F)

# Predictions
# -----------------------------------------------
# Run python feature calculation script first
# python3 bin/python_script/Predict_immunogenicity.py --file "data/02_feature_data/02_2_Benchmark_neoepitopes_calculated_features.tsv" --model Simple  --outfile "benchmark_data/04_Benchmark_neoepitopes_prediction.tsv" 

# -----------------------------------------------

