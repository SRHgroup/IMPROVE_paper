# Paper plots 20221112
# -----------------------------
# Annie Borch
# -------------------------------------------------------------
####                     Plots for paper #######
# -------------------------------------------------------------
# plot in generel new 
# load libraries 
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggbeeswarm)
library(ggrepel)
library(cowplot)
library(ggpubr)
# for roc curves 
#load necessary packages
library(ggplot2)

# core plot 
#install.packages("devtools")
#devtools::install_github("taiyun/corrplot", build_vignettes = TRUE)
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


# load data 
source("bin/R_script/04_plotting/00_Paper_plots_prep.R")

# ------------------------------------------------------------------------------
#                 Figure 1 - main 
# ------------------------------------------------------------------------------

source("bin/R_script/04_plotting/01_Figure1.R")


# ------------------------------------------------------------------------------
#                 Figure 2
# ------------------------------------------------------------------------------
source("bin/R_script/04_plotting/02_Figure2.R")

# --------------------------------------------
#               Figure 3  and Figure 4 
# --------------------------------------------


source("bin/R_script/04_plotting/03_Figure3_and_4.R")



# ------------------------------------------------------------------------------
#                 Figure 5
# ------------------------------------------------------------------------------


source("bin/R_script/04_plotting/05_Figure5.R")


# ------------------------------------------------------------------------------
#                 Figure 5
# ------------------------------------------------------------------------------


source("bin/R_script/04_plotting/06_1_benchmark_data_prep.R")
source("bin/R_script/04_plotting/06_2_benchmark_results.R")

