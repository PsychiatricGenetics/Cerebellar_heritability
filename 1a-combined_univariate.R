# Lachlan T Strike
# Direct estimates of variance components (ACE)
# Based on scripts shared by Hermine Maes (https://hermine-maes.squarespace.com)

rm(list = ls())


library(OpenMx)
library(tidyverse)
source("miFunctions.R")
source("1a-combined_univariate_function.R")
setwd("C:/GitHub/Cerebellar_heritability")
my.data <- readRDS("data/combined_cb_familywise.RDS")

# Single phenotype
result.single <- as_tibble(Univariate_ACDE(phenotype = "Right_VIZ", twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ", "cohort"), CorD = "C"))
result.single %>% select(Variable, AIC_BestFit, ACE_A_95CI, ACE_C_95CI, ACE_E_95CI, ACE_VA, ACE_VC, ACE_VE, AE_A_95CI, AE_E_95CI, AE_A, AE_E)

# Phenotype array
phenotype.array <- c("Left_I_III", "Right_I_III",
                     "Left_IV", "Right_IV",
                     "Left_V", "Right_V",
                     "Left_VI", "Right_VI",
                     "Left_Crus_I", "Right_Crus_I",
                     "Left_Crus_II", "Right_Crus_II",
                     "Left_VIIB", "Right_VIIB",
                     "Left_VIIIA", "Right_VIIIA",
                     "Left_VIIIB", "Right_VIIIB",
                     "Left_IX", "Right_IX",
                     "Left_X",  "Right_X",
                     "Vermis_VI", "Vermis_VII", "Vermis_VIII", "Vermis_IX", "Vermis_X",
                     "Corpus_Medullare", "Total_Cerebel_Vol")
phenotype.arrayZ <- paste0(phenotype.array, "Z")
result.multiple.ACE <- lapply(phenotype.arrayZ, Univariate_ACDE, twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ", "cohort"), CorD = "C") %>% bind_rows()
result.multiple.ACE <- as_tibble(result.multiple)

phenotype.array <- c("Left_I_III", "Right_I_III",
                     "Left_IV", "Right_IV",
                     "Left_V", "Right_V",
                     "Left_VI", "Right_VI",
                     "Left_Crus_I", "Right_Crus_I",
                     "Left_Crus_II", "Right_Crus_II",
                     "Left_VIIB", "Right_VIIB",
                     "Left_VIIIA", "Right_VIIIA",
                     "Left_VIIIB", "Right_VIIIB",
                     "Left_IX", "Right_IX",
                     "Left_X",  "Right_X",
                     "Vermis_VI", "Vermis_VII", "Vermis_VIII", "Vermis_IX", "Vermis_X",
                     "Corpus_Medullare", "Total_Cerebel_Vol")
phenotype.arrayZ <- paste0(phenotype.array, "Z")
result.multiple.ADE <- lapply(phenotype.arrayZ, Univariate_ACDE, twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ", "cohort"), CorD = "D") %>% bind_rows()
result.multiple.ADE <- as_tibble(result.multiple)

write.csv(result.multiple.ADE, "1a-combined_Univariate_ADE_output.csv", row.names = F)
