# Lachlan T Strike
# Direct estimates of variance components
# Based on scripts shared by Hermine Maes (https://hermine-maes.squarespace.com)

rm(list = ls())

library(OpenMx)
library(tidyverse)
setwd("C:/GitHub/Cerebellar_heritability")
source("miFunctions.R")
source("4-qtim_univariate_function.R")
my.data <- readRDS("data/qtim_cb_familywise.RDS")

# Single phenotype
result.single <- as_tibble(Univariate_ACDE(phenotype = "Right_VIIBZ", twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ"), CorD = "C"))
result.single %>% select(Variable, ACE_A_95CI, ACE_C_95CI, ACE_E_95CI, AE_A_95CI, AE_E_95CI, AE_A, AE_E)

#### ACE ####
variable.list <- c("Left_I_III", "Right_I_III",
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
variable.listZ <- paste0(variable.list, "Z")
result.multiple.ACE <- lapply(variable.listZ, Univariate_ACDE, twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ"), CorD = "C") %>% bind_rows()
write.csv(result.multiple.ACE, "4-qtim_univariate_ACE_output.csv", row.names = F)

#### ADE ####
result.multiple.ADE <- lapply(variable.listZ, Univariate_ACDE, twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ"), CorD = "D") %>% bind_rows()
# Tidy up results to match ADE model 
result.multiple.ADE <- 
  result.multiple.ADE %>% select(-CE_C_95CI, -CE_E_95CI, -CE_VA, -CE_VE, -CE_AIC)
result.multiple.ADE <- 
  result.multiple.ADE %>% rename("ADE_A_95CI"="ACE_A_95CI",
                                 "ADE_D_95CI"="ACE_C_95CI",
                                 "ADE_E_95CI"="ACE_E_95CI",
                                 "ADE_VA"="ACE_VA",
                                 "ADE_VD"="ACE_VC",
                                 "ADE_VE"="ACE_VE",
                                 "ADE_AIC"="ACE_AIC",
                                 "dropD"="dropC",
                                 "dropAD"="dropAC")
write.csv(result.multiple.ADE, "4-qtim_univariate_ADE_output.csv", row.names = F)
