# Lachlan T Strike
# Direct estimates of variance components
# Based on scripts shared by Hermine Maes (https://hermine-maes.squarespace.com)
# Two twins

rm(list=ls())

setwd("C:/GitHub/Cerebellar_heritability")
library(OpenMx)
library(tidyverse)
source("miFunctions.R")
source("8-qtab_bivariate_AE_function.R")

my.data <- readRDS("data/qtab_cb_familywise.RDS")
Bivariate_AE(phenotypes = c("Left_Crus_IZ", "Corpus_MedullareZ"), twin.data = my.data)

#### Run all pairwise correlations ####
list1 <- c("Left_I_IIIZ", "Right_I_IIIZ", "Left_IVZ", "Right_IVZ", "Left_VZ", 
           "Right_VZ", "Left_VIZ",  "Right_VIZ", "Left_Crus_IZ", "Right_Crus_IZ",
           "Left_Crus_IIZ", "Right_Crus_IIZ", "Left_VIIBZ", "Right_VIIBZ", "Left_VIIIAZ", 
           "Right_VIIIAZ", "Left_VIIIBZ", "Right_VIIIBZ", "Left_IXZ", "Right_IXZ", 
           "Left_XZ", "Right_XZ", "Vermis_VIZ", "Vermis_VIIZ", "Vermis_VIIIZ",
           "Vermis_IXZ", "Vermis_XZ", "Corpus_MedullareZ", "Total_Cerebel_VolZ")

for (i in list1){
  list2 <- list1[!list1==i]
  cb.vars1 <- crossing(i, list2)
  cb.vars2 <- as.list(as.data.frame(t(cb.vars1)))
  results <- lapply(cb.vars2, Bivariate_AE, twin.data = my.data) %>% bind_rows()
  write.csv(results, paste0("sfigure3/8-qtab_bivariate_AE_", i, "_output.csv"), row.names = F)
}