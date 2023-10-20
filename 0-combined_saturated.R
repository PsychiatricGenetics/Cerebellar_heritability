# Lachlan T Strike
# Saturated model used to check assumptions of twin data. Is setup for five zygosity groups, with opposite sex twin pairs coded M = 1, F = 0.
# Take rMZ/DZ & covariate estimates from model with equated means/variances, but separate covariances for MZ and DZ (i.e. CModel2)
# Sex - 1 = male, 0 = female
# Expects phenotype data to be Z-score, with "Z" appended to the end of the variable name (if not, comment out 'Get raw means & sd for table' section)
# Includes one singleton sibling of twins
# Based on scripts shared by Hermine Maes (https://hermine-maes.squarespace.com)

rm(list = ls())

setwd("C:/GitHub/Cerebellar_heritability")
library(OpenMx)
library(stringr)
library(dplyr)
source("miFunctions.R")
source("0-combined_saturated_function.R")

twin.data <- readRDS("data/combined_cb_familywise.RDS")

# Run for single phenotype
Saturated_covariate(phenotype = "Right_Crus_IIZ", twin.data = twin.data)

# Run for list of phenotypes
variable_list <- c("Left_I_III", "Right_I_III",
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
variable_list <- paste0(variable_list, "Z")
results.sat <- as_tibble(lapply(variable_list, Saturated_covariate, twin.data = twin.data) %>% bind_rows())

# Multiple Testing Correction
results.sat$H1m_pval_fdr <- p.adjust(p = results.sat$H1m_pval, method = "fdr")
results.sat$H2m_pval_fdr <- p.adjust(p = results.sat$H2m_pval, method = "fdr")
results.sat$H3m_pval_fdr <- p.adjust(p = results.sat$H3m_pval, method = "fdr")
results.sat$H4m_pval_fdr <- p.adjust(p = results.sat$H4m_pval, method = "fdr")
results.sat$H5m_pval_fdr <- p.adjust(p = results.sat$H5m_pval, method = "fdr")
results.sat$H1v_pval_fdr <- p.adjust(p = results.sat$H1v_pval, method = "fdr")
results.sat$H2v_pval_fdr <- p.adjust(p = results.sat$H2v_pval, method = "fdr")
results.sat$H3v_pval_fdr <- p.adjust(p = results.sat$H3v_pval, method = "fdr")
results.sat$H4v_pval_fdr <- p.adjust(p = results.sat$H4v_pval, method = "fdr")
results.sat$H5v_pval_fdr <- p.adjust(p = results.sat$H5v_pval, method = "fdr")
results.sat$H1c_pval_fdr <- p.adjust(p = results.sat$H1c_pval, method = "fdr")
results.sat$H2c_pval_fdr <- p.adjust(p = results.sat$H2c_pval, method = "fdr")
results.sat$H3c_pval_fdr <- p.adjust(p = results.sat$H3c_pval, method = "fdr")
results.sat$H4c_pval_fdr <- p.adjust(p = results.sat$H4c_pval, method = "fdr")
results.sat$NoAge_pval_fdr <- p.adjust(p = results.sat$NoAge_pval, method = "fdr")
results.sat$NoSex_pval_fdr <- p.adjust(p = results.sat$NoSex_pval, method = "fdr")
results.sat$NoeTIV_pval_fdr <- p.adjust(p = results.sat$NoeTIV_pval, method = "fdr")
results.sat$NoCohort_pval_fdr <- p.adjust(p = results.sat$NoCohort_pval, method = "fdr")

# Save output
write.csv(results.sat, "0-combined_saturated_output.csv", row.names = F)
