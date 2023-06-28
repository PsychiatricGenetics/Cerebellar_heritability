# Lachlan T Strike
# Direct estimates of variance components
# Based on scripts shared by Hermine Maes (https://hermine-maes.squarespace.com)

rm(list = ls())

#### Models ####
Univariate_ADE <- function(phenotype, twin.data, covariate) {
  nc <- length(covariate)
  # OpenMx does not tolerate missing values for definition variables.
  # Recode any missing definition variables as -999
  # BUT!!! Make sure there are not any cases of missing definition variables
  # with a phenotype present
  for (x in covariate) {
    twin01.missing <- twin.data[, paste0(phenotype, "_01")][is.na(twin.data[, paste0(x, "_01")])]
    twin02.missing <- twin.data[, paste0(phenotype, "_02")][is.na(twin.data[, paste0(x, "_02")])]
    stopifnot(is.na(twin01.missing))
    stopifnot(is.na(twin02.missing))
    twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
    twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
  }
  
  # Select variables
  selVars <- c(paste0(phenotype, "_01"), paste0(phenotype, "_02"))
  covVars <- c(paste0(covariate, "_01"), paste0(covariate, "_02"))
  useVars <- c(selVars, covVars)
  
  # Select data for analysis
  mzData <- subset(twin.data, zyg < 3, useVars)
  dzData <- subset(twin.data, zyg > 2, useVars)
  
  # Set Starting Values
  nt <- 2
  nv <- length(phenotype)
  ntv <- nv * nt
  svPa <- sqrt(var(c(twin.data[, paste0(phenotype, "_01")], twin.data[, paste0(phenotype, "_02")]), na.rm = T) / 3)
  svMe <- mean(unlist(twin.data[, selVars]), na.rm = T)
  
  # ------------------------------------------------------------------------------
  # PREPARE GENETIC MODEL
  # ------------------------------------------------------------------------------
  # Create Matrices for Covariates and linear Regression Coefficients
  defL <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02")), name = "defL")
  pathBl <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = 0.1, labels = c(paste0("beta", covariate)), name = "bl")
  
  # Create Algebra for expected Mean Matrices
  meanG <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, label = "mean", name = "meanG")
  expMean <- mxAlgebra(expression = meanG + bl %*% defL, name = "expMeanG")
  
  # Create Matrices for Variance Components
  covA <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" )
  covC <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VD11", name="VD" )
  covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VE11", name="VE" )
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP <- mxAlgebra( expression= VA+VD+VE, name="V" )
  covMZ <- mxAlgebra( expression= VA+VD, name="cMZ" )
  covDZ <- mxAlgebra( expression= 0.5%x%VA+0.25%x%VD, name="cDZ" )
  
  expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData( observed=mzData, type="raw" )
  dataDZ <- mxData( observed=dzData, type="raw" )
  
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal( covariance="expCovMZ", means="expMeanG", dimnames=selVars )
  expDZ <- mxExpectationNormal( covariance="expCovDZ", means="expMeanG", dimnames=selVars )
  funML <- mxFitFunctionML()
  
  # Create Model Objects for Multiple Groups
  pars <- list( pathBl, meanG, covA, covC, covE, covP )
  defs <- list( defL )
  modelMZ <- mxModel( pars, defs, expMean, covMZ, covDZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
  multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
  
  # Create Algebra for Variance Components
  rowUS <- rep('US',nv)
  colUS <- rep(c('VA','VD','VE','SA','SD','SE'),each=nv)
  estUS <- mxAlgebra( expression=cbind(VA,VD,VE,VA/V,VD/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
  
  # Create Confidence Interval Objects
  ciADE <- mxCI( "US[1,1:6]" )
  
  # Build Model with Confidence Intervals
  modelADE <- mxModel( "oneADEvca", pars, modelMZ, modelDZ, multi, estUS, ciADE )
  
  # ------------------------------------------------------------------------------
  fitADE <- mxTryHard(modelADE, intervals = T, extraTries = 50)
  sumADE <- summary(fitADE)
  
  # Run AE model
  modelAE <- mxModel( modelADE, name="oneAEvca" )
  modelAE <- omxSetParameters( modelAE, labels="VD11", free=FALSE, values=0 )
  fitAE <- mxTryHard( modelAE, intervals=T, extraTries = 50)
  sumAE <- summary(fitAE)
  
  # Run E model
  modelE <- mxModel( modelAE, name="oneEvca" )
  modelE <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0)
  fitE <- mxTryHard( modelE, intervals=T, extraTries = 50)
  sumE <- summary(fitE)
  
  # Model Comparison
  AW <- omxAkaikeWeights(models = list(fitADE, fitAE, fitE))
  model.comparison <- mxCompare( fitADE, nested <- list(fitAE, fitE))
  
  #### Create results table ####
  Results <- data.frame(
    Variable = phenotype,
    ObservedStatistics = sumADE$observedStatistics,
    ADE_A_95CI = paste0(sprintf("%.2f", round(sumADE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumADE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumADE$CI$ubound[4], 2)), ")"),
    ADE_D_95CI = paste0(sprintf("%.2f", round(sumADE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumADE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumADE$CI$ubound[5], 2)), ")"),
    ADE_E_95CI = paste0(sprintf("%.2f", round(sumADE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumADE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumADE$CI$ubound[6], 2)), ")"),
    AE_A_95CI = paste0(sprintf("%.2f", round(sumAE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[4], 2)), ")"),
    AE_E_95CI = paste0(sprintf("%.2f", round(sumAE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[6], 2)), ")"),
    ADE_AIC = sumADE$AIC.Mx,
    AE_AIC = sumAE$AIC.Mx,
    E_AIC = sumE$AIC.Mx,
    AIC_BestFit = AW$model[1],
    dropD = model.comparison$p[2],
    dropAD = model.comparison$p[3],
    AE_A = fitAE$US$result[4],
    AE_E = fitAE$US$result[6],
    stringsAsFactors = F)
    return(Results)
}

#### Run ####
library(OpenMx)
library(tidyverse)
setwd("C:/GitHub/Cerebellar_heritability")
my.data <- readRDS("qtab_cb_familywise.RDS")
source("miFunctions.R")

# Single phenotype
result.single <- as_tibble(Univariate_ADE(phenotype = "Left_I_IIIZ", twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ")))
result.single %>% select(Variable, AIC_BestFit, ADE_A_95CI, ADE_D_95CI, ADE_E_95CI, AE_A_95CI, AE_E_95CI, AE_A, AE_E)

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
result.multiple <- lapply(phenotype.arrayZ, Univariate_ADE, twin.data = my.data, covariate = c("Sex", "Age", "eTIVZ")) %>% bind_rows()
result.multiple <- as_tibble(result.multiple)

write.csv(result.multiple, "qtab_univariate_ADE_output.csv", row.names = F)
