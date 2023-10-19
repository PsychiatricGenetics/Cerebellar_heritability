# Definition variables(covariates) are hard coded
# OpenMx does not tolerate missing values for definition variables - recode any missing definition variables as -999
# Script will stop if missing definition variables are found (but a phenotype value exists) - avoids using -999 as an actual covariate value

Saturated_covariate <- function(phenotype, twin.data) {
  covariate <- c("Sex", "Age", "eTIVZ")
  nc <- length(covariate)
  for (x in covariate) {
    twin01.missing <- twin.data[, paste0(phenotype, "_01")][is.na(twin.data[, paste0(x, "_01")])]
    twin02.missing <- twin.data[, paste0(phenotype, "_02")][is.na(twin.data[, paste0(x, "_02")])]
    twin050.missing <- twin.data[, paste0(phenotype, "_050")][is.na(twin.data[, paste0(x, "_050")])]
    stopifnot(is.na(twin01.missing))
    stopifnot(is.na(twin02.missing))
    stopifnot(is.na(twin050.missing))
    twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
    twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
    twin.data[, paste0(x, "_050")][is.na(twin.data[, paste0(x, "_050")])] <- -999
  }
  
  # Select variables
  selVars <- c(paste0(phenotype, "_01"), paste0(phenotype, "_02"), paste0(phenotype, "_050"))
  covVars <- c(paste0(covariate, "_01"), paste0(covariate, "_02"), paste0(covariate, "_050"))
  useVars <- c(selVars, covVars)
  
  # Select data for analysis
  mzfData <- subset(twin.data, zyg == 1, useVars)
  mzmData <- subset(twin.data, zyg == 2, useVars)
  dzfData <- subset(twin.data, zyg == 3, useVars)
  dzmData <- subset(twin.data, zyg == 4, useVars)
  dzosData <- subset(twin.data, zyg == 5 | zyg == 6, useVars)
  
  # Set Starting Values
  nt <- 3
  nv <- length(phenotype)
  ntv <- nv * nt
  svMe <- mean(unlist(twin.data[, selVars]), na.rm = T)
  svVa <- var(unlist(twin.data[, selVars]), na.rm = T)
  svBe <- 0.01
  
  # Get raw means & sd for table (needed if using Z-score data)
  selVars.2 <- str_sub(string = phenotype, start = 1, end = -2)
  selVars.3 <- c(paste0(selVars.2, "_01"), paste0(selVars.2, "_02"), paste0(selVars.2, "_050"))
  mean.val <- mean(unlist(twin.data[, c(selVars.3[1], selVars.3[2], selVars.3[3])]), na.rm = T)
  sd.val <- sd(unlist(twin.data[, c(selVars.3[1], selVars.3[2], selVars.3[3])]), na.rm = T)
  
  # Create Matrices for Covariates and Linear Regression Coefficients
  defMZF <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02"), paste0("data.", covariate, "_050")), name = "defMZF")
  defMZM <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02"), paste0("data.", covariate, "_050")), name = "defMZM")
  defDZF <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02"), paste0("data.", covariate, "_050")), name = "defDZF")
  defDZM <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02"), paste0("data.", covariate, "_050")), name = "defDZM")
  defDZOS <- mxMatrix(type = "Full", nrow = nc, ncol = ntv, free = F, labels = c(paste0("data.", covariate, "_01"), paste0("data.", covariate, "_02"), paste0("data.", covariate, "_050")), name = "defDZOS")
  
  betaMZF <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaMZF")
  betaMZM <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaMZM")
  betaDZF <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaDZF")
  betaDZM <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaDZM")
  betaDZOS <- mxMatrix(type = "Full", nrow = 1, ncol = nc, free = TRUE, values = svBe, labels = c(paste0("beta", covariate)), name = "betaDZOS")
  
  # Algebra for expected Mean Matrices in MZ & DZ twins
  meanMZF <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mMZF1", "mMZF2", "mSib"), name = "meanMZF")
  meanMZM <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mMZM1", "mMZM2", "mSib"), name = "meanMZM")
  meanDZF <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mDZF1", "mDZF2", "mSib"), name = "meanDZF")
  meanDZM <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mDZM1", "mDZM2", "mSib"), name = "meanDZM")
  meanDZOS <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = c("mDZOS1", "mDZOS2", "mSib"), name = "meanDZOS")
  
  expMeanMZF <- mxAlgebra(expression = meanMZF + betaMZF %*% defMZF, name = "expMeanMZF")
  expMeanMZM <- mxAlgebra(expression = meanMZM + betaMZM %*% defMZM, name = "expMeanMZM")
  expMeanDZF <- mxAlgebra(expression = meanDZF + betaDZF %*% defDZF, name = "expMeanDZF")
  expMeanDZM <- mxAlgebra(expression = meanDZM + betaDZM %*% defDZM, name = "expMeanDZM")
  expMeanDZOS <- mxAlgebra(expression = meanDZOS + betaDZOS %*% defDZOS, name = "expMeanDZOS")
  
  # Twin Correlations
  expCorMZF <- mxAlgebra(cov2cor(expCovMZF), name = "expCorMZF")
  expCorMZM <- mxAlgebra(cov2cor(expCovMZM), name = "expCorMZM")
  expCorDZF <- mxAlgebra(cov2cor(expCovDZF), name = "expCorDZF")
  expCorDZM <- mxAlgebra(cov2cor(expCovDZM), name = "expCorDZM")
  expCorDZOS <- mxAlgebra(cov2cor(expCovDZOS), name = "expCorDZOS")
  
  # Create Algebra for expected Variance/Covariance Matrices
  expCovMZF <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vMZF1", "cMZF21", "cSib", "vMZF2", "cSib", "vSib"), name = "expCovMZF")
  expCovMZM <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vMZM1", "cMZM21", "cSib", "vMZM2", "cSib", "vSib"), name = "expCovMZM")
  expCovDZF <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vDZF1", "cDZF21", "cSib", "vDZF2", "cSib", "vSib"), name = "expCovDZF")
  expCovDZM <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vDZM1", "cDZM21", "cSib", "vDZM2", "cSib", "vSib"), name = "expCovDZM")
  expCovDZOS <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = c("vDZOS1", "cDZOS21", "cSib", "vDZOS2", "cSib", "vSib"), name = "expCovDZOS")
  
  # Data objects for Multiple Groups
  dataMZF <- mxData(observed = mzfData, type = "raw")
  dataMZM <- mxData(observed = mzmData, type = "raw")
  dataDZF <- mxData(observed = dzfData, type = "raw")
  dataDZM <- mxData(observed = dzmData, type = "raw")
  dataDZOS <- mxData(observed = dzosData, type = "raw")
  
  # Objective objects for Multiple Groups
  funML <- mxFitFunctionML()
  objMZF <- mxExpectationNormal(covariance = "expCovMZF", means = "expMeanMZF", dimnames = selVars)
  objMZM <- mxExpectationNormal(covariance = "expCovMZM", means = "expMeanMZM", dimnames = selVars)
  objDZF <- mxExpectationNormal(covariance = "expCovDZF", means = "expMeanDZF", dimnames = selVars)
  objDZM <- mxExpectationNormal(covariance = "expCovDZM", means = "expMeanDZM", dimnames = selVars)
  objDZOS <- mxExpectationNormal(covariance = "expCovDZOS", means = "expMeanDZOS", dimnames = selVars)
  
  # Combine Groups
  modelMZF <- mxModel("MZF", meanMZF, betaMZF, defMZF, expMeanMZF, expCorMZF, expCovMZF, dataMZF, objMZF, funML)
  modelMZM <- mxModel("MZM", meanMZM, betaMZM, defMZM, expMeanMZM, expCorMZM, expCovMZM, dataMZM, objMZM, funML)
  modelDZF <- mxModel("DZF", meanDZF, betaDZF, defDZF, expMeanDZF, expCorDZF, expCovDZF, dataDZF, objDZF, funML)
  modelDZM <- mxModel("DZM", meanDZM, betaDZM, defDZM, expMeanDZM, expCorDZM, expCovDZM, dataDZM, objDZM, funML)
  modelDZOS <- mxModel("DZOS", meanDZOS, betaDZOS, defDZOS, expMeanDZOS, expCorDZOS, expCovDZOS, dataDZOS, objDZOS, funML)
  minus2ll <- mxAlgebra(MZF.objective + MZM.objective + DZF.objective + DZM.objective + DZOS.objective, name = "minus2sumloglikelihood")
  obj <- mxFitFunctionAlgebra("minus2sumloglikelihood")
  ciCor <- mxCI(c("MZF.expCorMZF[2,1]", "DZF.expCorDZF[2,1]"))
  twinSatModel <- mxModel("twinSat", modelMZF, modelMZM, modelDZF, modelDZM, modelDZOS, minus2ll, obj, ciCor)
  
  # Run Saturated Model
  twinSatFit <- mxTryHard(twinSatModel, intervals = T, extraTries = 100)
  twinSatSumm <- summary(twinSatFit)
  
  #### Sub-models ####
  # MODELS TO TEST HETEROGENIETY OF MEANS
  # Constrain expected Means to be equal across twin order
  MModel1 <- twinSatModel
  MModel1 <- omxSetParameters(MModel1, label = c("mMZF1", "mMZF2"), newlabels = "mMZF")
  MModel1 <- omxSetParameters(MModel1, label = c("mMZM1", "mMZM2"), newlabels = "mMZM")
  MModel1 <- omxSetParameters(MModel1, label = c("mDZF1", "mDZF2"), newlabels = "mDZF")
  MModel1 <- omxSetParameters(MModel1, label = c("mDZM1", "mDZM2"), newlabels = "mDZM")
  MModel1Fit <- mxTryHard(MModel1, intervals = F, extraTries = 50)
  
  # Constrain expected Means to be equal across same sex groups
  MModel2 <- MModel1
  MModel2 <- omxSetParameters(MModel2, label = c("mMZF", "mDZF"), newlabels = "mFemale")
  MModel2 <- omxSetParameters(MModel2, label = c("mMZM", "mDZM"), newlabels = "mMale")
  MModel2Fit <- mxTryHard(MModel2, intervals = F, extraTries = 50)
  
  # Constrain expected Means to be equal across opposite sex groups. DZOS indid 1 = female, indid 2 = male
  MModel3 <- MModel2
  MModel3 <- omxSetParameters(MModel3, label = c("mDZOS1", "mFemale"), newlabels = "mF")
  MModel3 <- omxSetParameters(MModel3, label = c("mDZOS2", "mMale"), newlabels = "mM")
  MModel3Fit <- mxTryHard(MModel3, intervals = F, extraTries = 50)
  
  # Constrain expected Means to be equal across sex
  MModel4 <- MModel3
  MModel4 <- omxSetParameters(MModel4, label = c("mF", "mM"), newlabels = "mTwin")
  MModel4Fit <- mxTryHard(MModel4, intervals = F, extraTries = 50)
  
  # Constrain expected Means to be equal across twin/sibling
  MModel5 <- MModel4
  MModel5 <- omxSetParameters(MModel5, label = c("mTwin", "mSib"), newlabels = "m")
  MModel5Fit <- mxTryHard(MModel5, intervals = F, extraTries = 50)
  
  # MODELS TO TEST HETEROGENIETY OF VARIANCES
  # Constrain expected Variances to be equal across twin order
  VModel1 <- MModel5
  VModel1 <- omxSetParameters(VModel1, label = c("vMZF1", "vMZF2"), newlabels = "vMZF")
  VModel1 <- omxSetParameters(VModel1, label = c("vMZM1", "vMZM2"), newlabels = "vMZM")
  VModel1 <- omxSetParameters(VModel1, label = c("vDZF1", "vDZF2"), newlabels = "vDZF")
  VModel1 <- omxSetParameters(VModel1, label = c("vDZM1", "vDZM2"), newlabels = "vDZM")
  VModel1Fit <- mxTryHard(VModel1, intervals = F, extraTries = 50)
  
  # Constrain expected Variances to be equal across same sex groups
  VModel2 <- VModel1
  VModel2 <- omxSetParameters(VModel2, label = c("vMZF", "vDZF"), newlabels = "vFemale")
  VModel2 <- omxSetParameters(VModel2, label = c("vMZM", "vDZM"), newlabels = "vMale")
  VModel2Fit <- mxTryHard(VModel2, intervals = F, extraTries = 50)
  
  # Constrain expected Variances to be equal across opposite sex groups
  VModel3 <- VModel2
  VModel3 <- omxSetParameters(VModel3, label = c("vDZOS1", "vFemale"), newlabels = "vF")
  VModel3 <- omxSetParameters(VModel3, label = c("vDZOS2", "vMale"), newlabels = "vM")
  VModel3Fit <- mxTryHard(VModel3, intervals = F, extraTries = 50)
  
  # Constrain expected Variances to be equal across sex
  VModel4 <- VModel3
  VModel4 <- omxSetParameters(VModel4, label = c("vF", "vM"), newlabels = "vTwin")
  VModel4Fit <- mxTryHard(VModel4, intervals = F, extraTries = 50)
  
  # Constrain expected Variances to be equal across twin/sib
  VModel5 <- VModel4
  VModel5 <- omxSetParameters(VModel5, label = c("vTwin", "vSib"), newlabels = "v")
  VModel5Fit <- mxTryHard(VModel5, intervals = F, extraTries = 50)
  
  # MODELS TO TEST HETEROGENIETY OF COVARIANCES
  # Constrain expected covs to be equal within zygosity
  CModel1 <- VModel5
  CModel1 <- omxSetParameters(CModel1, label = c("cMZF21", "cMZM21"), values = 0, newlabels = "cMZ")
  CModel1 <- omxSetParameters(CModel1, label = c("cDZF21", "cDZM21"), values = 0, newlabels = "cDZ21")
  CModel1Fit <- mxTryHard(CModel1, intervals = F, extraTries = 50)
  
  # Constrain expected covs to be equal across same and & opposite sex groups
  CModel2 <- CModel1
  CModel2 <- omxSetParameters(CModel2, label = c("cDZOS21", "cDZ21", "cSib"), values = 0, newlabels = "cDZ")
  CModel2Fit <- mxTryHard(CModel2, intervals = T, extraTries = 50)
  CModel2Summ <- summary(CModel2Fit)
  
  # Constrain expected covs to be equal between MZ and DZ groups
  CModel3 <- CModel2
  CModel3 <- omxSetParameters(CModel3, label = c("cMZ", "cDZ"), values = 0, newlabels = "c")
  CModel3Fit <- mxTryHard(CModel3, intervals = F, extraTries = 50)
  
  # Drop cov to zero
  # H4c set all covariances to zero and tested for effects of familial aggregation
  CModel4 <- CModel3
  CModel4 <- omxSetParameters(CModel4, label = "c", free = FALSE, values = 0)
  CModel4Fit <- mxTryHard(CModel4, intervals = F, extraTries = 50)
  
  # Model fitting
  H1m <- mxCompare(twinSatFit, MModel1Fit)
  H2m <- mxCompare(MModel1Fit, MModel2Fit)
  H3m <- mxCompare(MModel2Fit, MModel3Fit)
  H4m <- mxCompare(MModel3Fit, MModel4Fit)
  H5m <- mxCompare(MModel4Fit, MModel5Fit)
  
  H1v <- mxCompare(MModel5Fit, VModel1Fit)
  H2v <- mxCompare(VModel1Fit, VModel2Fit)
  H3v <- mxCompare(VModel2Fit, VModel3Fit)
  H4v <- mxCompare(VModel3Fit, VModel4Fit)
  H5v <- mxCompare(VModel4Fit, VModel5Fit)
  
  H1c <- mxCompare(VModel5Fit, CModel1Fit)
  H2c <- mxCompare(CModel1Fit, CModel2Fit)
  H3c <- mxCompare(CModel2Fit, CModel3Fit)
  H4c <- mxCompare(CModel3Fit, CModel4Fit)
  
  #########################################################################################
  # Test Covariates
  # Drop Sex
  noSex <- CModel2
  noSex <- omxSetParameters(model = noSex, labels = "betaSex", free = FALSE, values = 0, name = "noSexModel")
  noSexFit <- mxTryHard(noSex, extraTries = 100, intervals = F)
  noSexSumm <- summary(noSexFit)
  
  # Drop Age
  noAge <- CModel2
  noAge <- omxSetParameters(model = noAge, labels = "betaAge", free = FALSE, values = 0, name = "noAgeModel")
  noAgeFit <- mxTryHard(noAge, extraTries = 100, intervals = F)
  noAgeSumm <- summary(noAgeFit)
  
  # Drop eTIV
  noETIV <- CModel2
  noETIV <- omxSetParameters(model = noETIV, labels = "betaeTIVZ", free = FALSE, values = 0, name = "noETIVModel")
  noETIVFit <- mxTryHard(noETIV, extraTries = 100, intervals = F)
  noETIVSumm <- summary(noETIVFit)
  
  # Model fitting
  noSexLRT <- mxCompare(CModel2Fit, noSexFit)
  noAgeLRT <- mxCompare(CModel2Fit, noAgeFit)
  noETIVLRT <- mxCompare(CModel2Fit, noETIVFit)

  # ACE or ADE model
  rMZ <- round(CModel2Fit$MZF$expCorMZF$result[2, 1],2)
  rDZ <- round(2*CModel2Fit$DZF$expCorDZF$result[2, 1],2)
  if (rMZ < rDZ) {
    model.select <- "ACE"
  } else
  {
    model.select <- "ADE"
  }
  
  Results <- data.frame(
    Variable = phenotype,
    ObservedStatistics = twinSatSumm$observedStatistics,
    Mean = mean.val,
    SD = sd.val,
    rMZF = VModel4Fit$MZF$expCorMZF$result[2, 1],
    rMZM = VModel4Fit$MZM$expCorMZM$result[2, 1],
    rDZF = VModel4Fit$DZF$expCorDZF$result[2, 1],
    rDZM = VModel4Fit$DZM$expCorDZM$result[2, 1],
    rDZOS = VModel4Fit$DZOS$expCorDZOS$result[2, 1],
    rMZ = CModel2Fit$MZF$expCorMZF$result[2, 1],
    rDZ = CModel2Fit$DZF$expCorDZF$result[2, 1],
    rMZ_95CI = paste0(sprintf("%.2f", round(CModel2Summ$CI$estimate[1], 2)), " (", sprintf("%.2f", round(CModel2Summ$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(CModel2Summ$CI$ubound[1], 2)), ")"),
    rDZ_95CI = paste0(sprintf("%.2f", round(CModel2Summ$CI$estimate[2], 2)), " (", sprintf("%.2f", round(CModel2Summ$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(CModel2Summ$CI$ubound[2], 2)), ")"),
    rTwin_model = model.select,
    sexEstimate = CModel2Summ$parameters$Estimate[2],
    ageEstimate = CModel2Summ$parameters$Estimate[3],
    eTIVEstimate = CModel2Summ$parameters$Estimate[4],
    NoSex_pval = noSexLRT$p[2],
    NoAge_pval = noAgeLRT$p[2],
    NoeTIV_pval = noETIVLRT$p[2],
    AIC_Sat = twinSatSumm$AIC.Mx,
    CODE_Sat = twinSatFit$output$status$code,
    H1m_pval = H1m$p[2],
    H2m_pval = H2m$p[2],
    H3m_pval = H3m$p[2],
    H4m_pval = H4m$p[2],
    H5m_pval = H5m$p[2],
    H1v_pval = H1v$p[2],
    H2v_pval = H2v$p[2],
    H3v_pval = H3v$p[2],
    H4v_pval = H4v$p[2],
    H5v_pval = H5v$p[2],
    H1c_pval = H1c$p[2],
    H2c_pval = H2c$p[2],
    H3c_pval = H3c$p[2],
    H4c_pval = H4c$p[2],
    stringsAsFactors = FALSE
  )
  return(Results)
}
