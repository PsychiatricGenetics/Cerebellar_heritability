# Definition variables(covariates) are hard coded
# OpenMx does not tolerate missing values for definition variables - recode any missing definition variables as -999
# Script will stop if missing definition variables are found (but a phenotype value exists) - avoids using -999 as an actual covariate value

Bivariate_AE <- function(phenotypes, twin.data){
  print(phenotypes)
  covariates = c("Age", "Sex", "eTIVZ")
  nc <- length(covariates)
  for (y in phenotypes){
    for (x in covariates) {
      twin01.missing <- twin.data[, paste0(y, "_01")][is.na(twin.data[, paste0(x, "_01")])]
      twin02.missing <- twin.data[, paste0(y, "_02")][is.na(twin.data[, paste0(x, "_02")])]
      twin050.missing <- twin.data[, paste0(y, "_050")][is.na(twin.data[, paste0(x, "_050")])]
      stopifnot(is.na(twin01.missing))
      stopifnot(is.na(twin02.missing))
      stopifnot(is.na(twin050.missing))
      twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
      twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
      twin.data[, paste0(x, "_050")][is.na(twin.data[, paste0(x, "_050")])] <- -999
    }
  }
  
  # Select variables
  selVars <- c(paste0(phenotypes, "_01"), paste0(phenotypes, "_02"), paste0(phenotypes, "_050"))
  covVars <- c(paste0(covariates, "_01"), paste0(covariates, "_02"), paste0(covariates, "_050"))
  nv <- length(phenotypes)
  nt <- 3
  ntv <- nv * nt
  
  # Covariate labels
  ageLabels <- paste(rep(c("Age"), times=nt), rep(c(1,2,50), each=nv), sep="_0")
  sexLabels <- paste(rep(c("Sex"), times=nt), rep(c(1,2,50), each=nv), sep="_0")
  eTIVLabels <- paste(rep(c("eTIVZ"), times=nt), rep(c(1,2,50), each=nv), sep="_0")
  
  # Beta labels
  beta_age_mean_labels <- rep(paste0("beta_age_mean_coef_", 1:nv), times=nt)
  beta_sex_mean_labels <- rep(paste0("beta_sex_mean_coef_", 1:nv), times=nt)
  beta_etiv_mean_labels <- rep(paste0("beta_etiv_mean_coef_", 1:nv), times=nt)
  
  # Select Data for Analysis
  useVars <- c(selVars, covVars)
  mzData <- subset(twin.data, zyg<=2, useVars)
  dzData <- subset(twin.data, zyg>=3, useVars)
  
  # Set Starting Values 
  svMe <- colMeans(twin.data[, selVars], na.rm = T)[1:nv]
  svPa <- var(twin.data[, selVars], na.rm = T)[1:nv]/3
  svPe <- svPa*2
  # ------------------------------------------------------------------------------
  
  # ------------------------------------------------------------------------------
  # Prepare Model
  # Intercept Coefficient for Means
  mean      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean", phenotypes), name="InterceptMean" )
  
  # Create Betas for Covariates
  betaAge <- mxMatrix( type = "Full", nrow = 1, ncol = nv*nt, free = TRUE, values = 0.01, labels = beta_age_mean_labels, name = "BetaAge" )
  betaSex <- mxMatrix( type = "Full", nrow = 1, ncol = nv*nt, free = TRUE, values = 0.01, labels = beta_sex_mean_labels, name = "BetaSex" )
  betaeTIV <- mxMatrix( type = "Full", nrow = 1, ncol = nv*nt, free = TRUE, values = 0.01, labels = beta_etiv_mean_labels, name = "BetaeTIV" )
  
  # Create Definition Variables for Covariates
  defAge <- mxMatrix( type = "Full", nrow = 1, ncol = nv*nt, free = F, labels = paste0(rep("data.", times=nv*nt), ageLabels), name = "DefAge" )
  defSex <- mxMatrix( type = "Full", nrow = 1, ncol = nv*nt, free = F, labels = paste0(rep("data.", times=nv*nt), sexLabels), name = "DefSex" )
  defeTIV <- mxMatrix( type = "Full", nrow = 1, ncol = nv*nt, free = F, labels = paste0(rep("data.", times=nv*nt), eTIVLabels), name = "DefeTIV" )
  
  # Expected Means
  expMean <- mxAlgebra(expression = InterceptMean + BetaSex*DefSex + BetaAge*DefAge + BetaeTIV*DefeTIV, name = "expMean")
  
  # Create Matrices for Variance Components
  covA <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VA", nv), name="VA" ) 
  covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VE", nv), name="VE" ) 
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins 
  covP <- mxAlgebra( expression= VA+VE, name="V" ) 
  covMZ <- mxAlgebra( expression= VA, name="cMZ" )
  covDZ <- mxAlgebra( expression= 0.5%x%VA, name="cDZ" )
  expCovMZ <- mxAlgebra( expression= rbind(cbind(V, cMZ, cDZ), 
                                           cbind(cMZ, V, cDZ),
                                           cbind(cDZ, cDZ, V)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind(cbind(V, cDZ, cDZ),
                                           cbind(cDZ, V, cDZ),
                                           cbind(cDZ, cDZ, V)), name="expCovDZ" )
  
  # Create Data Objects for Multiple Groups
  dataMZ    <- mxData( observed=mzData, type="raw" )
  dataDZ    <- mxData( observed=dzData, type="raw" )
  
  # Create Expectation Objects for Multiple Groups
  expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
  expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )
  funML     <- mxFitFunctionML()
  
  # Create Algebra for Standardization
  matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
  invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")
  
  # AE Estimates
  stVA <- mxAlgebra(VA/V, name = "StVA")
  stVE <- mxAlgebra(VE/V, name = "StVE")
  
  # Calculate correlations
  corV      <- mxAlgebra( expression=cov2cor(V), name ="rV" )
  corA      <- mxAlgebra( expression=cov2cor(VA), name ="rA" )
  corE      <- mxAlgebra( expression=cov2cor(VE), name ="rE" )
  
  # Create Model Objects for Multiple Groups
  pars      <- list(covA, covE, covP, stVA, stVE)
  modelMZ   <- mxModel( name="MZ", pars, mean, expMean, betaAge, defAge, betaSex, defSex, betaeTIV, defeTIV, covMZ, covDZ, expCovMZ, dataMZ, expMZ, funML )
  modelDZ   <- mxModel( name="DZ", pars, mean, expMean, betaAge, defAge, betaSex, defSex, betaeTIV, defeTIV, covDZ, expCovDZ, dataDZ, expDZ, funML )
  multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
  
  # Build Model with Confidence Intervals
  calc      <- list( matI, invSD, corV, corA, corE )
  ciAE <- mxCI(c("rA[2,1]", "rE[2,1]", "rV[2,1]"))
  modelAE  <- mxModel( "mulAEvc", pars, modelMZ, modelDZ, multi, calc, ciAE)
  # ------------------------------------------------------------------------------
  
  # Run AE Model
  fitAE <- mxTryHard(modelAE, intervals = F, extraTries = 50)
  sumAE <- summary(fitAE)
  # ------------------------------------------------------------------------------
  
  # Test rV significance
  dropRV <- mxModel( modelAE, name="no_rV" )
  dropRV <- omxSetParameters( dropRV, labels=c("VA21", "VE21"), free=FALSE, values=0 )
  fit_dropRV <- mxTryHard(dropRV, intervals = F, extraTries = 50)
  sum_dropRV <- summary(dropRV)
  # ------------------------------------------------------------------------------
  
  # Test rA significance
  dropRA <- mxModel( modelAE, name="no_rA" )
  dropRA <- omxSetParameters( dropRA, labels="VA21", free=FALSE, values=0 )
  fit_dropRA <- mxTryHard(dropRA, intervals = F, extraTries = 50)
  sum_dropRA <- summary(dropRA)
  # ------------------------------------------------------------------------------
  
  # Test rE significance
  dropRE <- mxModel( modelAE, name="no_rE" )
  dropRE <- omxSetParameters( dropRE, labels="VE21", free=FALSE, values=0 )
  fit_dropRE <- mxTryHard(dropRE, intervals = F, extraTries = 50)
  sum_dropRE <- summary(dropRE)
  # ------------------------------------------------------------------------------
  
  model.comparison <- mxCompare( fitAE, nested <- list(fit_dropRV, fit_dropRA, fit_dropRE))
  
  #### Create results table ####
  Results <- data.frame(
    Variable = paste0(phenotypes[1], "_corr_", phenotypes[2]),
    ObservedStatistics = sumAE$observedStatistics,
    AE_A11 = fitAE$MZ$StVA$result[1, 1],
    AE_E11 = fitAE$MZ$StVE$result[1, 1],
    AE_A22 = fitAE$MZ$StVA$result[2, 2],
    AE_E22 = fitAE$MZ$StVE$result[2, 2],
    AE_rV = fitAE$rV$result[2, 1],
    AE_rA = fitAE$rA$result[2, 1],
    AE_rE = fitAE$rE$result[2, 1],
    AE_rV_sig = model.comparison$p[2],
    AE_rA_sig = model.comparison$p[3],
    AE_rE_sig = model.comparison$p[4],
    stringsAsFactors = F)
  return(Results)
}
