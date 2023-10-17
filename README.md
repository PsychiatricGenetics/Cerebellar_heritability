# Heritability of Cerebellar Volumes in Adolescent and Young Adult Twins
Analysis code for Strike et al., (2023). Heritability of Cerebellar Subregion Volumes in Adolescent and Young Adult Twins.

## Saturated
Saturated models to test twin modelling assumptions and covariate effects and estimate correlations between twin pairs. 
The script expects five zygosity groups (MZF, MZM, DZF, DZM, DZOS). Combined and QTIM scripts additionally expect one singleton sibling of twins.
Covariate coefficients and twin correlations are obtained from saturated models with equal means and variances across birth order, zygosity, and twins/singleton siblings of twins.

## Univariate
Variance component estimates (ACE or ADE model). Uses direct estimates of variance components (i.e., negative estimates are permitted).

## Bivariate
Bivariate variance component estimates (AE model), genetic, environmental, and phenotypic correlations. Uses direct estimates of variance components (i.e., negative estimates are permitted).
