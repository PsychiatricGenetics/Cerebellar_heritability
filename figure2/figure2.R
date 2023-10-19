setwd("C:/GitHub/Cerebellar_heritability/figure2")

library(stringr)
library(tidyverse)
library(corrplot)
library(viridis)

#### Notes ####
# The upper & lower matrices should be similar, but not identical
# Each region pairing has two versions (i.e., the correlation between A & B
# has a model where the order is A-B, and another where the order is B-A)
# Check to see if there are any differences between the upper & lower matrices
################################################################################

rm(list=ls())

rV <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(rV) <- c("Phen1", "Phen2", "rV")
rV_A <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(rV_A) <- c("Phen1", "Phen2", "rV_A")
rV_E <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(rV_E) <- c("Phen1", "Phen2", "rV_E")

list1 <- c("Left_I_III", "Right_I_III", "Left_IV", "Right_IV", "Left_V", "Right_V",
           "Left_VI",  "Right_VI", "Left_Crus_I", "Right_Crus_I", "Left_Crus_II", "Right_Crus_II", 
           "Left_VIIB", "Right_VIIB", "Left_VIIIA", "Right_VIIIA", "Left_VIIIB", "Right_VIIIB",
           "Left_IX", "Right_IX", "Left_X", "Right_X", "Vermis_VI", "Vermis_VII", "Vermis_VIII",
           "Vermis_IX", "Vermis_X", "Corpus_Medullare", "Total_Cerebel_Vol")

for (i in list1){
  rMat <- read.csv(paste0("combined_bivariate_AE_", i, "Z_output.csv"))
  #rMat$rV_A <- rMat$AE_rA
  #rMat$rV_E <- rMat$AE_rE
  rMat$rV_A = sqrt(rMat$AE_A11) * rMat$AE_rA * sqrt(rMat$AE_A22)
  rMat$rV_E = sqrt(rMat$AE_E11) * rMat$AE_rE * sqrt(rMat$AE_E22)
  rMat$Phen1 <- str_split(rMat$Variable , "_corr_", simplify=TRUE)[,1]
  rMat$Phen2 <- str_split(rMat$Variable , "_corr_", simplify=TRUE)[,2]
  rMat$Phen1 <- str_split(rMat$Phen1, "covariate", simplify=TRUE)[,1]
  rMat$Phen2 <- str_split(rMat$Phen2, "covariate", simplify=TRUE)[,1]
  rMat$Phen1 <- substr(rMat$Phen1, start = 1, nchar(rMat$Phen1)-1)
  rMat$Phen2 <- substr(rMat$Phen2, start = 1, nchar(rMat$Phen2)-1)
  
  # Create the diag correlations
  rMat[29, "Phen1"] <- rMat$Phen1[1]
  rMat[29, "Phen2"] <- rMat$Phen1[1]
  rMat[29, "AE_rV"] <- 1
  rMat[29, "rV_A"] <- 1
  rMat[29, "rV_E"] <- 1
  
  rV <- rbind(rV, rMat[, c('Phen1', 'Phen2', 'AE_rV', "AE_rV_sig")])
  rV_A <- rbind(rV_A, rMat[, c('Phen1', 'Phen2', 'rV_A', "AE_rA_sig")])
  rV_E <- rbind(rV_E, rMat[, c('Phen1', 'Phen2', 'rV_E', "AE_rE_sig")])
  
}

rV.2 <- rV
rV_A.2 <- rV_A
rV_E.2 <- rV_E

rV_pvals <- rV
rV <-
  rV %>% select(Phen1, Phen2, AE_rV)
rV_pvals$AE_rV_sig_fdr <- p.adjust(rV_pvals$AE_rV_sig, method = "BH")
rV_pvals %>% filter(AE_rV_sig_fdr < 0.05)
rV <-
  rV %>% spread(Phen1, AE_rV)
rV_pvals <-
  rV_pvals %>% select(Phen1, Phen2, AE_rV_sig_fdr)
rV_pvals <-
  rV_pvals %>% spread(Phen1, AE_rV_sig_fdr)

rV_A_pvals <- rV_A
rV_A <-
  rV_A %>% select(Phen1, Phen2, rV_A)
rV_A_pvals$AE_rA_sig_fdr <- p.adjust(rV_A_pvals$AE_rA_sig, method = "BH")
rV_A_pvals %>% filter(AE_rA_sig_fdr < 0.05)
rV_A <-
  rV_A %>% spread(Phen1, rV_A)
rV_A_pvals <-
  rV_A_pvals %>% select(Phen1, Phen2, AE_rA_sig_fdr)
rV_A_pvals <-
  rV_A_pvals %>% spread(Phen1, AE_rA_sig_fdr)

rV_E_pvals <- rV_E
rV_E <-
  rV_E %>% select(Phen1, Phen2, rV_E)
rV_E_pvals$AE_rE_sig_fdr <- p.adjust(rV_E_pvals$AE_rE_sig, method = "BH")
rV_E_pvals %>% filter(AE_rE_sig_fdr < 0.05)
rV_E <-
  rV_E %>% spread(Phen1, rV_E)
rV_E_pvals <-
  rV_E_pvals %>% select(Phen1, Phen2, AE_rE_sig_fdr)
rV_E_pvals <-
  rV_E_pvals %>% spread(Phen1, AE_rE_sig_fdr)

#### Check for differences between upper & lower matrices ####
rV.2 <-
  rV.2 %>% select(Phen1, Phen2, AE_rV)
rV.2 <-
  rV.2 %>% spread(Phen2, AE_rV)

corr.comparison <- matrix(nrow = 29, ncol = 29, data = NA)
colnames(corr.comparison) <- list1
row.names(corr.comparison) <- list1
for (i in list1){
  corr.comparison[, i] <- abs(rV[, i] - rV.2[, i])
}

rV_A.2 <-
  rV_A.2 %>% select(Phen1, Phen2, rV_A)
rV_A.2 <-
  rV_A.2 %>% spread(Phen2, rV_A)

corr.comparison.rA <- matrix(nrow = 29, ncol = 29, data = NA)
colnames(corr.comparison.rA) <- list1
row.names(corr.comparison.rA) <- list1
for (i in list1){
  corr.comparison.rA[, i] <- abs(rV_A[, i] - rV_A.2[, i])
}

rV_E.2 <-
  rV_E.2 %>% select(Phen1, Phen2, rV_E)
rV_E.2 <-
  rV_E.2 %>% spread(Phen2, rV_E)

corr.comparison.rE <- matrix(nrow = 29, ncol = 29, data = NA)
colnames(corr.comparison.rE) <- list1
row.names(corr.comparison.rE) <- list1
for (i in list1){
  corr.comparison.rE[, i] <- abs(rV_E[, i] - rV_E.2[, i])
}

corr.comparison[which(corr.comparison > 0.00001)]
corr.comparison.rA[which(corr.comparison.rA > 0.00001)]
corr.comparison.rE[which(corr.comparison.rE > 0.00001)]
# No differences at 0.1e-04

# Row/col order
rc.order <- read.table("combined_bivariate_AE_matrices_order.txt", sep = '\t', header = T)

#### Phenotypic Correlation ####
rV.mat1 <- left_join(rV, rc.order, by=c("Phen2"="Region"))
rV.mat1$Phen2 <- rV.mat1$Region2
rV.mat1 <- 
  rV.mat1 %>% relocate(all_of(c("Phen2", list1)))
rV.mat1 <- 
  rV.mat1 %>% arrange(Order)
rV.mat1 <- 
  rV.mat1 %>% select(-Order, -Region2)
colnames(rV.mat1) <- c("Phen2", rV.mat1$Phen2)

row.names(rV.mat1) <- rV.mat1[,1]
rV.mat1 <- rV.mat1[,-1]

rV_matrix <- matrix(data = NA, nrow = dim(rV.mat1)[1], ncol = dim(rV.mat1)[2])
colnames(rV_matrix) <- colnames(rV.mat1)
rownames(rV_matrix) <- rownames(rV.mat1)

for (i in 1:dim(rV.mat1)[2]) {
  rV_matrix[,i] <- c(as.numeric(rV.mat1[[i]]))
}
diag(rV_matrix) <- 1

#### Phenotypic Correlation p-values ####
rV.pvals.mat1 <- left_join(rV_pvals, rc.order, by=c("Phen2"="Region"))
rV.pvals.mat1$Phen2 <- rV.pvals.mat1$Region2
rV.pvals.mat1 <- 
  rV.pvals.mat1 %>% relocate(all_of(c("Phen2", list1)))
rV.pvals.mat1 <- 
  rV.pvals.mat1 %>% arrange(Order)
rV.pvals.mat1 <- 
  rV.pvals.mat1 %>% select(-Order, -Region2)
colnames(rV.pvals.mat1) <- c("Phen2", rV.pvals.mat1$Phen2)

row.names(rV.pvals.mat1) <- rV.pvals.mat1[,1]
rV.pvals.mat1 <- rV.pvals.mat1[,-1]

rV_pvals_matrix <- matrix(data = NA, nrow = dim(rV.pvals.mat1)[1], ncol = dim(rV.pvals.mat1)[2])
colnames(rV_pvals_matrix) <- colnames(rV.pvals.mat1)
rownames(rV_pvals_matrix) <- rownames(rV.pvals.mat1)

for (i in 1:dim(rV.pvals.mat1)[2]) {
  rV_pvals_matrix[,i] <- c(as.numeric(rV.pvals.mat1[[i]]))
}
diag(rV_pvals_matrix) <- 0
rV_pvals_matrix <- as.matrix(rV_pvals_matrix)

#### Genetic Contribution ####
rV_A.mat1 <- left_join(rV_A, rc.order, by=c("Phen2"="Region"))
rV_A.mat1$Phen2 <- rV_A.mat1$Region2
rV_A.mat1 <- 
  rV_A.mat1 %>% relocate(all_of(c("Phen2", list1)))
rV_A.mat1 <- 
  rV_A.mat1 %>% arrange(Order)
rV_A.mat1 <- 
  rV_A.mat1 %>% select(-Order, -Region2)
colnames(rV_A.mat1) <- c("Phen2", rV_A.mat1$Phen2)

row.names(rV_A.mat1) <- rV_A.mat1[,1]
rV_A.mat1 <- rV_A.mat1[,-1]

rA_matrix <- matrix(data = NA, nrow = dim(rV_A.mat1)[1], ncol = dim(rV_A.mat1)[2])
colnames(rA_matrix) <- colnames(rV_A.mat1)
rownames(rA_matrix) <- rownames(rV_A.mat1)

for (i in 1:dim(rV_A.mat1)[2]) {
  rA_matrix[,i] <- c(as.numeric(rV_A.mat1[[i]]))
}

diag(rA_matrix) <- 0

#### Genetic Contribution p-values ####
rV_A.pvals.mat1 <- left_join(rV_A_pvals, rc.order, by=c("Phen2"="Region"))
rV_A.pvals.mat1$Phen2 <- rV_A.pvals.mat1$Region2
rV_A.pvals.mat1 <- 
  rV_A.pvals.mat1 %>% relocate(all_of(c("Phen2", list1)))
rV_A.pvals.mat1 <- 
  rV_A.pvals.mat1 %>% arrange(Order)
rV_A.pvals.mat1 <- 
  rV_A.pvals.mat1 %>% select(-Order, -Region2)
colnames(rV_A.pvals.mat1) <- c("Phen2", rV_A.pvals.mat1$Phen2)

row.names(rV_A.pvals.mat1) <- rV_A.pvals.mat1[,1]
rV_A.pvals.mat1 <- rV_A.pvals.mat1[,-1]

rV_A_pvals_matrix <- matrix(data = NA, nrow = dim(rV_A.pvals.mat1)[1], ncol = dim(rV_A.pvals.mat1)[2])
colnames(rV_A_pvals_matrix) <- colnames(rV_A.pvals.mat1)
rownames(rV_A_pvals_matrix) <- rownames(rV_A.pvals.mat1)

for (i in 1:dim(rV_A.pvals.mat1)[2]) {
  rV_A_pvals_matrix[,i] <- c(as.numeric(rV_A.pvals.mat1[[i]]))
}
diag(rV_A_pvals_matrix) <- 0
rV_A_pvals_matrix <- as.matrix(rV_A_pvals_matrix)

#### Environmental Contribution ####
rV_E.mat1 <- left_join(rV_E, rc.order, by=c("Phen2"="Region"))
rV_E.mat1$Phen2 <- rV_E.mat1$Region2
rV_E.mat1 <- 
  rV_E.mat1 %>% relocate(all_of(c("Phen2", list1)))
rV_E.mat1 <- 
  rV_E.mat1 %>% arrange(Order)
rV_E.mat1 <- 
  rV_E.mat1 %>% select(-Order, -Region2)
colnames(rV_E.mat1) <- c("Phen2", rV_E.mat1$Phen2)

row.names(rV_E.mat1) <- rV_E.mat1[,1]
rV_E.mat1 <- rV_E.mat1[,-1]

rE_matrix <- matrix(data = NA, nrow = dim(rV_E.mat1)[1], ncol = dim(rV_E.mat1)[2])
colnames(rE_matrix) <- colnames(rV_E.mat1)
rownames(rE_matrix) <- rownames(rV_E.mat1)

for (i in 1:dim(rV_E.mat1)[2]) {
  rE_matrix[,i] <- c(as.numeric(rV_E.mat1[[i]]))
}

diag(rE_matrix) <- 0

#### Environmental Contribution p-values ####
rV_E.pvals.mat1 <- left_join(rV_E_pvals, rc.order, by=c("Phen2"="Region"))
rV_E.pvals.mat1$Phen2 <- rV_E.pvals.mat1$Region2
rV_E.pvals.mat1 <- 
  rV_E.pvals.mat1 %>% relocate(all_of(c("Phen2", list1)))
rV_E.pvals.mat1 <- 
  rV_E.pvals.mat1 %>% arrange(Order)
rV_E.pvals.mat1 <- 
  rV_E.pvals.mat1 %>% select(-Order, -Region2)
colnames(rV_E.pvals.mat1) <- c("Phen2", rV_E.pvals.mat1$Phen2)

row.names(rV_E.pvals.mat1) <- rV_E.pvals.mat1[,1]
rV_E.pvals.mat1 <- rV_E.pvals.mat1[,-1]

rV_E_pvals_matrix <- matrix(data = NA, nrow = dim(rV_E.pvals.mat1)[1], ncol = dim(rV_E.pvals.mat1)[2])
colnames(rV_E_pvals_matrix) <- colnames(rV_E.pvals.mat1)
rownames(rV_E_pvals_matrix) <- rownames(rV_E.pvals.mat1)

for (i in 1:dim(rV_E.pvals.mat1)[2]) {
  rV_E_pvals_matrix[,i] <- c(as.numeric(rV_E.pvals.mat1[[i]]))
}
diag(rV_E_pvals_matrix) <- 0
rV_E_pvals_matrix <- as.matrix(rV_E_pvals_matrix)

#### Final Matrices ####
rA_E_matrix <- rA_matrix
rA_E_matrix[upper.tri(rA_E_matrix)] <- rE_matrix[upper.tri(rE_matrix)]
rV_AE_pvals_matrix <- rV_A_pvals_matrix
rV_AE_pvals_matrix[upper.tri(rV_AE_pvals_matrix)] <- rV_E_pvals_matrix[upper.tri(rV_E_pvals_matrix)]

#### Plots ####
pdf(file = "Figure2_prelim.pdf", width = 20, height = 20)
par(mfrow=c(1,2))
corrplot(rV_matrix, type = "lower", method = "shade", col = viridis(n = 100, option = "D"), cl.pos = "n", cl.cex = 0.75, tl.cex = 1, diag = T, tl.col = 'black', pch.col = "black", pch.cex = 1, p.mat = rV_pvals_matrix) |> corrRect(c(1, 7, 15, 21, 23, 27))
corrplot(rA_E_matrix, method = "shade", col = viridis(n = 100, option = "D"), diag = F, tl.cex = 1, tl.col = 'black', cl.pos = "n", cl.cex = 0.75, pch.col = "black", pch.cex = 1, p.mat = rV_AE_pvals_matrix) |> corrRect(c(1, 7, 15, 21, 23, 27))
dev.off()

pdf(file = "Figure2_legend.pdf", width = 20, height = 20)
corrplot(rA_E_matrix, method = "shade", col = viridis(n = 100, option = "D"), diag = F, tl.cex = 1, tl.col = 'black', cl.pos = "b", cl.cex = 0.75, pch.col = "black", pch.cex = 1, p.mat = rV_AE_pvals_matrix) |> corrRect(c(1, 7, 15, 21, 23, 27))
dev.off()

# Manually combine the correlation plots & the colorbar