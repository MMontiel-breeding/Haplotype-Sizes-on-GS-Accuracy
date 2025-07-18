###############################################################################
# Script: Genomic Predictions and Accurcay Caluclation using MP2 and MP6-8
# Training and Validation Sets
# Author: Maria Guadalupe Montiel, Ph.D.
# Email: mmontiel@horizonseed.com
# Paper: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#        An Empirical Study in Rice
# 
###############################################################################
# Clear workspace and load libraries
rm(list=ls())
library(sommer)
library(ggplot2)

#### Training set
pheno.tr <- read.csv("001 Phenotypic Data/BLUEs MP/21-22 MP6-8_blups.csv",header = T)

#### Validation Set
pheno.v0 <- read.csv("001 Phenotypic Data/BLUEs MP/20-21 MP2_blups.csv",header = T)
pheno.tr2 <- pheno.tr[!pheno.tr$germplasmName %in% pheno.v0$germplasmName,]
pheno.v0$germplasmName2 <- pheno.v0$germplasmName
pheno.tr2$germplasmName2 <- pheno.tr2$germplasmName

### Gebvs for blups.sp (yield)

pheno.v0$Chalk_Adjusted_Value=NA
pheno.all <- rbind(pheno.tr2,pheno.v0)

#### Genotype and GM
geno_mp2 <- as.matrix(read.csv("002 Genotypic Data/004 Data Processed/Agriplex MP2 550_Numeric_wright.csv", row.names = 1))
geno_mp68 <- as.matrix(read.csv("002 Genotypic Data/004 Data Processed/Agriplex MP6-8_Numeric_wright.csv", row.names = 1))

setdiff(colnames(geno_mp2), colnames(geno_mp68))  # Columns in MP2 but not in MP6-8
setdiff(colnames(geno_mp68), colnames(geno_mp2))  # Columns in MP6-8 but not in MP2
common_markers <- intersect(colnames(geno_mp2), colnames(geno_mp68))
geno_mp2 <- geno_mp2[, common_markers]
geno_mp68 <- geno_mp68[, common_markers]

geno <- rbind(geno_mp2, geno_mp68)
geno <- geno - 1
G_mat <- A.mat(X=as.matrix(geno), min.MAF = .05)

# Model using only additive effects
mm.gp <- mmer(fixed = Chalk_Adjusted_Value ~ 1, 
              random = ~vsr(germplasmName,Gu=G_mat),
              rcov = ~ vsr(units),
              data = droplevels(pheno.all[pheno.all$germplasmName %in% rownames(G_mat),])
)

pheno.v0 <- read.csv("001 Phenotypic Data/BLUEs MP/20-21 MP2_blups.csv",header = T)
gebv.blup <- merge(pheno.v0, 
                   data.frame(germplasmName=names(mm.gp$U$`u:germplasmName`$Chalk_Adjusted_Value),
                              gebv=as.numeric(mm.gp$U$`u:germplasmName`$Chalk_Adjusted_Value)),
                   by = "germplasmName")
gebv.blup <- na.omit(gebv.blup)

#####  Calculating the prediction accuracy
cor(gebv.blup$gebv, gebv.blup$Chalk_Adjusted_Value)
