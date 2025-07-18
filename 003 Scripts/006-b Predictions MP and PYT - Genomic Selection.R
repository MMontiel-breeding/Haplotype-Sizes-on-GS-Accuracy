###############################################################################
# Script: Genomic Predictions and Accuracy Calculation using MP2 and MP6-8
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
library(snpReady)
library(dplyr)
library(tidyr)

#### Training set
#Combined
pheno.tr <- read.csv("001 Phenotypic Data/BLUEs MP/2022_MP68_blups.csv", header = T)[, c("germplasmName", "Yield_Adjusted_Value")]

#By Year
pheno.tr <- read.csv("001 Phenotypic Data/BLUEs MP/mp2_2021_1_blues.csv ", header = T) #[, c("germplasmName", "Yield_Adjusted_Value")]

# Process the data
pheno.tr <- pheno.tr %>%
  # Remove unnecessary columns if they exist, including "germplasmName"
  select(-any_of(c("germplasmName", "adjusted_value.colNumber", "adjusted_value.rowNumber", 
                   "adjusted_value.colf", "adjusted_value.rowf")), everything()) %>%
  # Pivot to wide format
  pivot_wider(
    names_from = trait,
    values_from = adjusted_value.predicted.values
  ) %>%
  # Group by adjusted_value.germplasmName and summarize to remove duplicate entries
  group_by(adjusted_value.germplasmName) %>%
  summarize(
    Yield.LSU_01.0000138 = first(na.omit(`Yield.LSU_01.0000138`))
  ) %>%
  ungroup()

# Rename columns: change adjusted_value.germplasmName to germplasmName
pheno.tr <- pheno.tr %>%
  rename(germplasmName = adjusted_value.germplasmName,
         Yield_Adjusted_Value = `Yield.LSU_01.0000138`)

# Print the first few rows to verify changes
head(pheno.tr)

#### Validation Set
pyt1 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_21_blups.csv", header = T)[, c("germplasmName", "Yield_Adjusted_Value")]
pyt2 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_22_blups.csv", header = T)[, c("germplasmName", "Yield_Adjusted_Value")]
pheno.v0 <- rbind(pyt21, pyt22)

#for combined
pheno.v0 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_22_blups.csv", header = T)[, c("germplasmName", "Yield_Adjusted_Value")]
pheno.tr2 <- pheno.tr[!pheno.tr$germplasmName %in% pheno.v0$germplasmName,]

### Gebvs for blups.sp (yield)

pheno.v0$Yield_Adjusted_Value=NA
pheno.all <- rbind(pheno.tr2,pheno.v0)

#### Genotype and GM
geno_mp <- read.csv("002 Genotypic Data/004 Data Processed/Agriplex MP6-8_Numeric_notimputed.csv", row.names = 1)
#geno_py <- read.csv(unique("002 Genotypic Data/004 Data Processed/Agriplex PYT 550_Numeric_wright.csv", row.names = 1))
geno_py <- read.csv("002 Genotypic Data/004 Data Processed/Agriplex PYT 550_Numeric_notimputed.csv")
geno_py[, 1] <- make.unique(as.character(geno_py[, 1]))
rownames(geno_py) <- geno_py[, 1]
geno_py <- geno_py[, -1]

# choosing common markers
common_markers <- intersect(colnames(geno_mp), colnames(geno_py))
geno_mp <- geno_mp[, common_markers]
geno_py <- geno_py[, common_markers]
geno.all <- rbind(geno_mp, geno_py)
geno <- geno.all
geno <- geno.all[rownames(geno.all) %in% pheno.all$germplasmName, ]
geno[geno == 99] <- NA

geno <- geno - 1
geno[1:5, 1:5]

G_mat <- A.mat(X=as.matrix(geno), min.MAF = .05)
G_mat[1:5, 1:5]

# Mixed Model 
mm.gp <- mmer(fixed = Yield_Adjusted_Value ~ 1, 
              random = ~vsr(germplasmName,Gu=G_mat),
              rcov = ~ vsr(units),
              data = droplevels(pheno.all[pheno.all$germplasmName %in% rownames(G_mat),])
)

#pheno.v0 <- rbind(pyt21, pyt22)
pheno.v0 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_22_blups.csv",header = T)[, c("germplasmName", "Yield_Adjusted_Value")]
gebv.blup <- merge(pheno.v0, 
                   data.frame(germplasmName=names(mm.gp$U$`u:germplasmName`$Yield_Adjusted_Value),
                              gebv=as.numeric(mm.gp$U$`u:germplasmName`$Yield_Adjusted_Value)),
                   by = "germplasmName")


#####  Calculating the prediction accuracy
cor(gebv.blup$gebv, gebv.blup$Yield_Adjusted_Value)
