###############################################################################
# Script: Heritability and Variance Components Calculation - By trial and trait
# Author: Maria Guadalupe Montiel, Ph.D.
# Email: mmontiel@horizonseed.com
# Paper: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#        An Empirical Study in Rice
#
# Description: 
# This script calculates narrow-sense heritability and variance components 
# using additive (A) and additive + epistatic (A+I) genetic models. 
###############################################################################

rm(list=ls())
library(sommer)
library(SpATS)
library(lme4)
library(car)


### Heritability for each model

pheno.v0 <- read.csv("001 Phenotypic Data/MP materials/project_20mp2.csv") # Replace for each Trial

# Genotype -------------------------------------------------------------------

#Ga <- as.matrix(readRDS("002 Genotypic Data/Clean Genotype/Ga_mp2_matrix.rds", row.names = 1))
#Ga <- readRDS("002 Genotypic Data/Matrix Calculations/MP6-8_A_mat_sommer.rds")
Geno <- as.matrix(read.csv("002 Genotypic Data/MP2_Numeric_wright.csv", row.names = 1))
#Gaa <- readRDS("002 Genotypic Data/Matrix Calculations/MP6-8_E_mat_sommer.rds") 

str(Geno)
Geno <- Geno-1
Ga <- A.mat(X=as.matrix(Geno))  
Ga[1:5, 1:5]

Gaa <- E.mat(X=as.matrix(Geno))  
Gaa[1:5, 1:5]


traits <- c("Yield.LSU_01.0000138", "Days.to.heading.LSU_01.0000130", 
            "Whole.milling.percentage.LSU_01.0000100", "Chalk.impact.LSU_01.0000108")

#geno.v <- Geno[rownames(Geno) %in% pheno.v0$germplasmName, ]

#### Model 1: additive effects - Calculating h2 and Gblups####------------------
pheno.v1 <- pheno.v0[pheno.v0$germplasmName %in% rownames(Ga), ]
pheno.v1$germplasmName2 <- pheno.v1$germplasmName # epistasis accounted with a cloned id factor

pheno.v1$rowf <- as.factor(pheno.v1$rowNumber)
pheno.v1$colf <- as.factor(pheno.v1$colNumber)
pheno.v1$repf <- as.factor(pheno.v1$replicate)

## Grain Yield
mm.sp.a <- mmer(fixed = Yield.LSU_01.0000138 ~ 1, 
                random = ~vsr(germplasmName,Gu=Ga) +
                  vsr(rowf) + vsr(colf) 
                #+ vsr (repf) # in case needed
                + spl2Da(rowNumber,colNumber),
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.a)$varcomp
vpredict(mm.sp.a, h2 ~ (V1) / ( V1+V5)) # h2 = (Va/Va+Ve)
mm.sp.a$U$`u:germplasmName`
gebv.blup <- merge(pheno.v0, 
                   data.frame(germplasmName=names(mm.gp$U$`u:germplasmName`$Days_Adjusted_Value),
                              gebv=as.numeric(mm.gp$U$`u:germplasmName`$Days_Adjusted_Value)),
                   by = "germplasmName")


## Days to Heading
mm.sp.a <- mmer(fixed = Days.to.heading.LSU_01.0000130 ~ 1,
                random = ~vsr(germplasmName,Gu=Ga) +
                  vsr(rowf) + vsr(colf) 
                #+ vsr (repf) # in case needed
                + spl2Da(rowNumber,colNumber),
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.a)$varcomp
vpredict(mm.sp.a, h2 ~ (V1) / ( V1+V5)) # h2 = (Va/Va+Ve)

## Chalk %
mm.sp.a <- mmer(fixed = Chalk.impact.LSU_01.0000108 ~ 1,
                random = ~vsr(germplasmName,Gu=Ga) +
                  vsr(rowf) + vsr(colf) 
               # + vsr (repf) # in case needed
                + spl2Da(rowNumber,colNumber),
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.a)$varcomp
vpredict(mm.sp.a, h2 ~ (V1) / ( V1+V5)) # h2 = (Va/Va+Ve)

## Whole Milling %
mm.sp.a <- mmer(fixed = Whole.milling.percentage.LSU_01.0000100 ~ 1,
                random = ~vsr(germplasmName,Gu=Ga) +
                  vsr(rowf) + vsr(colf) 
                #+ vsr (repf) # in case needed
                + spl2Da(rowNumber,colNumber),
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.a)$varcomp
vpredict(mm.sp.a, h2 ~ (V1) / ( V1+V5)) # h2 = (Va/Va+Ve)

#### Model 2: epistatic and additive  effects ####------------------------------

## Grain Yield
mm.sp.i <- mmer(fixed = Yield.LSU_01.0000138 ~ 1, 
                random = ~vsr(germplasmName,Gu=Ga)+ vsr(germplasmName2, Gu=Gaa)+
                  vsr(rowf) + vsr(colf) 
                #+ vsr (repf) 
                + spl2Da(rowNumber,colNumber), 
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.i)$varcomp
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V6)) # h2 =  (Va/Va+Ve)
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V2+V6)) # h2 = (Va/Va+VI+Ve)
#vpredict(mm.sp.i, h2 ~ (V1+V2) / ( V1+V2+V6)) # h2 = (Va+VI/Va+VI+Ve)

## Days to Heading
mm.sp.i <- mmer(fixed = Days.to.heading.LSU_01.0000130 ~ 1, 
                random = ~vsr(germplasmName,Gu=Ga)+ vsr(germplasmName2, Gu=Gaa)+
                  vsr(rowf) + vsr(colf) 
               # + vsr (repf) 
                + spl2Da(rowNumber,colNumber), 
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.i)$varcomp
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V6)) # h2 =  (Va/Va+Ve)
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V2+V6)) # h2 = (Va/Va+VI+Ve)
#vpredict(mm.sp.i, h2 ~ (V1+V2) / ( V1+V2+V6)) # h2 = (Va+VI/Va+VI+Ve)

## Chalk %
mm.sp.i <- mmer(fixed = Chalk.impact.LSU_01.0000108 ~ 1, 
                random = ~vsr(germplasmName,Gu=Ga)+ vsr(germplasmName2, Gu=Gaa)+
                  vsr(rowf) + vsr(colf) 
                #+ vsr (repf) 
                + spl2Da(rowNumber,colNumber), 
                rcov = ~ vsr(units),
                data = pheno.v1
)

summary(mm.sp.i)$varcomp
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V6)) # h2 =  (Va/Va+Ve)
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V2+V6)) # h2 = (Va/Va+VI+Ve)
#vpredict(mm.sp.i, h2 ~ (V1+V2) / ( V1+V2+V6)) # h2 = (Va+VI/Va+VI+Ve)

## Whole Milling %
mm.sp.i <- mmer(fixed =Whole.milling.percentage.LSU_01.0000100~ 1, 
                random = ~vsr(germplasmName,Gu=Ga)+ vsr(germplasmName2, Gu=Gaa)+
                  vsr(rowf) + vsr(colf) 
                #+ vsr (repf) 
                + spl2Da(rowNumber,colNumber), 
                rcov = ~ vsr(units),
                data = pheno.v1
)
summary(mm.sp.i)$varcomp
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V6)) # h2 =  (Va/Va+Ve)
vpredict(mm.sp.i, h2 ~ (V1) / ( V1+V2+V6)) # h2 = (Va/Va+VI+Ve)
#vpredict(mm.sp.i, h2 ~ (V1+V2) / ( V1+V2+V6)) # h2 = (Va+VI/Va+VI+Ve)


## ==========================The End ==========================================## 

