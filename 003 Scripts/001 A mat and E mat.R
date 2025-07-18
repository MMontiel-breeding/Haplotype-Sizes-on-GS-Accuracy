##############################################################################
# Script:    001 A mat and E mat.R
# Author:    Maria Guadalupe Montiel, Ph.D.
# Email:     mmontiel@horizonseed.com
# Paper:     The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#            An Empirical Study in Rice
#
# Description:
# This script calculates the additive (A) and epistatic (E) genomic relationship
# matrices for MP2 and MP6-8 populations using the sommer package. It reads in 
# SNP dosage matrices (coded as 0, 1, 2), filters markers with a minor allele 
# frequency (MAF) < 0.05, and generates both A and E matrices. The resulting 
# matrices are saved as RDS files for downstream genomic prediction analyses.
##############################################################################
rm(list=ls())
library(sommer)

#### MP2 Matrix Calculations ---------------------------------------------------

SNP_matrix <- as.matrix(read.csv("002 Genotypic Data/MP2_Numeric_wright.csv", header = TRUE, row.names=1))
SNP_Matrix<-as.numeric(SNP_matrix)
str(SNP_matrix)
SNP_matrix <- SNP_matrix - 1

# A_Mat by Sommer package
G_mat <- A.mat(X=as.matrix(SNP_matrix), 
               min.MAF = .05
               )

G_mat[1:5, 1:5]

E_mat <- E.mat(X=as.matrix(SNP_matrix),min.MAF = .05)
E_mat [1:5, 1:5]

saveRDS(G_mat, "002 Genotypic Data/Matrix Calculations/MP2_A_mat.rds")
saveRDS(E_mat, "002 Genotypic Data/Matrix Calculations/MP2_E_mat.rds")

#### MP6-8 Matrix Calculations -------------------------------------------------

SNP_matrix <- as.matrix(read.csv("002 Genotypic Data/004 Jomar's Data Processed/Agriplex MP6-8_Numeric_wright.csv", header = TRUE, row.names=1))
SNP_Matrix<-as.numeric(SNP_matrix)
str(SNP_matrix)
SNP_matrix <- SNP_matrix - 1

# A_Mat by Sommer package
G_mat <- A.mat(X=as.matrix(SNP_matrix), min.MAF = .05)
G_mat[1:5, 1:5]
saveRDS(G_mat, "002 Genotypic Data/Matrix Calculations/MP6-8_A_mat_sommer.rds")

# E mat by sommer package
E_mat <- E.mat(X=as.matrix(SNP_matrix),min.MAF = .05)
E_mat [1:5, 1:5]
saveRDS(E_mat, "002 Genotypic Data/Matrix Calculations/MP6-8_E_mat_sommer.rds")

#######------- The end-------------------------------------------------------###
