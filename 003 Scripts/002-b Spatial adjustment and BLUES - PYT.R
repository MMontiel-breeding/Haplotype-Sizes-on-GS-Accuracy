###############################################################################
# Script:    002-b Spatial adjustment and BLUES PYT.R
# Author:    Maria Guadalupe Montiel, Ph.D.
# Email:     mmontiel@horizonseed.com
# Paper:     The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#            An Empirical Study in Rice
#
# Description:
# This script performs spatial adjustment and calculates Best Linear Unbiased 
# Estimates (BLUES) for traits measured in Preliminary Yield Trials (PYT) 
# across the 2020–2022 field seasons. Using the SpATS package (via sommer), 
# the script accounts for spatial trends within each trial to improve the 
# precision of phenotypic estimates. The workflow includes data cleaning, 
# model fitting using PSANOVA, and extraction of adjusted BLUES, which are 
# exported as CSV files for downstream genomic analysis.
#
# Traits analyzed:
# - Grain Yield
# - Days to Heading
# - Whole Milling Percentage
# - Chalk
###############################################################################

# Clear workspace and load libraries
rm(list=ls())

# Load necessary libraries
library(car)
library(SpATS)
library(dplyr)

# Load phenotype data for the PYT
PYT_CN_2020 <- read.csv("raw/Preliminary Yield Trials/PY_CONV_LONG_RRS_2020_phenotypes.csv", header = TRUE)
PYT_CL_2020 <- read.csv("raw/Preliminary Yield Trials/PY_CL_LONG_RRS_2020_phenotypes.csv", header = TRUE)
PYT_PV_2020 <- read.csv("raw/Preliminary Yield Trials/PY_PV_LONG_RRS_2020_phenotypes.csv", header = TRUE)

PYT_CN_2021 <- read.csv("raw/Preliminary Yield Trials/21PYL_phenotypes.csv", header = TRUE)
PYT_CL_2021_1 <- read.csv("raw/Preliminary Yield Trials/21CLPYL-T1_phenotypes.csv", header = TRUE)
PYT_CL_2021_2 <- read.csv("raw/Preliminary Yield Trials/21CLPYL-T2_phenotypes.csv", header = TRUE)
PYT_PV_2021 <- read.csv("raw/Preliminary Yield Trials/21PVPY_phenotypes.csv", header = TRUE)

PYT_CN_2022_1 <- read.csv("raw/Preliminary Yield Trials/22_PYL_RRS_phenotypes.csv", header = TRUE)
PYT_CN_2022_2 <- read.csv("raw/Preliminary Yield Trials/22_CLPYL_RRSL_phenotypes.csv", header = TRUE)
PYT_CL_2022_1 <- read.csv("raw/Preliminary Yield Trials/22_CLPYL_RRS_phenotypes.csv", header = TRUE)
PYT_CL_2022_2 <- read.csv("raw/Preliminary Yield Trials/22_CLPYL_RRSL_phenotypes.csv", header = TRUE)
PYT_PV_2022_1 <- read.csv("raw/Preliminary Yield Trials/22_PVPY_RRS_phenotypes.csv", header = TRUE)
PYT_PV_2022_2 <- read.csv("raw/Preliminary Yield Trials/22_PVPY_RRSL_phenotypes.csv", header = TRUE)

# Traits for analysis
trait_columns <- c(
  "Yield.LSU_01.0000138"
)

# Function to process each trial and save results
process_trial <- function(data, trial_name) 
  {
  # Outlier detection and replacement with NA
  #for (trait in trait_columns) {
   # fit <- lm(as.formula(paste(trait, "~ germplasmName")), data = data)
    #outlier <- names(outlierTest(fit)$p)
    #if (length(outlier) > 0) {
     # data[outlier, trait] <- NA
    #}
  #}
  
  # Prepare for spatial adjustment
  data <- data %>%
    mutate(
      rowf = as.factor(rowNumber),
      colf = as.factor(colNumber),
      repf = as.factor(replicate)
    )
  
  nrow <- max(data$rowNumber, na.rm = TRUE)
  ncol <- max(data$colNumber, na.rm = TRUE)
  
  # Number of segments for PSANOVA
  nseg.row <- nrow
  nseg.col <- ncol
  
  # Initialize a list to store adjusted results for each trait
  adjusted_results <- list()
  
  # Loop through each trait for spatial adjustment and BLUES extraction
  for (trait in trait_columns) {
    fitR <- SpATS(
      response = trait,
      fixed = ~ 1,
      random = ~ colf + rowf+ repf, # Add `repf` if needed
      spatial = ~ PSANOVA(colNumber, rowNumber, nseg = c(nseg.col, nseg.row)),
      genotype = "germplasmName",
      genotype.as.random = FALSE,
      data = data
    )
    
    # Extract variance-covariance matrix and weights
    vcov.mme <- fitR$vcov$C11_inv
    w <- diag(vcov.mme) # in case we want to see the weights, but we can ommit
    
    # Extract BLUES
    blues <- predict.SpATS(fitR, which = "germplasmName")
    adjusted <- data.frame(germplasmName = rownames(blues), adjusted_value = blues, weight = w)
    
    
    # Add trait name and year to results
    adjusted$trait <- trait
    adjusted$year <- data$year[1]  # Add year from the dataset
    
    adjusted_results[[trait]] <- adjusted
  }
  
  # Combine all adjusted results
  combined_results <- do.call(rbind, adjusted_results)
  
  # Save results to a CSV file
  output_file <- paste0(trial_name, "_blues.csv")
  write.csv(combined_results, output_file, row.names = FALSE)
  
  # Confirmation message
  cat("Processed", trial_name, "- Results saved to", output_file, "\n")
}

# Process each PYT trial
process_trial(PYT_CN_2020, "PYT_CN_2020")
process_trial(PYT_CL_2020, "PYT_CL_2020")
process_trial(PYT_PV_2020, "PYT_PV_2020")

process_trial(PYT_CN_2021, "PYT_CN_2021")
process_trial(PYT_CL_2021_1, "PYT_CL_2021_1")
process_trial(PYT_CL_2021_2, "PYT_CL_2021_2")
process_trial(PYT_PV_2021, "PYT_PV_2021")

process_trial(PYT_CN_2022_1, "PYT_CN_2022_1")
process_trial(PYT_CN_2022_2, "PYT_CN_2022_2")
process_trial(PYT_CL_2022_1, "PYT_CL_2022_1")
process_trial(PYT_CL_2022_2, "PYT_CL_2022_2")
process_trial(PYT_PV_2022_1, "PYT_PV_2022_1")
process_trial(PYT_PV_2022_2, "PYT_PV_2022_2")
##===========================================================================

PYT_2020_CN <- read.csv("PYT_CN_2020_blues.csv")
PYT_2020_CL <- read.csv("PYT_CL_2020_blues.csv")
PYT_2020_PV <- read.csv("PYT_PV_2020_blues.csv")

PYT_2020 <- rbind(PYT_2020_CN, PYT_2020_CL, PYT_2020_PV)
write.csv(PYT_2020, "2020 Preliminary Yield Trials (2).csv")
PYT_2020 <- read.csv("2020 Preliminary Yield Trials (2).csv")

PYT_2020_old <- read.csv("2020 Preliminary Yield Trials.csv")

cor(PYT_2020$adjusted_value.predicted.values, PYT_2020_old$blue)
