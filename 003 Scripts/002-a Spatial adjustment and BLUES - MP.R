###############################################################################
# Script: Calculation of BLUES and Spatial Adjustment using SpATS
# Author: Maria Guadalupe Montiel, Ph.D.
# Email: mmontiel@horizonseed.com
# Paper: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#        An Empirical Study in Rice
#
# Description:
# This script calculates Best Linear Unbiased Estimates (BLUEs) for various 
# traits by adjusting the means for spatial components using the SpATS package 
# from the sommer suite. The analysis includes outlier detection and 
# replacement, spatial adjustment using PSANOVA, and extraction of BLUES for 
# each trial and trait. The adjusted results are saved to CSV files for further 
# analysis.
#
# Traits analyzed:
# - Yield
# - Days to Heading
# - Whole Milling Percentage
# - Chalk Impact
###############################################################################

# Clear workspace and load libraries
rm(list=ls())

# Load necessary libraries
library(car)
library(SpATS)
library(dplyr)
##==========MP2 materials======================================================#

# Load phenotype data for the MP2 materials
mp2_2020 <- read.csv("raw/project_20mp2.csv", header = TRUE)
mp2_2021_1 <- read.csv("raw/project_21mp2_1.csv", header = TRUE)
mp2_2021_2 <- read.csv("raw/project_21mp2_2.csv", header = TRUE)

# Add a "year" column to each dataset
mp2_2020$year <- 2020
mp2_2021_1$year <- 2021
mp2_2021_2$year <- 2021

# Traits for analysis
trait_columns <- c(
  "Yield.LSU_01.0000138",
  "Days.to.heading.LSU_01.0000130",
  "Whole.milling.percentage.LSU_01.0000100",
  "Chalk.impact.LSU_01.0000108"
)

# Function to process each trial and save results
process_trial <- function(data, trial_name) {
  # Outlier detection and replacement with NA
  for (trait in trait_columns) {
    fit <- lm(as.formula(paste(trait, "~ germplasmName")), data = data)
    outlier <- names(outlierTest(fit)$p)
    if (length(outlier) > 0) {
      data[outlier, trait] <- NA
    }
  }
  
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
      random = ~ colf + rowf + repf, # Add `repf` if needed
      spatial = ~ PSANOVA(colNumber, rowNumber, nseg = c(nseg.col, nseg.row)),
      genotype = "germplasmName",
      genotype.as.random = TRUE,
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

# Process each trial
process_trial(mp2_2020, "mp2_2020")
process_trial(mp2_2021_1, "mp2_2021_1")
process_trial(mp2_2021_2, "mp2_2021_2")

##==================MP6-8 materilas ==========================================#

# Load phenotype data for the MP6-8 materials
mp68_2021_1 <- read.csv("raw/project_21mp68_1.csv", header = TRUE)
mp68_2021_2 <- read.csv("raw/project_21mp68_2.csv", header = TRUE)
mp68_2022_1 <- read.csv("raw/project_22mp68_1.csv", header = TRUE)
mp68_2022_2 <- read.csv("raw/project_22mp68_2.csv", header = TRUE)

# Add a "year" column to each dataset
mp68_2021_1$year <- 2021
mp68_2021_2$year <- 2021
mp68_2022_1$year <- 2022
mp68_2022_2$year <- 2022

# Traits for analysis
trait_columns <- c(
  "Yield.LSU_01.0000138",
  "Days.to.heading.LSU_01.0000130",
  "Whole.milling.percentage.LSU_01.0000100"#,
  #"Chalk.impact.LSU_01.0000108"
)

# Function to process each trial and save results
process_trial <- function(data, trial_name) {
  # Outlier detection and replacement with NA
  for (trait in trait_columns) {
    fit <- lm(as.formula(paste(trait, "~ germplasmName")), data = data)
    outlier <- names(outlierTest(fit)$p)
    if (length(outlier) > 0) {
      data[outlier, trait] <- NA
    }
  }
  
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
      random = ~ colf + rowf,
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

# Process each trial
process_trial(mp68_2021_1, "mp68_2021_1")
process_trial(mp68_2021_2, "mp68_2021_2")
process_trial(mp68_2022_1, "mp68_2022_1")
process_trial(mp68_2022_2, "mp68_2022_2")

### ---------------The end -------------------------------------------------#
