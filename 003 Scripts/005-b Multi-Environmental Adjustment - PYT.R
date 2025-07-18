###############################################################################
# Script: Calculation of BLUPs for different field PYT trials
# Author: Maria Guadalupe Montiel, Ph.D.
# Email: mmontiel@horizonseed.com
# Paper: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#        An Empirical Study in Rice
#
# Description:
# This script calculates the Best Linear Unbiased Predictions (BLUPs) for several traits
# based on field trial data. For each trait, a multienvironmental model (MET) is used,
# first adjusting for each planting date and then adjusting by year for the MP2 and MP6-8
# populations. The model is fitted using mixed-effects models (MERM) with fixed effects for 
# the trial, and random effects for the germplasmName. 
# The traits analyzed in this study are as follows:
# - Yield 
# - Days to Heading 
# - Whole Milling Percentage 
# - Chalk Impact 
###############################################################################
# Clear the workspace
rm(list = ls())

# Load required libraries
library(sommer)
library(SpATS)
library(dplyr) # For data manipulation
library(tidyr) # For reshaping data

##### Loading the data #####
# Load BLUP1 and add trial identifier
blup1 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_CN_2022_1_blues.csv", header = TRUE) %>%
  mutate(trial = "PYT CN 2022 1") %>%
  select(-germplasmName) %>% # Remove the original germplasmName column
  rename(germplasmName = adjusted_value.germplasmName) %>% # Rename adjusted_value.germplasmName to germplasmName
  

# Load BLUP2 and add trial identifier
blup2 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_CN_2022_2_blues.csv", header = TRUE) %>%
  mutate(trial = "PYT CN 2022 2") %>%
  select(-germplasmName) %>% # Remove the original germplasmName column
  rename(germplasmName = adjusted_value.germplasmName) %>% # Rename adjusted_value.germplasmName to germplasmName


##### Process the combined BLUP data #####
#Trial 1
data_processed_trial1 <- blup1 %>%
  # Remove unnecessary columns (customize as needed)
  select(-adjusted_value.colNumber, -adjusted_value.rowNumber, 
         -adjusted_value.colf, -adjusted_value.rowf) %>%
  # Pivot the data to wide format
  pivot_wider(
    names_from = trait,
    values_from = adjusted_value.predicted.values
  ) %>%
  # Group by germplasmName and summarize to remove duplicate entries
  group_by(germplasmName, trial) %>%
  summarize(
    Yield.LSU_01.0000138 = first(na.omit(`Yield.LSU_01.0000138`))
    
  ) %>%
  ungroup()

#Trial 2
data_processed_trial2 <- blup2 %>%
  # Remove unnecessary columns (customize as needed)
  select(-adjusted_value.colNumber, -adjusted_value.rowNumber, 
         -adjusted_value.colf, -adjusted_value.rowf) %>%
  # Pivot the data to wide format
  pivot_wider(
    names_from = trait,
    values_from = adjusted_value.predicted.values
  ) %>%
  # Group by germplasmName and summarize to remove duplicate entries
  group_by(germplasmName, trial) %>%
  summarize(
    Yield.LSU_01.0000138 = first(na.omit(`Yield.LSU_01.0000138`))
    
  ) %>%
  ungroup()

#data_processed_trial2$Chalk.impact.LSU_01.0000108 <- NA

# View the processed data
head(data_processed_trial1)
head(data_processed_trial2)

# Combine BLUEs data
#blup <- rbind (data_processed_trial1, data_processed_trial2)
blup <- rbind (blup1, blup2)
str(blup)

#### Adjust for planting dates within each year ============================#

# Fit model and extract BLUPs for Yield
met_yield <- mmer(fixed = Yield_Adjusted_Value ~ 1 + trial,  
                  random = ~ germplasmName,
                  rcov = ~ vsr(units),
                  data = blup)

blups_yield <- predict.mmer(met_yield, classify = c("germplasmName"), D = "germplasmName")
yield_blups <- blups_yield$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Yield_Adjusted_Value = predicted.value)


# Check the final result structure
print(head(yield_blups))

# Save the results to a CSV file
write.csv(yield_blups, "001 Phenotypic Data/BLUEs PYT/PYT_CN_2022_blups.csv", row.names = FALSE)

# Join all

CL_2022 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_CL_2022_blups.csv", header = TRUE) %>%
  mutate(trial = "PYT CL 2022")
CN_2022 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_CN_2022_blups.csv", header = TRUE) %>%
  mutate(trial = "PYT CN 2022")
PV_2022 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_PV_2022_blups.csv", header = TRUE) %>%
  mutate(trial = "PYT PV 2022")

blups_2022 <- rbind(CL_2022, CN_2022, PV_2022)
write.csv(blups_2022, "001 Phenotypic Data/BLUEs PYT/PYT_22_blups.csv")


CL_2021 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_CL_2021_blues.csv", header = TRUE) %>%
  mutate(trial = "PYT CL 2021")
CN_2021 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_CN_2021_blues.csv", header = TRUE) %>%
  mutate(trial = "PYT CN 2021")
PV_2021 <- read.csv("001 Phenotypic Data/BLUEs PYT/PYT_PV_2021_blues.csv", header = TRUE) %>%
  mutate(trial = "PYT PV 2021")

blups_2021 <- rbind(CL_2021, CN_2021, PV_2021)
write.csv(blups_2022, "001 Phenotypic Data/BLUEs PYT/PYT_21_blups.csv")

Blups_22_21 <- rbind(blups_2021, blups_2022)
write.csv(Blups_22_21, "001 Phenotypic Data/BLUEs PYT/PYT_21-22_blups.csv")
