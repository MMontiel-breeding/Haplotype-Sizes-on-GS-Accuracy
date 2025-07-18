###############################################################################
# Script: Calculation of BLUPs for different field trials
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
blup1 <- read.csv("../BLUEs/mp68_2022_1_blues.csv", header = TRUE) %>%
  select(-germplasmName) %>% # Remove the original germplasmName column
  rename(germplasmName = adjusted_value.germplasmName) %>% # Rename adjusted_value.germplasmName to germplasmName
  mutate(trial = "MP68 2022 1")

# Load BLUP2 and add trial identifier
blup2 <- read.csv("../BLUEs/mp68_2022_2_blues.csv", header = TRUE) %>%
  select(-germplasmName) %>% # Remove the original germplasmName column
  rename(germplasmName = adjusted_value.germplasmName) %>% # Rename adjusted_value.germplasmName to germplasmName
  mutate(trial = "MP68 2022 2")

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
    Yield.LSU_01.0000138 = first(na.omit(`Yield.LSU_01.0000138`)),
    Days.to.heading.LSU_01.0000130 = first(na.omit(`Days.to.heading.LSU_01.0000130`)),
    Whole.milling.percentage.LSU_01.0000100 = first(na.omit(`Whole.milling.percentage.LSU_01.0000100`))
    # Uncomment if needed
    , Chalk.impact.LSU_01.0000108 = first(na.omit(`Chalk.impact.LSU_01.0000108`))
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
    Yield.LSU_01.0000138 = first(na.omit(`Yield.LSU_01.0000138`)),
    Days.to.heading.LSU_01.0000130 = first(na.omit(`Days.to.heading.LSU_01.0000130`)),
    Whole.milling.percentage.LSU_01.0000100 = first(na.omit(`Whole.milling.percentage.LSU_01.0000100`))
    # Uncomment if needed
    #, Chalk.impact.LSU_01.0000108 = first(na.omit(`Chalk.impact.LSU_01.0000108`))
  ) %>%
  ungroup()

#data_processed_trial2$Chalk.impact.LSU_01.0000108 <- NA

# View the processed data
head(data_processed_trial1)
head(data_processed_trial2)

# Combine BLUEs data
blup <- rbind (data_processed_trial1, data_processed_trial2)
str(blup)

#### Adjust for planting dates within each year ============================#

# Fit model and extract BLUPs for Yield
met_yield <- mmer(fixed = Yield.LSU_01.0000138 ~ 1 + trial,  
                  random = ~ germplasmName,
                  rcov = ~ vsr(units),
                  data = blup)

blups_yield <- predict.mmer(met_yield, classify = c("germplasmName"), D = "germplasmName")
yield_blups <- blups_yield$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Yield_Adjusted_Value = predicted.value)

# Fit model and extract BLUPs for Days to Heading
met_days <- mmer(fixed = Days.to.heading.LSU_01.0000130 ~ 1 + trial,  
                 random = ~ germplasmName,
                 rcov = ~ vsr(units),
                 data = blup)

blups_days <- predict.mmer(met_days, classify = c("germplasmName"), D = "germplasmName")
days_blups <- blups_days$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Days_Adjusted_Value = predicted.value)

# Fit model and extract BLUPs for Whole Milling Percentage
met_milling <- mmer(fixed = Whole.milling.percentage.LSU_01.0000100 ~ 1 + trial,  
                    random = ~ germplasmName,
                    rcov = ~ vsr(units),
                    data = blup)

blups_milling <- predict.mmer(met_milling, classify = c("germplasmName"), D = "germplasmName")
milling_blups <- blups_milling$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Milling_Adjusted_Value = predicted.value)

# Fit model and extract BLUPs for Chalk Impact
met_chalk <- mmer(fixed = Chalk.impact.LSU_01.0000108 ~ 1 + trial,  
                  random = ~ germplasmName,
                  rcov = ~ vsr(units),
                  data = blup)

blups_chalk <- predict.mmer(met_chalk, classify = c("germplasmName"), D = "germplasmName")
chalk_blups <- blups_chalk$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Chalk_Adjusted_Value = predicted.value)

# Merge all results by germplasmName
final_results <- yield_blups %>%
  left_join(days_blups, by = "germplasmName") %>%
  left_join(milling_blups, by = "germplasmName") #%>%
  left_join(chalk_blups, by = "germplasmName")

# Check the final result structure
print(head(final_results))

# Save the results to a CSV file
write.csv(final_results, "../BLUEs/2022_MP68_blups.csv", row.names = FALSE)

#### Adjust for year effects having one value per line ======================#
# Clear the workspace
rm(list = ls())
##### Loading the data #####
# Load the dataset (adjust the path to your file)
blup1 <- read.csv("../BLUEs/mp2_2020_blues.csv", header = TRUE)
blup2 <- read.csv("../BLUEs/2021_MP2_blups.csv", header = TRUE)

# Add trial year information based on your filenames
blup1 <- blup1 %>%
  mutate(Year = as.character(2020), trial = "MP6-8 2021")

blup2 <- blup2 %>%
  mutate(Year = as.character(2021), trial = "MP6-8 2022")

# Combine both datasets into one
blup <- bind_rows(blup1, blup2)

# Check the structure of the combined dataset
head(blup)

# Perform MET for Yield
met_yield <- mmer(fixed = Yield_Adjusted_Value ~ 1 + Year,  
                  random = ~ germplasmName,
                  rcov = ~ vsr(units),
                  data = blup)

blups_yield <- predict.mmer(met_yield, classify = c("germplasmName"), D = "germplasmName")
yield_blups <- blups_yield$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Yield_Adjusted_Value = predicted.value)

# Perform MET for Days to Heading
met_days <- mmer(fixed = Days_Adjusted_Value ~ 1 + Year,  
                 random = ~ germplasmName,
                 rcov = ~ vsr(units),
                 data = blup)

blups_days <- predict.mmer(met_days, classify = c("germplasmName"), D = "germplasmName")
days_blups <- blups_days$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Days_Adjusted_Value = predicted.value)

# Perform MET for Milling Percentage
met_milling <- mmer(fixed = Milling_Adjusted_Value ~ 1 + Year,  
                    random = ~ germplasmName,
                    rcov = ~ vsr(units),
                    data = blup)

blups_milling <- predict.mmer(met_milling, classify = c("germplasmName"), D = "germplasmName")
milling_blups <- blups_milling$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Milling_Adjusted_Value = predicted.value)

# Perform MET for Chalk Percentage
met_chalk <- mmer(fixed = Chalk_Adjusted_Value ~ 1 + Year,  
          random = ~ germplasmName,         
          rcov = ~ vsr(units),              
          data = blup)    

blups_chalk <- predict.mmer(met_chalk, classify = c("germplasmName"), D = "germplasmName")
chalk_blups <- blups_chalk$pvals %>%
  select(germplasmName, predicted.value) %>%
  rename(Chalk_Adjusted_Value = predicted.value)

# Merge the BLUPs results by germplasmName, including Chalk Impact
final_results <- yield_blups %>%
  left_join(days_blups, by = "germplasmName") %>%
  left_join(milling_blups, by = "germplasmName") %>%
  left_join(chalk_blups, by = "germplasmName")

# View the merged final results
head(final_results)

# Save the results to a CSV file
write.csv(final_results, "../BLUEs/20-21 MP2_blups.csv", row.names = FALSE)
