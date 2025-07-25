

Dataset: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis: An Empirical Study in Rice
Author: Maria Montiel
Years: 2020–2022
Last Updated: May 2025

================================================================================
Description
================================================================================
This repository contains all raw data, processed datasets, and analytical R scripts used in the empirical study titled: 
“The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis in Rice.”

The study investigates how haplotype block size impacts genomic selection (GS) accuracy and the detection of epistatic interactions in rice breeding populations. 
It includes phenotypic and genotypic datasets across three years and utilizes various statistical models and genomic prediction pipelines.

================================================================================
Repository Structure
================================================================================
The repository is organized into three main folders:

├── 001 Phenotypic Data/
│   ├── MP materials/
│   └── Preliminary Yield Trials/
├── 002 Genotypic Data/
├── 003 Scripts/
└── README.txt

================================================================================
001 Phenotypic Data
================================================================================
Contains raw phenotypic data from:
- Multi-parent populations (MP6-8)
- Bi-parental population (MP2)
- Preliminary Yield Trials (PYT) evaluated in 2020, 2021, and 2022.

Subfolders:
- MP materials/
- Preliminary Yield Trials/

Traits included:
- Grain Yield
- Days to Heading
- Chalk
- Whole Milling Rate

Example file: project_20mp2.csv

================================================================================
002 Genotypic Data
================================================================================
Genotypic datasets in dosage matrix format (0, 1, 2) for:
- MP2
- MP4
- MP6-8
- Preliminary Yield Trials (2020–2022)

File format: .csv

Example file: MP2_Numeric_wright.csv

================================================================================
003 Scripts
================================================================================
This folder contains all R scripts used in the study for:
- Data processing
- Spatial adjustment
- Genomic prediction
- Epistasis modeling
- Linkage Disequilibrium analysis

List of R scripts and descriptions:

001 A mat and E mat.R
Calculates the additive (A) and epistatic (E) genomic relationship matrices for MP2 and MP6-8 rice populations using the sommer package. Inputs are SNP dosage matrices (coded as 0, 1, 2) located in the 002 Genotypic Data/ folder. SNP matrices are centered by subtracting 1, and markers with minor allele frequency (MAF) < 0.05 are filtered out. Uses A.mat() to compute the additive relationship matrix and E.mat() for the epistatic matrix. Outputs are saved as .rds files in the 002 Genotypic Data/Matrix Calculations/ folder: MP2_A_mat.rds and MP2_E_mat.rds for the MP2 population MP6-8_A_mat_sommer.rds and MP6-8_E_mat_sommer.rds for the MP6-8 population
These matrices are required for downstream genomic prediction and variance component analyses.

002-a Spatial adjustment and BLUES MP.R
- Performs spatial adjustment and calculates BLUES for MP materials using the 'sommer' package.

002-b Spatial adjustment and BLUES PYT.R
- Performs spatial adjustment and calculates BLUES for Preliminary Yield Trial materials.

(Note: script 002-b for PYT materials follows the same structure and logic.)
Performs spatial adjustment and calculates Best Linear Unbiased Estimates (BLUES) for Multi-parent (MP2 and MP6-8) rice populations. Uses the SpATS package for spatial modeling, applying a PSANOVA surface to correct for spatial variability across rows and columns in the field design. Processes raw phenotypic datasets from multiple years and trials (project_20mp2.csv, project_21mp2_1.csv, etc.).
Traits analyzed include: Yield, Days to Heading, Whole Milling Percentage, Chalk Impact.
Steps performed for each trial: Detects and replaces outliers in each trait using linear models. Converts row, column, and replicate into factors for modeling. Fits a spatial model with genotype as a random (MP2) or fixed (MP6-8) effect. Extracts BLUES and associated weights. Saves adjusted results to .csv files named per trial (e.g., mp2_2020_blues.csv, mp68_2022_2_blues.csv). These outputs are used in heritability estimation and genomic prediction steps.

003 Heritability and Variance Components by Trial.R
This script calculates narrow-sense heritability (h²) and variance components for multiple traits using two genetic models:
Model 1 (Additive only): Uses an additive relationship matrix (A matrix).
Model 2 (Additive + Epistatic): Uses both additive (A) and epistatic (A×A) relationship matrices.
Key Steps: Loads phenotypic data for a single trial. Loads additive (Ga) and epistatic (Gaa) genomic relationship matrices. Performs spatially-adjusted mixed model analysis for each trait. Extracts variance components and calculates h² under both genetic models.
Traits Analyzed: Yield, Days to Heading, Whole Milling Percentage, Chalk Impact

004 K-Fold Cross Validation with A and A+I Models.R
This script performs K-fold cross-validation (CV) to evaluate the prediction accuracy of two genomic selection models across multiple rice trials. The models assessed are:
Additive Model (A): Uses only the genomic relationship matrix (Ga) to model additive genetic effects.
Additive + Epistatic Model (A+I): Includes both the additive relationship matrix (Ga) and the epistatic relationship matrix (Gaa) to capture gene-gene interactions.

Data Preparation: Phenotypic data from rice trials is loaded, cleaned, and reshaped to focus on traits like grain yield, days to heading, and milling percentage. The dataset is matched with genotype data, which includes both additive and epistatic relationship matrices for genomic prediction.

K-Fold Cross-Validation: CV is performed 25 times (cyc = 25) using 5-fold validation (k = 5). For each trial and trait, prediction accuracy is calculated for both model.
Prediction accuracies for each model and trial are computed and saved.
Results include trait-specific prediction accuracies across different materials and models.

005 Multi-Environmental Adjustment.R
- Adjusts phenotypic data across environments to account for genotype-by-environment interactions.

006 Predictions MP and PYT Genomic Selection.R
- Performs genomic predictions for MP and PYT populations. Set up a mixed model to estimate genomic breeding values (GEBVs).
- Fit the model using the mmer function from the sommer package.
- Extract Genomic Breeding Values (GEBVs):
- Extract the GEBVs for each germplasm from the fitted model.
- Merge GEBVs with Phenotypic Data: Combine the GEBVs with the phenotypic data (BLUPs) for the validation set (MP2). Calculate Prediction Accuracy: Calculate the correlation between the predicted GEBVs and the observed phenotypic values to determine the accuracy of the genomic predictions.

007 LD Calculation.R
- Calculates Linkage Disequilibrium (LD) and generates supporting plots.

================================================================================
Abbreviations
================================================================================
MP2   : Bi-parental population 2  
MP6-8 : Multi-parent populations 6 through 8  
PYT   : Preliminary Yield Trials  
BLUES: Best Linear Unbiased Estimates  
LD    : Linkage Disequilibrium  
A Model   : Additive genetic model  
A+I Model : Additive + Epistatic interaction model  
K-Fold    : Cross-validation method

================================================================================
Data Formats
================================================================================
Phenotypic Data : .csv  
Genotypic Data  : .csv  
Scripts         : .R  

Steps to Reproduce
## 1. **Set Up the Environment**
   - Ensure you have R installed (version 4.0 or higher).
   - Install necessary R packages:
     ```r
     install.packages(c("sommer", "ggplot2", "SpATS"))
     ```

## 2. **Prepare the Data**
   - The phenotypic and genotypic data are stored in the following folders:
     - `001 Phenotypic Data/`: Contains phenotypic data from multi-parent populations (MP2, MP6-8) and preliminary yield trials (PYT).
     - `002 Genotypic Data/`: Contains genotypic data in dosage matrix format (0, 1, 2) for each population.
   - Example files:
     - Phenotypic: `project_20mp2.csv`
     - Genotypic: `Genotype_raw_MP2.csv`

## 3. **Data Processing**
   - **Spatial Adjustment and BLUES Calculation**:
     Use the script `002-a Spatial adjustment and BLUES MP.R` to perform spatial adjustment and calculate Best Linear Unbiased Estimates (BLUES) for MP materials, and `002-b Spatial adjustment and BLUES PYT.R` for PYT materials.
     - These steps account for spatial variability and correct the phenotypic data for each trial.

## 4. **Calculate Genetic Relationship Matrices**
   - Run the script `001 A mat and E mat.R` to calculate the additive (A) and epistatic (E) genomic relationship matrices for MP2 and MP6-8 populations.
     - The `A.mat()` function computes the additive relationship matrix.
     - The `E.mat()` function computes the epistatic relationship matrix.
     - These matrices are stored as `.rds` files in `002 Genotypic Data/Matrix Calculations/`.

## 5. **Estimate Heritability and Variance Components**
   - Use the script `003 Heritability and Variance Components by Trial.R` to calculate narrow-sense heritability (h²) and variance components for multiple traits (e.g., grain yield, days to heading) using both additive and additive + epistatic models.

## 6. **Perform K-Fold Cross-Validation**
   - The script `004 K-Fold Cross Validation with A and A+I Models.R` performs K-fold cross-validation (with 5 folds) to evaluate the prediction accuracy of the two genomic selection models:
     - **A Model**: Additive model using the additive genetic relationship matrix.
     - **A+I Model**: Additive + epistatic model.
   - This script computes prediction accuracy for different trials and traits, and results are saved for analysis.

## 7. **Genomic Prediction**
   - Run the script `006 Predictions MP and PYT Genomic Selection.R` to perform genomic predictions:
     - The script uses a mixed model to estimate Genomic Breeding Values (GEBVs).
     - The `mmer` function from the `sommer` package is used to fit the model.
     - GEBVs are extracted and merged with the phenotypic data (BLUPs) for the validation set (MP2).
     - The accuracy of predictions is assessed by calculating the correlation between predicted GEBVs and observed phenotypic values.

## 8. **Calculate Linkage Disequilibrium (LD)**
   - Use the script `007 LD Calculation.R` to calculate and plot linkage disequilibrium (LD) across the markers.

## 9. **Output and Results**
   - The outputs of the scripts include adjusted phenotypic data (e.g., BLUES), genomic relationship matrices, heritability estimates, prediction accuracies, and LD calculations.
   - All results are saved in respective `.csv` and `.rds` files within the repository.

By following these steps, you can reproduce the analysis and genomic predictions used in this study.

================================================================================
Citation
================================================================================
If you use this dataset or scripts in your research, please cite:

Montiel, M. (2025). 
The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis: An Empirical Study in Rice.

