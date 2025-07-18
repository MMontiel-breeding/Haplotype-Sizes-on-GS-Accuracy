###############################################################################
# Script: K-Fold Cross Validation of Genomic Prediction Accuracy by Trial
# Models: Additive (A) Model and Additive + Epistatic (A+I) Model
# Author: Maria Guadalupe Montiel, Ph.D.
# Email: mmontiel@horizonseed.com
# Paper: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#        An Empirical Study in Rice
#
# Description:
# This script evaluates the predictive performance of genomic selection models 
# across multiple rice trials using K-fold cross-validation (CV). Two genomic 
# prediction models are compared:
#   1. **Additive Model (A)** – incorporates the genomic relationship matrix (G).
#   2. **Additive + Epistatic Model (A+I)** – includes both the additive (G) and 
#      epistatic (G×G) relationship matrices.
###############################################################################

# Clear workspace and load libraries
rm(list=ls())

## Reshaping the BLUEs data frames to proceed CV
# Load necessary library
library(dplyr)
library(tidyr)

# Load the data
data <- read.csv("001 Phenotypic Data/BLUEs MP/mp68_2022_2_blues.csv") # Replace here by the trial under consideration

# Process the data
data_processed <- data %>%
  # Remove the first column
  select(-germplasmName) %>%
  # Rename the second column
  rename(germplasmName = adjusted_value.germplasmName) %>%
  # Remove columns 3 to 6
  select(-adjusted_value.colNumber, -adjusted_value.rowNumber, 
         -adjusted_value.colf, -adjusted_value.rowf) %>%
  # Pivot the data to wide format
  pivot_wider(
    names_from = trait, 
    values_from = adjusted_value.predicted.values
  ) %>%
  # Group by germplasmName and summarize to remove duplicate entries
  group_by(germplasmName) %>%
  summarize(
    Yield.LSU_01.0000138 = first(na.omit(`Yield.LSU_01.0000138`)),
    Days.to.heading.LSU_01.0000130 = first(na.omit(`Days.to.heading.LSU_01.0000130`)),
    Whole.milling.percentage.LSU_01.0000100 = first(na.omit(`Whole.milling.percentage.LSU_01.0000100`))
    #,Chalk.impact.LSU_01.0000108 = first(na.omit(`Chalk.impact.LSU_01.0000108`))
  ) %>%
  ungroup()

# View the processed data
head(data_processed)

# Save the processed data to a new CSV file
write.csv(data_processed, "001 Phenotypic Data/BLUEs MP/processed_blues_data_mp68_2022_2.csv", row.names = FALSE, na = "")
pheno <- read.csv("001 Phenotypic Data/BLUEs MP/processed_blues_data_mp68_2022_2.csv")
str(pheno)

# Read Genetic Matrix
Ga <- as.matrix(readRDS("002 Genotypic Data/Matrix Calculations/MP6-8_A_mat_sommer.rds"))
Gaa <- as.matrix(readRDS("002 Genotypic Data/Matrix Calculations/MP6-8_E_mat_sommer.rds"))
pheno$germplasmName2 <- pheno$germplasmName # create a new column for the epistatic effects

# Checking Dimensions: 
common_names <- intersect(rownames(Ga), rownames(Gaa))
common_names <- intersect(common_names, pheno$germplasmName)

# Filter phenotype and genotype matrices
pheno <- pheno[pheno$germplasmName %in% common_names, ]
Ga <- Ga[common_names, common_names]
Gaa <- Gaa[common_names, common_names]

ID <- colnames(pheno)[1]
res1 <- data.frame()
cyc=25
k=5

# K- Fold CrossValidation using the A Model Only

for (r in 1:25){
  df. <- pheno
  df. <- df.[sample(1:nrow(df.),nrow(df.)),]
  G_mat_CV <- Ga
  G_mat.v <- Ga
  pred.loo <- data.frame()
  for (tr in c("Yield.LSU_01.0000138", "Days.to.heading.LSU_01.0000130", 
               "Whole.milling.percentage.LSU_01.0000100"
               , "Chalk.impact.LSU_01.0000108"
  )){
    acc.k=NULL
    flds <- caret::createFolds(seq(1:nrow(df.)),
                               k = k, list = TRUE, returnTrain = FALSE)
    pred.i <- data.frame()
    for(fd in 1:k){
      df <- df.
      pred.k <- data.frame()
      geno.k <- df[flds[[fd]],"germplasmName"]
      df[df[,ID]%in%geno.k,tr] <- NA    #NA masking of the lines in the k-fold to be predicted
      df[,ID] <- as.factor(df[,ID])
      frm.cv <- paste(tr,"~ 1  ",collapse = " ")
      fix_frm.cv <- formula(paste(frm.cv, collapse = " "))
      mm.gb <-mmer(fixed = fix_frm.cv,
                   random =  ~  vsr(germplasmName ,Gu=G_mat.v) ,
                   rcov = ~vsr(units),
                   #tolParInv = 100000,
                   data = droplevels(df[df$germplasmName %in% rownames(G_mat.v),])
      )
      pred <- data.frame(germplasmName=names(mm.gb$U[[1]][[tr ]]),
                         gebv=as.numeric(mm.gb$U[[1]][[tr ]]))
      
      rownames(pred) <- pred[,1]
      pred.k <- rbind(pred.k,pred[geno.k ,])
      pred <- merge(pred.k,df.[,c(ID,tr)],by=ID)
      acc.k <- c(acc.k,cor(pred[,2],pred[,3],use = "complete.obs"))
      pred.i <- rbind(pred.i,pred.k)
    }
    pred.i <- merge(pred.i,df.[,c(ID,tr)],by=ID)
    
    res1 <- rbind(res1,data.frame(
      trait=tr,n=nrow(df.),
      acc.k=mean(acc.k,na.rm = T),
      acc.i=cor(pred.i[,2],pred.i[,3],use = "complete.obs"))
    )
    
  }
}

# save results
write.csv(res1, "001 Phenotypic Data/CV results/MP6-8_2022_2_CV_A model.csv") # Replace by the name of the trial


# K- Fold Cross Validation using the A+I Model  - Include epistatic kernel

for (r in 1:25){
  df. <- pheno
  df. <- df.[sample(1:nrow(df.),nrow(df.)),]
  G_mat_CV <- Ga
  G_mat.v <- Ga
  pred.loo <- data.frame()
  for (tr in c("Yield.LSU_01.0000138", "Days.to.heading.LSU_01.0000130", 
               "Whole.milling.percentage.LSU_01.0000100"
               , "Chalk.impact.LSU_01.0000108"
  )){
    acc.k=NULL
    flds <- caret::createFolds(seq(1:nrow(df.)),
                               k = k, list = TRUE, returnTrain = FALSE)
    pred.i <- data.frame()
    for(fd in 1:k){
      df <- df.
      pred.k <- data.frame()
      geno.k <- df[flds[[fd]],"germplasmName"]
      df[df[,ID]%in%geno.k,tr] <- NA    #NA masking of the lines in the k-fold to be predicted
      df[,ID] <- as.factor(df[,ID])
      frm.cv <- paste(tr,"~ 1  ",collapse = " ")
      fix_frm.cv <- formula(paste(frm.cv, collapse = " "))
      mm.gb <-mmer(fixed = fix_frm.cv,
                   random =  ~  vsr(germplasmName ,Gu=G_mat.v) + vsr(germplasmName2 ,Gu=Gaa)  ,
                   rcov = ~vsr(units),
                   #tolParInv = 10000,
                   data = droplevels(df[df$germplasmName %in% rownames(G_mat.v),])
      )
      pred <- data.frame(germplasmName=names(mm.gb$U[[1]][[tr ]]),
                         gebv=as.numeric(mm.gb$U[[1]][[tr ]] + mm.gb$U[[2]][[tr ]]))
      
      rownames(pred) <- pred[,1]
      pred.k <- rbind(pred.k,pred[geno.k ,])
      pred <- merge(pred.k,df.[,c(ID,tr)],by=ID)
      acc.k <- c(acc.k,cor(pred[,2],pred[,3],use = "complete.obs"))
      pred.i <- rbind(pred.i,pred.k)
    }
    pred.i <- merge(pred.i,df.[,c(ID,tr)],by=ID)
    
    res1 <- rbind(res1,data.frame(
      trait=tr,n=nrow(df.),
      acc.k=mean(acc.k,na.rm = T),
      acc.i=cor(pred.i[,2],pred.i[,3],use = "complete.obs"))
    )
    
  }
}

# save results
write.csv(res1, "001 Phenotypic Data/CV results/MP6-8_2022_2_CV_AI model.csv") # Replace by the name of the trial

### Figure 3 ------------------------------------------------------------------#

# Clear workspace 
rm(list=ls())

# Load necessary libraries
library(tidyverse)

# MP2 CV Data
mp2_2020_A <- read.csv("MP2_2020_CV_A model.csv", header = TRUE)
mp2_2020_A$Material <- factor("MP2")
mp2_2020_A$Trial <- factor("MP2 2020")
mp2_2020_A$Model <- factor("Additive")

mp2_2021_1_A <- read.csv("MP2_2021_1_CV_A model.csv", header = TRUE)
mp2_2021_1_A$Material <- factor("MP2")
mp2_2021_1_A$Trial <- factor("MP2 2021_1")
mp2_2021_1_A$Model <- factor("Additive")

mp2_2021_2_A <- read.csv("MP2_2021_2_CV_A model.csv", header = TRUE)
mp2_2021_2_A$Material <- factor("MP2")
mp2_2021_2_A$Trial <- factor("MP2 2021_2")
mp2_2021_2_A$Model <- factor("Additive")

mp2_2020_AI <- read.csv("MP2_2020_CV_AI model.csv", header = TRUE)
mp2_2020_AI$Material <- factor("MP2")
mp2_2020_AI$Trial <- factor("MP2 2020")
mp2_2020_AI$Model <- factor("Additive + Epistasis")

mp2_2021_1_AI <- read.csv("MP2_2021_1_CV_AI model.csv", header = TRUE)
mp2_2021_1_AI$Material <- factor("MP2")
mp2_2021_1_AI$Trial <- factor("MP2 2021_1")
mp2_2021_1_AI$Model <- factor("Additive + Epistasis")

mp2_2021_2_AI <- read.csv("MP2_2021_2_CV_AI model.csv", header = TRUE)
mp2_2021_2_AI$Material <- factor("MP2")
mp2_2021_2_AI$Trial <- factor("MP2 2021_2")
mp2_2021_2_AI$Model <- factor("Additive + Epistasis")

# MP6-8 CV Data
mp68_2021_1_A <- read.csv("MP6-8_2021_1_CV_A model.csv", header = TRUE)
mp68_2021_1_A$Material <- factor("MP6-8")
mp68_2021_1_A$Trial <- factor("MP6-8 2021_1")
mp68_2021_1_A$Model <- factor("Additive")

mp68_2021_2_A <- read.csv("MP6-8_2021_2_CV_A model.csv", header = TRUE)
mp68_2021_2_A$Material <- factor("MP6-8")
mp68_2021_2_A$Trial <- factor("MP6-8 2021_2")
mp68_2021_2_A$Model <- factor("Additive")

mp68_2022_1_A <- read.csv("MP6-8_2022_1_CV_A model.csv", header = TRUE)
mp68_2022_1_A$Material <- factor("MP6-8")
mp68_2022_1_A$Trial <- factor("MP6-8 2022_1")
mp68_2022_1_A$Model <- factor("Additive")

mp68_2022_2_A <- read.csv("MP6-8_2022_2_CV_A model.csv", header = TRUE)
mp68_2022_2_A$Material <- factor("MP6-8")
mp68_2022_2_A$Trial <- factor("MP6-8 2022_2")
mp68_2022_2_A$Model <- factor("Additive")

mp68_2021_1_AI <- read.csv("MP6-8_2021_1_CV_AI model.csv", header = TRUE)
mp68_2021_1_AI$Material <- factor("MP6-8")
mp68_2021_1_AI$Trial <- factor("MP6-8 2021_1")
mp68_2021_1_AI$Model <- factor("Additive + Epistasis")

mp68_2021_2_AI <- read.csv("MP6-8_2021_2_CV_AI model.csv", header = TRUE)
mp68_2021_2_AI$Material <- factor("MP6-8")
mp68_2021_2_AI$Trial <- factor("MP6-8 2021_2")
mp68_2021_2_AI$Model <- factor("Additive + Epistasis")

mp68_2022_1_AI <- read.csv("MP6-8_2022_1_CV_AI model.csv", header = TRUE)
mp68_2022_1_AI$Material <- factor("MP6-8")
mp68_2022_1_AI$Trial <- factor("MP6-8 2022_1")
mp68_2022_1_AI$Model <- factor("Additive + Epistasis")

mp68_2022_2_AI <- read.csv("MP6-8_2022_2_CV_AI model.csv", header = TRUE)
mp68_2022_2_AI$Material <- factor("MP6-8")
mp68_2022_2_AI$Trial <- factor("MP6-8 2022_2")
mp68_2022_2_AI$Model <- factor("Additive + Epistasis")

# Combine all datasets
res <- rbind(
  mp2_2020_A, mp2_2021_1_A, mp2_2021_2_A, mp2_2020_AI, mp2_2021_1_AI, mp2_2021_2_AI,
  mp68_2021_1_A, mp68_2021_2_A, mp68_2022_1_A, mp68_2022_2_A,
  mp68_2021_1_AI, 
  mp68_2021_2_AI, mp68_2022_1_AI, mp68_2022_2_AI
)

res <- res %>% rename(Accuracy = acc.k)

# Calculate summary statistics
summary_data <- res %>%
  group_by(Model, Material, trait, Trial) %>%
  summarize(
    Mean_Accuracy = mean(Accuracy, na.rm = TRUE),
    SD_Accuracy = sd(Accuracy, na.rm = TRUE)
  )

write.csv(summary_data, "CV averages.csv")
write.csv(res, "CV all.csv")

res <- read.csv("CV all.csv")
# Create a mapping of traits to new names
trait_mapping <- c(
  "Yield.LSU_01.0000138" = "Grain Yield",
  "Days.to.heading.LSU_01.0000130" = "Days to Heading",
  "Whole.milling.percentage.LSU_01.0000100" = "Whole Milling %",
  "Chalk.impact.LSU_01.0000108" = "Chalk %"
)

# Replace trait names using the mapping
res$trait <- trait_mapping[res$trait]


# Generate boxplot (see also, figure 3 script)

ggplot(res, aes(x = Material, y = Accuracy, fill = Model)) +
  geom_boxplot(color = "black", position = position_dodge(preserve = "single")) +
  geom_jitter(color = "black", 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              alpha = 0.4, size = 0.5)+
  #scale_fill_manual(values = palette1) +
  facet_wrap(~trait) +
  theme_bw(16) +
  ggtitle("Cross Validation Prediction Accuracy\nof the MP Materials using the A or the A+I Models") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(y = "Prediction Accuracy", x = "Materials") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2))

ggsave("Figure 03 Epistasis Amat x Dmat with dots.jpg", width = 8, height = 10, dpi = 300)
