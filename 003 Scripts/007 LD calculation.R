###############################################################################
# Script: LD decay caluculatrion and Haplotype sizes
# Author: Maria Guadalupe Montiel, Ph.D.
# Email: mmontiel@horizonseed.com
# Paper: The Effect of Haplotype Size on Genomic Selection Accuracy and Epistasis:
#        An Empirical Study in Rice
# Description: This script analyzes LD decay across different rice populations (MP2, MP4, MP6-8). 
# It calculates LD decay based on SNP distances and visualizes the results using ggplot2.
# 
###############################################################################
# Clear workspace and load libraries
rm(list=ls())
library(sommer)
library(ggplot2)
library(gridExtra)
library(snpReady)
library(ggplot2)
library(paletteer)


#### Reading files and cleaning the data ####
# MP2
mp2_2020 <- read.csv("001 Phenotypic Data/MP materials/project_20mp2.csv", header = TRUE)
mp2_2021_1 <- read.csv("001 Phenotypic Data/MP materials/project_21mp2_1.csv", header = TRUE)
mp2_2021_2 <- read.csv("001 Phenotypic Data/MP materials/project_21mp2_2.csv", header = TRUE)
mp2 <- rbind(mp2_2020, mp2_2021_1, mp2_2021_2)

# MP4
mp4<- read.csv("001 Phenotypic Data/MP materials/project_mp4_pheno.csv", header = TRUE)

# MP6-8
mp68_2021_1 <- read.csv("001 Phenotypic Data/MP materials/project_21mp68_1.csv", header = TRUE)
mp68_2021_2 <- read.csv("001 Phenotypic Data/MP materials/project_21mp68_2.csv", header = TRUE)
mp68_2022_1 <- read.csv("001 Phenotypic Data/MP materials/project_22mp68_1.csv", header = TRUE)
mp68_2022_2 <- read.csv("001 Phenotypic Data/MP materials/project_22mp68_2.csv", header = TRUE)
mp68 <- rbind(mp68_2021_1, mp68_2021_2, mp68_2022_1, mp68_2022_2)

#### Prepare and Load genotype data ####
# Perform quality control for each genotype dataset

# MP2 QC
geno_mp2 <- as.matrix(read.csv("002 Genotypic Data/MP2_Numeric_wright.csv", row.names = 1))
geno_mp2 <- geno_mp2[rownames(geno_mp2) %in% mp2$germplasmName, ]
geno_mp2[1:5, 1:5]

# MP4 QC
geno_mp4 <- as.matrix(read.csv("002 Genotypic Data/MP4_Numeric_wright.csv  ", row.names = 1))
geno_mp4 <- geno_mp4[rownames(geno_mp4) %in% mp4$germplasmName, ]
dim(geno_mp4)
geno_mp4[1:5, 1:5]

# MP6-8 QC
geno_mp68 <- as.matrix(read.csv("002 Genotypic Data/MP6-8_Numeric_wright.csv", row.names = 1))
geno_mp68 <- geno_mp68[rownames(geno_mp68) %in% mp68$germplasmName, ]
dim(geno_mp68)
geno_mp68[1:5, 1:5]

#### Get the marker matrix ####
CPgeno_mp2 <- geno_mp2  # MP2 genotypic data
#CPgeno_mp2 <- CPgeno_mp2 - 1

CPgeno_mp68 <- geno_mp68  # MP68 genotypic data
#CPgeno_mp68 <- CPgeno_mp68 - 1

CPgeno_mp4 <- geno_mp4  # MP68 genotypic data
#CPgeno_mp4 <- CPgeno_mp4 - 1

#### Get the map ####
#### Reading the SNP distance in cM ####
distance <- read.csv("002 Genotypic Data/LDdistance.csv", header = TRUE)
str(distance)
distance$Position <- as.numeric(distance$Position)

mapCP <- distance; head(mapCP)
names(mapCP) <- c("Locus","Position","LG")
mapCP <- mapCP[which(mapCP$LG <= 12),]

#### LD decay , r value calculation for the MP2, MP4, and MP6-8 ####
res_mp2 <- LD.decay(CPgeno_mp2, mapCP)
res_mp4 <- LD.decay(CPgeno_mp4, mapCP)
res_mp68 <- LD.decay(CPgeno_mp68, mapCP)

# Calculate mean and median r2 before filtering
mean(res_mp2$all.LG$r2)
mean(res_mp4$all.LG$r2)
mean(res_mp68$all.LG$r2)

# Subset only markers with significant LD for all pops
res_mp2$all.LG <- res_mp2$all.LG[which(res_mp2$all.LG$p < .001),]
res_mp4$all.LG <- res_mp4$all.LG[which(res_mp4$all.LG$p < .001),]
res_mp68$all.LG <- res_mp68$all.LG[which(res_mp68$all.LG$p < .001),]

######## Metrics MP2
#Haplotypes
mean(res_mp2$all.LG$d)
median(res_mp2$all.LG$d)
quantile(res_mp2$all.LG$d, probs = c(0.1, 0.25, 0.75, 0.9))

#LD decay in r2
plot(density(res_mp2$all.LG$r2))
plot(density(res_mp2$all.LG$d))
hist(res_mp2$all.LG$r2)

mean(res_mp2$all.LG$r2)
median(res_mp2$all.LG$r2)

###### Metrics MP4
#Haplotypes
mean(res_mp4$all.LG$d)
median(res_mp4$all.LG$d)
quantile(res_mp4$all.LG$d, probs = c(0.1, 0.25, 0.75, 0.9))

#LD decay and r2
mean(res_mp4$all.LG$r2)
median(res_mp4$all.LG$r2)

###### Metrics MP6-8
#Haplotypes
mean(res_mp68$all.LG$d)
median(res_mp68$all.LG$d)
quantile(res_mp68$all.LG$d, probs = c(0.1, 0.25, 0.75, 0.9))

#LD decay and r2
mean(res_mp68$all.LG$r2)
median(res_mp68$all.LG$r2)

#### Plot LD decay for each population ####
res_mp2 <- LD.decay(CPgeno_mp2, mapCP)
res_mp4 <- LD.decay(CPgeno_mp4, mapCP)
res_mp68 <- LD.decay(CPgeno_mp68, mapCP)

plot_ld_decay <- function(data, title) {
  avg_r2 <- mean(data$r2)
  ggplot(data, aes(x = d, y = r2)) +
    geom_point(color = "cadetblue4", alpha = 0.5) +
    #geom_smooth(method = "loess", se = FALSE, color = "blue") +
    geom_hline(yintercept = avg_r2, linetype = "dashed", color = "red") +
    labs(title = title, x = "Physical distance (Kb)", y = "r2") +
    theme(text = element_text(size = 25)) +
    theme_light()
}


p1 <- plot_ld_decay(res_mp2$all.LG, "a) MP2 - Recombination cycles = 2")
p2 <- plot_ld_decay(res_mp4$all.LG, "b) MP4 - Recombination cycles = 5")
p3 <- plot_ld_decay(res_mp68$all.LG,"c) MP6-8 - Recombination cycles = 6")

grid_plot <- grid.arrange(p1, p2, p3, ncol = 1)
ggsave("LD_decay_plots.png", grid_plot, width = 10, height = 15)

