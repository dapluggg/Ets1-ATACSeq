# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("sinhrks/ggfortify")
library(ggplot2)
library(plotly)
library(stringr)
library(readr)
library(plotly)
library(ggfortify)
print("START")
# Make 4 plots from the TPM and log2FC data: 
# 1. PCA of TPMs 
# 2. Heatmap of TPMs
# 3. PCA of log2FCs
# 4. Heatmap of log2FCs

# Load in TPM data 
TPM.matrix <- read.table("counts/tpm_counts.txt", header = T, sep = "\t")
#str(TPM.matrix)

# Load in metadata 
meta <- read.table("Tsuda_MegaExp.txt", header = T, sep = "\t")

# # Remove long sample names 
# long.names <- colnames(TPM.matrix)
# short.names <- str_replace(long.names, "X.blue.m.pipkin.s.nagaraja.", "")
# short.names <- str_replace(short.names, "X.blue.m.pipkin.stsuda.", "")
# short.names <- str_replace(short.names, ".bams_noDups", "")
# colnames(TPM.matrix) <- short.names

# Transpose dataframe for PCA 
# TPM.matrix <- t(TPM.matrix)

# Ets1 KO Experiments
samples <- c(meta[(meta$Exp == "Ets1KO_invitro"), c("TPM_column")])
Ets1.KO <- TPM.matrix[,samples]
str(Est1.KO)


# TPM.PCA <- prcomp(t(TPM.matrix[,7:ncol(TPM.matrix)]), center = T)
# summary(TPM.PCA)
# # Save PCA Plot to PDF
# # plot PCA of TPMs
# pdf("counts/TPM_PCA.pdf")
# autoplot(TPM.PCA, data = t(TPM.matrix[,7:ncol(TPM.matrix)]), label = TRUE)
# dev.off()