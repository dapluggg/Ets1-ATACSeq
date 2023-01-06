#!/usr/bin/env Rscript

# Deseq2 optimized for paralell excecution of contrast
library(DESeq2, quietly=TRUE)
library(stringr)
# library(ashr)
# library(EnhancedVolcano)
# library(scales)



###################////////////////////////
###################////////////////////////
# # First Run this block of code to create the dds object
# # This will require a pretty beefy node: >128 cores, ~1TB RAM
# # ^ No that is not a typo. I mean 1 terabyte of RAM.
# cts = read.table("counts/counts.mtx", header=TRUE, row.names=1)
# ssheet = read.table("Tsuda_MegaExp.txt", header=TRUE, stringsAsFactors=TRUE)
# results_dir = "deseq2"

# ####
# coldata <- ssheet[, c("SampleName", "Status")]
# rownames(coldata) <- ssheet[, "SampleName"]
# # print(coldata)

# # Fix featurecounts output names; replace paths with sample names
# name.vec <- c("Chr", "Start", "End", "Strand", "Length")
# for (i in 6:ncol(cts)) {
#     name.vec <- append(name.vec, as.character(ssheet$SampleName[i-5]))
# }
# names(cts) <- name.vec


# mat = cts[,6:ncol(cts)]
# mat[is.na(mat)] <- 0

# # print(str(mat))
# # print(head(mat))

# # print(str(coldata))
# # print(head(coldata))

# dds <- DESeqDataSetFromMatrix(countData = mat,
#                               colData = coldata,
#                               design= ~ Status)
# dds <- DESeq(dds, parallel=TRUE)
# saveRDS(dds, "dds.rds")

<<<<<<< HEAD
# print("RDS has been writtwn to dds.rds")
# # resultsNames(dds)
# ###################////////////////////////
# ###################////////////////////////




#################**********************
###################**********************
# Then comment the above block out and run the below block 
print("Parsing Args")
args = commandArgs(trailingOnly=TRUE)
print("Parsed")
cts = read.table(args[1], header=TRUE)
print("Counts")
ssheet = read.table(args[2], header=TRUE, stringsAsFactors=TRUE)
print("SSheet")
results_dir = args[3]
print("results dir")
control = args[4]
print("control")
treatment = args[5]
print("treatment")
# Assumes the results_dir has already been created by python 
# Assumes that the order of samples is the exact same in cts and ssheet
# To enable parallelization, this script only does one contrast per job.
# You must submit this script separately for each contrast.
dds <- readRDS("dds.rds")

print(paste("Contrast: ", treatment, " vs ", control))
res <- results(dds, contrast=c("Status", treatment, control),
               independentFiltering=TRUE, parallel=TRUE)
print(paste0(results_dir, "/", treatment, "_vs_", control, ".csv"))
write.csv(as.data.frame(res), paste0(results_dir, "/", treatment, "_vs_", control, ".csv"))

##################**********************
##################**********************
=======
dds <- readRDS("/blue/m.pipkin/s.nagaraja/Tsuda_Ets1_MegaExp/dds.rds")
resultsNames(dds)
#########
# To enable parallelization, this script only does one contrast per job.
# You must submit this script separately for each contrast.
print(paste("Contrast: ", control, " vs ", treatment))
res <- results(dds, contrast=c("Status", control, treatment),
               independentFiltering=TRUE, parallel=TRUE)
print(paste0(results_dir, "/", control, "_vs_", treatment, ".csv"))
write.csv(as.data.frame(res), paste0(results_dir, "/", control, "_vs_", treatment, ".csv"))



#########
# for (control in unique(ssheet$Status)) {
#     for (treatment in unique(ssheet$Status)) {
#         if (control == treatment){
#             next
#         } else {
#             print(paste("Contrast: ", control, " vs ", treatment))
#             res <- results(dds, contrast=c("Status", control, treatment),
#                         independentFiltering=TRUE, parallel=TRUE)
#             # res <- lfcShrink(dds, contrast=c("Status", control, treatment), res=res, type="ashr")
#             print(paste0(results_dir, "/", control, "_vs_", treatment, ".csv"))
#             write.csv(as.data.frame(res), paste0(results_dir, "/", control, "_vs_", treatment, ".csv"))

#             # print(paste0(results_dir, "/", control, "_vs_", treatment, ".pdf"))
#             # print(paste(control, "vs.", treatment))
#             # volcano <- EnhancedVolcano(res,
#             #                 lab = rownames(res),
#             #                 x = 'log2FoldChange',
#             #                 y = 'pvalue',
#             #                 title = paste(control, "vs.", treatment),
#             #                 pCutoff = 0.05,
#             #                 FCcutoff = 1.2
#             #                 )
#             # ggsave(paste0(results_dir, "/", control, "_vs_", treatment, ".pdf"), volcano)

#             # # # Make pairwse peak/count plots for each sample
#             # control_samples <- as.character(ssheet[ssheet$Status == control,"SampleName"])
#             # treatment_samples <- as.character(ssheet[ssheet$Status == treatment,"SampleName"])
#             # # rowMeans behaes extremely stupidly and doesn't work on a dataframe with only one column 
#             # # This is why R sucks...
#             # if (length(control_samples) > 1) {
#             #     control_counts <- rowMeans(mat[, control_samples])
#             # } else {
#             #     control_counts <- mat[, control_samples]
#             # }
#             # if (length(treatment_samples) > 1) {
#             #     treatment_counts <- rowMeans(mat[, treatment_samples])
#             # } else {
#             #     treatment_counts <- mat[, treatment_samples]
#             # }
#             # ppc_df <- as.data.frame(cbind(control_counts, treatment_counts))
#             # ppc <- ggplot(ppc_df, aes(x=control_counts,
#             #                         y=treatment_counts),
#             #             ) +
#             #     geom_density_2d() +
#             #     geom_abline(intercept=0, slope=1, color="Black") +
#             #     scale_y_continuous(trans='log10', 
#             #                         breaks = trans_breaks('log10', function(x) 10^x),
#             #                         labels = trans_format('log10', math_format(10^.x))
#             #                         ) +
#             #     scale_x_continuous(trans='log10', 
#             #                         breaks = trans_breaks('log10', function(x) 10^x),
#             #                         labels = trans_format('log10', math_format(10^.x))
#             #                         ) +
#             #     theme_minimal() +
#             #     xlab(paste0(control)) +
#             #     ylab(paste0(treatment)) +
#             #     ggtitle("ATAC-seq mean counts/peak")
#             # ###### These plots may also be used
#             # # ppc <- ggplot(ppc_df, aes(x=control_counts,
#             # #                           y=treatment_counts),
#             # #               ) +
#             # #        geom_hex() +
#             # #        geom_abline(intercept=0, slope=1, color="Black") +
#             # #        scale_fill_continuous(type = "viridis") +
#             # #        scale_y_continuous(trans='log10') +
#             # #        scale_x_continuous(trans='log10')
#             # # ppc <- ggplot(ppc_df, aes(x=control_counts,
#             # #                           y=treatment_counts),
#             # #               ) +
#             # #        geom_point(alpha=0.05) +
#             # #        geom_abline(intercept=0, slope=1, color="Black") +
#             # #        scale_y_continuous(trans='log10') +
#             # #        scale_x_continuous(trans='log10')
#             # #######
#             # ggsave(paste0(results_dir, "/", control, "_vs_", treatment, "_ppc",".pdf"), ppc)
#         }
#     }
# }

#########
main <- function(){
}
>>>>>>> 3b2999f8eb69916458f4b9fa662c322abde97b2d
