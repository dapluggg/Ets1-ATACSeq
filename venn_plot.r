library(BioVenn)
set.seed(123)

args = commandArgs(trailingOnly = TRUE)
# Read narrowpeak into list 
narrowpeak_to_list <- function(narrowpeak){
    infile <- read.table(narrowpeak, header=FALSE, sep="\t")
    peakset <- paste(infile$V1, infile$V2, infile$V3, sep = "_")

    return(peakset)
}

draw_venn <- function(peakset_A, samplename_A, peakset_B, samplename_B){
    list_A <- narrowpeak_to_list(peakset_A)
    list_B <- narrowpeak_to_list(peakset_B)

    filename_constructor <- paste0("/blue/m.pipkin/s.nagaraja/Tsuda_Ets1_MegaExp/venn/", samplename_A, "_", samplename_B, ".svg")
    print(filename_constructor)
    svg(filename=filename_constructor)
    draw.venn(list_A, list_B, NULL, 
              xtitle = samplename_A, ytitle = samplename_B,
              title = paste0(samplename_A, "\n", " and ", "\n", samplename_B),
              t_s = 0.5, xt_s = 0.5, yt_s = 0.5,
              subtitle = "",
              t_f = "arial", xt_f = "arial", yt_f = "arial", nr_f = "arial"
    )
    dev.off()
    return()
}

main <- function(){
    peakset_A <- args[1]
    samplename_A <- args[2]
    peakset_B <- args[3]
    samplename_B <- args[4]
    draw_venn(peakset_A, samplename_A, peakset_B, samplename_B)
}

main()