#! /usr/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=15:00:00
#SBATCH --output=slurm_outputs/%j.out
module load bedtools;
bedtools slop -i /blue/m.pipkin/s.nagaraja/Tsuda_MegaExp_v2/merged/merged_peaks_sorted.bed -g /home/s.nagaraja/ATAC-seqPipeline/core/GCF_000001635.27_GRCm39_genomic.size -b 100000 > merged/genes_100kb_up_down.bed;
bedtools intersect -a merged/genes_100kb_up_down.bed -b /blue/m.pipkin/s.nagaraja/MusRef/GCF_000001635.27_GRCm39_genomic.gtf -wa -wb > merged/genes_100kb_up_down.bed.genes;
