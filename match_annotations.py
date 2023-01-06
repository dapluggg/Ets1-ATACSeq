import pandas as pd 
import sys 

logFC_file = sys.argv[1]
bedtools_anno_file = sys.argv[2]
out_filename = sys.argv[3]

def match_annotations(logFC_file, bedtools_anno_file):
    '''
    Match the annotations from bedtols closest and apend them to the DESeq2 output. 
    ''' 
    # Read in BED file
    bed_colnames = ['chr', 'start', 'end', 'name', 'score', 'strand', 
                    'closest_chr', 'closest_gene_start', 'closest_gene_end', 'gene_name', 'closest_score', 'closest_strand',
                    'Gene_Database', 'type', 'unknown', 'match_string', 'distance']
    closest = pd.read_csv(bedtools_anno_file, sep='\t', header=None, names=bed_colnames)
    # Read in the DESeq2 output
    logFC_colnames = ["Region","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
    deseq2 = pd.read_csv(logFC_file, sep=',', header=None, names=logFC_colnames)
    # Merge the two dataframes
    merged = pd.merge(closest, deseq2, left_on='name', right_on='Region').filter(["Region","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj", "gene_name", "distance"])
    # Write the merged dataframe to a file
    merged.to_csv(f'{out_filename}', sep=',', index=False)

if __name__ == '__main__':
    match_annotations(logFC_file, bedtools_anno_file)
