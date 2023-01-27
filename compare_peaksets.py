import pandas as pd 
import os, sys
import itertools
from tqdm import tqdm
import time
import glob
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import pyranges as pr
print("IMPORTED")

# Compare peaksets in a pairwise comparison
# First, merge all peaks using HOMER 
# Then plot venn diagrams 
# Maybe draw a heatmap? 

def merge(peaksets_meta):
    '''
    Concatenate all peaksets into one file and merge. 
    '''
    if not os.path.exists('venn/'):
        os.makedirs('venn/')
    if not os.path.exists('merged/'):
        os.makedirs('merged/')

    peaksets_df = pd.read_csv(peaksets_meta, sep='\t')
    all_peaksets = (peaksets_df['path'] + '/macs2/' + peaksets_df['SampleName_inExp'] + '_peaks.narrowPeak').tolist()

    # concatenate peaks into one file
    with open('merged/concat_peaks.narrowPeak', 'w') as outfile: 
        for narrowpeak in tqdm(all_peaksets): 
            with open(narrowpeak) as infile:
                for line in infile: 
                    outfile.write(line)

    cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
module load bedtools;
bedtools sort -i merged/concat_peaks.narrowPeak > merged/merged_peaks_sorted.bed;
bedtools merge -i merged/merged_peaks_sorted.bed > merged/merged_peaks.bed;
    '''
    f = open(f'merged/merge_peaks.sh', 'w+')
    f.write(cmd)
    f.close()
    os.system(f'sbatch merged/merge_peaks.sh')

    return None 

def replace_peaks(peaksets_meta):
    '''
    Replace peaks in all peaksets with merged peaks unsing bedtools intersect. 
    '''
    # create directory for replacing peaks
    if not os.path.exists('replace_peaks/'):
        os.makedirs('replace_peaks/')

    peaksets_df = pd.read_csv(peaksets_meta, sep='\t')
    all_peaksets = (peaksets_df['path'] + '/macs2/' + peaksets_df['SampleName_inExp'] + '_peaks.narrowPeak').tolist()
    
    for index, sample in peaksets_df.iterrows():
        sample_name = sample['SampleName']
        narrowpeak = sample['path'] + '/macs2/' + sample['SampleName_inExp'] + '_peaks.narrowPeak'

        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
module load bedtools;
bedtools intersect -wb -a {narrowpeak} -b merged/merged_peaks.bed > replace_peaks/{sample_name}.narrowPeak;
cut -f 11,12,13 replace_peaks/{sample_name}.narrowPeak > replace_peaks/{sample_name}.bed;
        '''

        f = open(f'replace_peaks/{sample_name}.sh', 'w+')
        f.write(cmd)
        f.close()
        os.system(f'sbatch replace_peaks/{sample_name}.sh')

    return None

def draw_venn(peaksets_meta):
    '''
    Draw venn diagrams for all peaksets. 
    '''
    # create directory for venn diagrams
    if not os.path.exists('venn/'):
        os.makedirs('venn/')

    peaksets_df = pd.read_csv(peaksets_meta, sep='\t')
    all_samples = peaksets_df['SampleName'].tolist()
    
    combs = list(itertools.combinations(all_samples, 2))
    print(f'Submitting {len(combs)} jobs.')
    
    for comb in combs:
        sample_A = f'replace_peaks/{comb[0]}.bed'
        sample_B = f'replace_peaks/{comb[1]}.bed'

        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=2:00:00
conda activate ATACseq_env;
python3 /blue/m.pipkin/s.nagaraja/Tsuda_Ets1_MegaExp/venn_diagrams_modified.py \
-A {sample_A} -aname {comb[0]} \
-B {sample_B} -bname {comb[1]} \
-O /blue/m.pipkin/s.nagaraja/Tsuda_Ets1_MegaExp/venn;
    '''
        f = open(f'venn/{comb[0]}_{comb[1]}.sh', 'w+')
        f.write(cmd)
        time.sleep(0.1)
        f.close()
        time.sleep(0.1)
        os.system(f'sbatch venn/{comb[0]}_{comb[1]}.sh')
        time.sleep(0.1)

    return None

def bed2gtf(infile, outfile):
    '''
    Convert merged BED peaks  to GTF for use in feature Counts. Stealing this function from ATACseqPipeline.  
    '''
    import bed2gtf
    bed2gtf.run(infile, outfile)

    return None 

def count_reads(merged_peaks, peaksets_meta, counts):
    '''
    Use FeatureCounts to count reads in merged bed (now converted to GTF)
    merged_peaks: GTF file of merged peaks. 
    peaksets_meta: TXT file of samples. 
    counts: Directory in which output count matrices should be written. 
    '''
    if not os.path.exists('counts/'):
        os.makedirs('counts/')

    peaksets_df = pd.read_csv(peaksets_meta, sep='\t')

    endedness_dict = {
        'PE' : peaksets_df[peaksets_df['Endedness']=='PE'],
        'SE' : peaksets_df[peaksets_df['Endedness']=='SE']
    }
    endedness_cmd_option = ' -p '
    for ended, peakset in endedness_dict.items():
        bam_files = peakset['path'] + '/bams_noDups/' + peakset['SampleName_inExp'] + '.bam'
        bam_files = ' '.join(bam_files)
        if ended == 'PE':
            pass
        elif ended == 'SE':
            endedness_cmd_option = ' '
        else: 
            print('Unknown ended type.')
            return False
            
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output=slurm_outputs/slurm-%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load subread;
featureCounts -T 16 -a {merged_peaks} -t 'peak' -g 'peak_id'{endedness_cmd_option}-o {counts}/counts_{ended}.mtx {bam_files};
        '''
        f = open(f'counts/peak_counts_{ended}.sh', 'w+')
        f.write(cmd)
        f.close()
        os.system(f'sbatch counts/peak_counts_{ended}.sh')

    return None 

def differential_peaks(peaksets_meta):
    '''
    Run DESeq2 to test for differential peaks. Use marked samples as controls. 
    '''
    if not os.path.exists('deseq2/'):
        os.makedirs('deseq2/')

    peaksets_df = pd.read_csv(peaksets_meta, sep='\t')
    
    # Iterate through all permutations of samples. 
    perms = list(itertools.permutations(peaksets_df['Status'].tolist(), 2))
    print(f'Submitting {len(perms)} jobs.')

    for perm in tqdm(perms): 
        control, treatment = perm[0], perm[1]

        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=4:00:00
#SBATCH --output=slurm_outputs/slurm-%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load R;
Rscript differential_exp.r {'counts/counts.mtx'} {peaksets_meta} {'deseq2/'} {control} {treatment};
        '''
        f = open(f'deseq2/differential_analysis_{control}_vs_{treatment}.sh', 'w+')
        f.write(cmd)
        f.close()
        # os.system(f'sbatch deseq2/differential_analysis.sh')


    return None

def filter_up_down(peaksets_dir, fc_cutoff, pval_cutoff):
    '''
    Assign peaks as up or downregulated based on fold change and p-value. Then draw pairwise venn diagrams. 
    '''
    if not os.path.exists('venn_up_down/'):
        os.makedirs('venn_up_down/')
    if not os.path.exists('venn_up_down/filtered/'):
        os.makedirs('venn_up_down/filtered/')
    
    for filepath in tqdm(glob.iglob(f'{peaksets_dir}/*.csv')):
        file_string = filepath.split('/')[-1]
        file_string = file_string.replace('.csv', '')
        # control = file_string.split('_vs_')[0]
        # treatment = file_string.split('_vs_')[1]
        # print(f'{control} ------ {treatment}')

        df = pd.read_csv(filepath)
        
        # Filter by fold change and p-value.
        up_peaks = df[(df['log2FoldChange'] > fc_cutoff[0]) & (df['pvalue'] < pval_cutoff)]
        down_peaks = df[(df['log2FoldChange'] < fc_cutoff[1]) & (df['pvalue'] < pval_cutoff)]

        up_peaks.to_csv(f'venn_up_down/filtered/{file_string}_up.csv', index=False)
        down_peaks.to_csv(f'venn_up_down/filtered/{file_string}_down.csv', index=False)
        
    return None 

def venn_up_down(peaksets_dir): 
    '''
    Draw Venn Diagram for up and downregulated peaks. DO this for every pairwise comparison.
    This is not the cleanest code here, but it works. 
    '''
    filtered_up = sorted(glob.glob(f'{peaksets_dir}/filtered/*_up.csv'))
    filtered_down = sorted(glob.glob(f'{peaksets_dir}/filtered/*_down.csv'))

    def get_peakslist_from_csv(csv_file):
        df = pd.read_csv(csv_file)
        df.rename( columns={'Unnamed: 0':'peak'}, inplace=True )
        peakslist = df['peak'].tolist()
        return set(peakslist)

    def rename_sample(csv_file):
        with_control = csv_file.split('/')[-1]
        with_control = with_control.replace('_up.csv', '')
        with_control = with_control.replace('_down.csv', '')
        without_control = with_control.split('_vs_')[1]
        return without_control

    combs_up = list(itertools.combinations(filtered_up, 2))
    combs_down = list(itertools.combinations(filtered_down, 2))

    for comb in tqdm(combs_up): 
        sample_A_name, sample_B_name = rename_sample(comb[0]), rename_sample(comb[1])
        peaks_list_A, peaks_list_B = get_peakslist_from_csv(comb[0]), get_peakslist_from_csv(comb[1])
        
        # Compute intersection and difference between A and B
        A_only = set(peaks_list_A) - set(peaks_list_B)
        intersection = set(peaks_list_A) & set(peaks_list_B)
        B_only = set(peaks_list_B) - set(peaks_list_A)
        total = len(A_only) + len(B_only) + len(intersection)

        venn2(subsets=(len(A_only), len(B_only), len(intersection)), set_labels=(sample_A_name, sample_B_name), alpha=0.25)
        venn2_circles(subsets=(len(A_only), len(B_only), len(intersection)), linewidth=0.25)
        plt.title('Regions of Increased Accessibility')
        plt.savefig(f'{peaksets_dir}/{sample_A_name}_and_{sample_B_name}_up.svg', bbox_inches = "tight")
        plt.close()

    for comb in tqdm(combs_down): 
        sample_A_name, sample_B_name = rename_sample(comb[0]), rename_sample(comb[1])
        peaks_list_A, peaks_list_B = get_peakslist_from_csv(comb[0]), get_peakslist_from_csv(comb[1])
        
        # Compute intersection and difference between A and B
        A_only = set(peaks_list_A) - set(peaks_list_B)
        intersection = set(peaks_list_A) & set(peaks_list_B)
        B_only = set(peaks_list_B) - set(peaks_list_A)
        total = len(A_only) + len(B_only) + len(intersection)

        venn2(subsets=(len(A_only), len(B_only), len(intersection)), set_labels=(sample_A_name, sample_B_name), alpha=0.25)
        venn2_circles(subsets=(len(A_only), len(B_only), len(intersection)), linewidth=0.25)
        plt.title('Regions of Decreased Accessibility')
        plt.savefig(f'{peaksets_dir}/{sample_A_name}_and_{sample_B_name}_down.svg', bbox_inches = "tight")
        plt.close()

    return None

def venn_up_down_3way():
    '''
    Draw three-way venn diagram for up and down regulated peaks. Only do this for a set of supplied conditions from a txt file. 
    '''
    
    return None

def find_motifs(logFC_dir, motif_dir, motif_length, genome):
    '''
    Find motifs in up and downregulated peaks using HOMER. DESeq2 output is assumed to be pre-filtered on p-value (0.05) and fold change (< -1 or > 1).
    Assumes mm10 genome.
    '''
    if not os.path.exists('HOMER/'):
        os.makedirs('HOMER/')
    if not os.path.exists('HOMER/bed/'):
        os.makedirs('HOMER/bed/')
    if not os.path.exists('HOMER/submissions/'):
        os.makedirs('HOMER/submissions/')
    if not os.path.exists('HOMER/motifs/'):
        os.makedirs('HOMER/motifs/')
    if not os.path.exists('HOMER/preparsed/'):
        os.makedirs('HOMER/preparsed/')
    
    def recontstruct_bed(logFC_file):
        '''
        Reconstruction of bed file from DESeq2 output.
        '''
        sample_name = logFC_file.split('/')[-1]
        sample_name = sample_name.replace('.csv', '')
        df = pd.read_csv(logFC_file)
        # print(df.columns)
        # Rename the first column to 'Region' irrespective of what it is called in the input file
        df.rename( columns={df.columns[0]:'Region'}, inplace=True )
        # df.rename( columns={'Unnamed: 0':'peak'}, inplace=True )
        # df.rename( columns={'x':'peak'}, inplace=True )
        df['chr'] = df['Region'].str.split('_').str[0] + '_' + df['Region'].str.split('_').str[1]
        df['start'] = df['Region'].str.split('_').str[2]
        df['end'] = df['Region'].str.split('_').str[3]
        df['strand'] = '+'
        df['name'] = df['Region']
        df = df[['name', 'chr', 'start', 'end', 'strand']]
        df.to_csv(f'HOMER/bed/{sample_name}.bed', sep='\t', index=False, header=False)
        return None

    for filepath in tqdm(glob.iglob(f'{logFC_dir}/*.csv')):
        
        file_string = filepath.split('/')[-1]
        file_string = file_string.replace('.csv', '')
        # Create all the folders 
        if not os.path.exists(f'{motif_dir}/{file_string}'):
            os.makedirs(f'{motif_dir}/{file_string}')
        if not os.path.exists(f'HOMER/preparsed/{file_string}'):
            os.makedirs(f'HOMER/preparsed/{file_string}')
        #
        recontstruct_bed(filepath)
        bed_path = f'HOMER/bed/{file_string}.bed'
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=4:00:00
#SBATCH --output=slurm_outputs/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load homer;
findMotifsGenome.pl {bed_path} {genome} {motif_dir}/{file_string} -size {motif_length} -preparsedDir HOMER/preparsed/{file_string} -p 16 -mis 3;
        '''
        f = open(f'HOMER/submissions/{file_string}.sh', 'w+')
        f.write(cmd)
        f.close()
        os.system(f'sbatch HOMER/submissions/{file_string}.sh')
    return None

def counts_to_tpm(counts):
    '''
    Convert raw counts (from featureCounts) to TPM.
    '''
    path = os.path.dirname(counts)
    counts = pd.read_table(counts, sep='\t', header=0, index_col=False, comment='#')
    matrix = counts.iloc[:, 6:]
    lengths = counts.iloc[:, 5]
    print(matrix.shape)
    matrix = matrix.div(lengths, axis=0)
    scaling_factors = matrix.sum(axis=0) / 1e6
    print(scaling_factors.shape)
    matrix = matrix.div(scaling_factors, axis=1)
    matrix = matrix.round(3)
    print(matrix.shape)
    tpm_counts = pd.concat([counts.iloc[:, 0:6], matrix], axis=1)
    tpm_counts.to_csv(f'{path}/tpm_counts.txt', sep='\t', index=False)

    return None

def annotate_logFC(logFC_dir, genes_bed):
    '''
    Annotate the logFC files with gene names and gene IDs. This is some 
    gnarly code but it works, so whatever. 
    Here is the idea: 

    Take BED file -> convert to HOMER format -> annotate with bedtools closest ->
    cut -f 2,3,4,5,6 merged_peaks_HOMER.bed | tail -n +2 | sort -k1,1V -k2,2n -k3,3n > merged_peaks_HOMER.formatted.bed;
    bedtools closest -a merged_peaks_HOMER.formatted.bed -b /blue/m.pipkin/s.nagaraja/MusRef/Mouse39Genes.bed -D -t first > matched_genes.txt;
    '''
    if not os.path.exists('annotated_logFC/'):
        os.makedirs('annotated_logFC/')
    if not os.path.exists('annotated_logFC/bed/'):
        os.makedirs('annotated_logFC/bed/')
    if not os.path.exists('annotated_logFC/sorted_bed/'):
        os.makedirs('annotated_logFC/sorted_bed/')
    if not os.path.exists('annotated_logFC/annotated_bed/'):
        os.makedirs('annotated_logFC/annotated_bed/')
    if not os.path.exists('annotated_logFC/submissions/'):
        os.makedirs('annotated_logFC/submissions/')
    
    for filepath in tqdm(glob.iglob(f'{logFC_dir}/*.csv')):
        sample_name = filepath.split('/')[-1]
        sample_name = sample_name.replace('.csv', '')
        
        bed_file = f'annotated_logFC/bed/{sample_name}.bed'
        sorted_bed = f'annotated_logFC/sorted_bed/{sample_name}.sorted.bed'
        annotated_bed = f'annotated_logFC/annotated_bed/{sample_name}.annotated.bed'
        sub_file = f'annotated_logFC/submissions/annotate_{sample_name}.sh'

        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=4:00:00
#SBATCH --output=slurm_outputs/%j.out
source ~/.bashrc;
conda activate ATACseq_env;
module load bedtools;
python3 reconstruct_bed.py {filepath};
bedtools sort -i {bed_file} > {sorted_bed};
bedtools closest -a {sorted_bed} -b {genes_bed} -D b -t first > {annotated_bed};
    '''
        f = open(f'{sub_file}', 'w+')
        f.write(cmd)
        f.close()
        os.system(f'sbatch {sub_file}')
    return None 

def match_anno_to_logFC(logFC_dir, anno_dir):
    '''
    Match the annotated BED files to the logFC files. 
    '''
    if not os.path.exists('annotated_logFC/'):
        os.makedirs('annotated_logFC/')

    for filepath in tqdm(glob.iglob(f'{logFC_dir}/*.csv')):
        sample_name = filepath.split('/')[-1]
        sample_name = sample_name.replace('.csv', '')

        anno_file = f'{anno_dir}/{sample_name}.annotated.bed'
        sub_file = f'annotated_logFC/submissions/match_{sample_name}.sh'
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=4:00:00
#SBATCH --output=slurm_outputs/%j.out
source ~/.bashrc;
conda activate ATACseq_env;
python3 match_annotations.py {filepath} {anno_file} annotated_logFC/{sample_name}.annotated.csv;
        '''
        f = open(f'{sub_file}', 'w+')
        f.write(cmd)
        f.close()
        os.system(f'sbatch {sub_file}')
        
    return None

def intersection_motifs(logFC_dir, motif_dir, motif_length, genome):
    '''
    Find motifs from an intersection of regions. This function is almost the same as find_motifs, except the
    column naming convention and the feature padding. Padding is assumed to be in bp. 
    '''
    if not os.path.exists('venn_intersections/HOMER/'):
        os.makedirs('venn_intersections/HOMER/')
    if not os.path.exists('venn_intersections/HOMER/bed/'):
        os.makedirs('venn_intersections/HOMER/bed/')
    if not os.path.exists('venn_intersections/HOMER/submissions/'):
        os.makedirs('venn_intersections/HOMER/submissions/')
    if not os.path.exists('venn_intersections/HOMER/motifs/'):
        os.makedirs('venn_intersections/HOMER/motifs/')
    if not os.path.exists('venn_intersections/HOMER/preparsed/'):
        os.makedirs('venn_intersections/HOMER/preparsed/')
    if not os.path.exists('venn_intersections/HOMER/slopped/'):
        os.makedirs('venn_intersections/HOMER/slopped/')
    
    def recontstruct_bed(logFC_file):
        '''
        Reconstruction of bed file from DESeq2 output.
        '''
        sample_name = logFC_file.split('/')[-1]
        sample_name = sample_name.replace('.csv', '')
        df = pd.read_csv(logFC_file)
        df.rename( columns={'x':'peak'}, inplace=True )
        df['chr'] = df['peak'].str.split('_').str[0] + '_' + df['peak'].str.split('_').str[1]
        df['start'] = df['peak'].str.split('_').str[2]
        df['end'] = df['peak'].str.split('_').str[3]
        df['strand'] = df['peak'].str.split('_').str[4]
        df['score'] = 0
        df['name'] = df['peak']
        df = df[['chr', 'start', 'end', 'name', 'score', 'strand']]
        # print(df.head())
        df.to_csv(f'venn_intersections/HOMER/bed/{sample_name}.bed', sep='\t', index=False, header=False)
        return None

    for filepath in tqdm(glob.iglob(f'{logFC_dir}/*.csv')):
        
        file_string = filepath.split('/')[-1]
        file_string = file_string.replace('.csv', '')
        # Create all the folders 
        if not os.path.exists(f'{motif_dir}/{file_string}'):
            os.makedirs(f'{motif_dir}/{file_string}')
        if not os.path.exists(f'venn_intersections/HOMER/preparsed/{file_string}'):
            os.makedirs(f'venn_intersections/HOMER/preparsed/{file_string}')
        #
        recontstruct_bed(filepath)
        bed_path = f'venn_intersections/HOMER/bed/{file_string}.bed'
        # slopped_path = f'venn_intersections/HOMER/slopped/{file_string}.bed'
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=48:00:00
#SBATCH --output=slurm_outputs/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load homer;
findMotifsGenome.pl {bed_path} {genome} {motif_dir}{file_string} -size {motif_length} -preparsedDir venn_intersections/HOMER/preparsed/{file_string} -p 16 -mis 3 -preparse;
        '''
        f = open(f'venn_intersections/HOMER/submissions/{file_string}.sh', 'w+')
        f.write(cmd)
        f.close()
        os.system(f'sbatch venn_intersections/HOMER/submissions/{file_string}.sh')
    return None

def annotate_genes_near_peaks(merged_peaks, genome_sizes, genome_features, outfile):
    '''
    Annotate genes that are within 100kb of peaks. 
    1. Use bedtools slop to extend peaks by 100kb
    2. intersect the extended peaks with the gene features
    '''    
    cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=15:00:00
#SBATCH --output=slurm_outputs/%j.out
module load bedtools;
bedtools slop -i {merged_peaks} -g {genome_sizes} -b 100000 > merged/{outfile};
bedtools intersect -a merged/{outfile} -b {genome_features} -wa -wb > merged/{outfile}.genes;
'''
    f = open(f'submissions/{outfile.split("/")[-1]}.sh', 'w+')
    f.write(cmd)
    f.close()
    os.system(f'sbatch submissions/{outfile.split("/")[-1]}.sh')

    return None

    # how do I find genes that are within a bed file? 
