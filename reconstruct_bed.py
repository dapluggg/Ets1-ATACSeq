import sys 
import pandas as pd

deseq2_output = sys.argv[1]    
    
def recontstruct_bed(logFC_file):
    '''
    Reconstruction of bed file from DESeq2 output.
    '''
    sample_name = logFC_file.split('/')[-1]
    sample_name = sample_name.replace('.csv', '')
    df = pd.read_csv(logFC_file)
    df.rename( columns={'Unnamed: 0':'peak'}, inplace=True )
    df['chr'] = df['peak'].str.split('_').str[0] + '_' + df['peak'].str.split('_').str[1]
    df['start'] = df['peak'].str.split('_').str[2]
    df['end'] = df['peak'].str.split('_').str[3]
    df['strand'] = '+'
    df['name'] = df['peak']
    df['score'] = 0
    df = df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    df.to_csv(f'annotated_logFC/bed/{sample_name}.bed', sep='\t', index=False, header=False)
    return None

if __name__ == '__main__':
    recontstruct_bed(deseq2_output)