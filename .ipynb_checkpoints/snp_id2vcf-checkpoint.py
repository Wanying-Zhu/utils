# Create an input file (VCF) of VEP annotation
# Given a list of SNP IDs (chr:pos:ref:alt), create a VCF4.0 file
# Check here for details: http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#input

'''
Example run:
python /data100t1/home/wanying/lab_code/utils/snp_id2vcf.py \
    --input /data100t1/home/wanying/CCHC/lipidomics/code/post_gwas/test.txt \
    --output test.converted \
    --id_col_name SNP

python /data100t1/home/wanying/lab_code/utils/snp_id2vcf.py \
    --input /data100t1/home/wanying/CCHC/lipidomics/code/post_gwas/test.txt \
    --output test.converted \
    --id_col_name SNP \
    --clean_id_col     # Set --clean_id_col to True
'''

import argparse
import pandas as pd
import os

def process_args():
    '''
    Process arguments and sanity check
    '''
    parser = argparse.ArgumentParser(prog='SNP IDs to VCF file',
                                     description='Create an input file (VCF) of VEP annotation')
    
    parser.add_argument('-i', '--input', type=str,
                        help='Input file should contain one SNP per row. Should contain a column with the name specified by --id_col_name')
    parser.add_argument('-o', '--output_prefix',  type=str, help='Prefix of output VCF')
    parser.add_argument('--id_col_name', type=str, help='Column name of the SNP ID, in the format of chr:pos:ref:alt')
    parser.add_argument('--clean_id_col', action='store_true', help='Default fill the ID column in the output VCF with ".". If true, save the --id_col_name column as ID column')
    # parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    
    # Sanity check
    if not os.path.isfile(args.input):
        print('# ERROR: Input file not found')
        print(f'# - {args.input}')
        exit()
    print('# Arguments used:')
    for arg in vars(args):
        print('# - %s: %s' % (arg, getattr(args, arg)))
    return args


def snp_id2vcf(data, output_prefix, id_col_name, clean_id_col=False):
    '''
    Create an input file (VCF) of VEP annotation
    Given a list of SNP IDs (chr:pos:ref:alt), create a VCF4.0 file
    Check here for details: http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#input
    
    Params:
    - data: a pandas DataFrame object, OR an input file name (string)
            Input file should contain one SNP per row. Should contain a column with the name specified by id_col_name. Other columns are allowed but will not be used
    - output_prefix: Prefix of output VCF
    - id_col_name: Column name of the SNP ID, in the format of chr:pos:ref:alt
    - clean_id_col: Default fill the ID column in the output VCF with ".". If True, save the --id_col_name column as ID column

    Output:
    - Save result to a VCF4.0 with name {output_prefix}.vcf

    Example run:
    python /data100t1/home/wanying/lab_code/utils/snp_id2vcf.py \
        --input /data100t1/home/wanying/CCHC/lipidomics/code/post_gwas/test.txt \
        --output test.converted \
        --id_col_name SNP
    
    python /data100t1/home/wanying/lab_code/utils/snp_id2vcf.py \
        --input /data100t1/home/wanying/CCHC/lipidomics/code/post_gwas/test.txt \
        --output test.converted \
        --id_col_name SNP \
        --clean_id_col     # Set --clean_id_col to True
    '''
    if type(data) is str:
        if data.endswith('csv'): df = pd.read_csv(input_fn)
        else: df = pd.read_csv(input_fn, sep='\t')
    else:
        df = data
            
    # Expected columns in the output file: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    df['#CHROM'] = df[id_col_name].apply(lambda x: x.split(':')[0].upper().split('CHR')[-1]) # Just need the numerical value, remove "CHR" or "chr"
    df['POS'] = df[id_col_name].apply(lambda x: x.split(':')[1])
    df['REF'] = df[id_col_name].apply(lambda x: x.split(':')[2])
    df['ALT'] = df[id_col_name].apply(lambda x: x.split(':')[3])
    
    if clean_id_col: df['ID'] = '.'
    else: df['ID'] = df[id_col_name]
        
    for col in ['QUAL', 'FILTER', 'INFO', 'FORMAT']:
        df[col] = '.'
    
    cols_to_save = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    df[cols_to_save].to_csv(output_prefix+'.vcf', sep='\t', index=False)

if __name__ == "__main__":
    print('# Convert SNP IDs to VCF for VEP annotation')
    args = process_args()
    snp_id2vcf(data=args.input,
               output_prefix=args.output_prefix,
               id_col_name=args.id_col_name,
               clean_id_col=args.clean_id_col)
