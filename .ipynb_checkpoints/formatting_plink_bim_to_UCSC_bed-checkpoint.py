# This code convert plink .bim file to UCSC BED file
# This is necessary for UCSC liftover program (convert between different genome assemblies)
# UCSC BED file has this format:
#       chrom     chromStart      chromEnd      SNP_id(Can be omitted)
#       chr3        12345           12346       rs1111(Can be omitted)

import pandas as pd

# Convert plink .bim file to UCSC BED format
def plink_bim_to_ucsc_bed(input_fn, output_fn):
    output_fh = open(output_fn, 'w')
    with open(input_fn, 'r') as input_fh:
        line = input_fh.readline().strip()
        count = 0
        while line != '':
            lst_tmp = line.split()
            chromosome = 'chr' + lst_tmp[0]
            start_pos = lst_tmp[3]
            end_pos = int(start_pos) + 1

            # Add SNP id (lst_tmp[1]) so that UCSC BED file can be converted back to plink bim file later
            output_line = chromosome+'\t'+start_pos+'\t'+str(end_pos)+'\t'+lst_tmp[1]+'\n'
            output_fh.write(output_line)

            count = count + 1
            if count%20000 == 0: print('.', end='', flush=True)   # Keep console busy
            line = input_fh.readline().strip()
        print('\nDone')

# Convert UCSC BED file back to plink bim file, based on SNP id and old plink .bim file
# (This bed file should be processed by liftover already)
# Parameters:
#   - input_bed_fn: UCSC BED file, it is the converted result from Liftover program
#   - input_bim_fn: The original plink bim file
#   - output_fn: output file name
def ucsc_bed_to_plink_bim(input_bed_fn, input_bim_fn, output_fn):
    print('Read in files ...')
    df_bed = pd.read_csv(input_bed_fn, sep='\t', dtype='str', header=None)
    df_bim = pd.read_csv(input_bim_fn, sep='\t', dtype='str', header=None)
    df_bed.columns = ['chr', 'start', 'end', 'SNP']
    df_bim.columns = ['chr', 'SNP', 'pos', 'coordinate', 'ALT', 'REF']

    print('Converting ...')
    df_merge = df_bim.merge(df_bed, left_on='SNP', right_on='SNP')
    df_merge.dropna(axis=0, inplace=True)

    output_columns = ['chr_x','SNP', 'pos', 'start','ALT','REF']
    df_merge[output_columns].to_csv(output_fn, sep='\t', header=None, index=False)
    print('Done')





# Change input and output file names as needed
input_fn = '/data100t1/home/wanying/hla_analysis/data_and_doc/chr6_TopMed_imputed/chr6_merged_raw.bim'
output_fn = '/data100t1/home/wanying/hla_analysis/data_and_doc/chr6_TopMed_imputed/converted_UCSC_bed.bed'

# plink_bim_to_ucsc_bed(input_fn, output_fn)  # Convert plink .bim file to UCSC BED format
# Use this after processing
# liftOver /data100t1/home/wanying/hla_analysis/data_and_doc/chr6_TopMed_imputed/converted_UCSC_bed.bed hg38ToHg19.over.chain.gz output.bed unlifted.bed

input_bed_fn = '/data100t1/home/wanying/hla_analysis/data_and_doc/chr6_TopMed_imputed/output.bed'
input_bim_fn = '/data100t1/home/wanying/hla_analysis/data_and_doc/chr6_TopMed_imputed/chr6_merged_raw.bim'
output_fn = '/data100t1/home/wanying/hla_analysis/data_and_doc/chr6_TopMed_imputed/chr6_merged_raw_converted_hg19.bim'
ucsc_bed_to_plink_bim(input_bed_fn, input_bim_fn, output_fn)