'''
Find GWAS region based on lead SNP (SNP with the smallest P value)
Use region_index in the output file to avoid overlap regions
Only use distance to look for SNPs.
Check the --clump function in plink if need to consider LD

Example call:

python /data100t1/home/wanying/lab_code/utils/post_gwas/S01_find_signal_regions/01_find_GWAS_regions.py \
--output_path /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/gwas_regions \
--output_prefix test \
--gwas_result "/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_noadj_BMI_AGE2/TG-O-54:4-_[NL-18:2].fastGWA" \
--delim tab \
--pval_threshold 5e-8 \
--window_size 1000000 \
--colname_chr CHR
'''

import pandas as pd
import numpy as np
import logging
from step02_merge_overlapping_regions import merge_regions

import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils')
from setup_logger_and_parse_args import process_args

def find_a_single_region(df_lead_snps, df_gwas_result, count):
    '''
    Takes a list of lead SNPs, find the SNP with the smallest pval
    Params
    - df_lead_snps: a dataframe of lead SNPs filtered by pval threshold from the GWAS result
    - df_gwas_result: a dataframe of GWAS result
    Return
    - df_region: a dataframe of GWAS region
    '''
    top_snp = df_lead_snps[df_lead_snps[args.colname_pval]==df_lead_snps[args.colname_pval].min()].iloc[0, :]
    chromosome, position = top_snp[args.colname_chr], top_snp[args.colname_pos]
    pos_mask = (df_gwas_result[args.colname_pos]>=position-args.window_size/2) & (df_gwas_result[args.colname_pos]<=position+args.window_size/2)
    region_mask = (df_gwas_result[args.colname_chr]==chromosome) & pos_mask
    df_region = df_gwas_result[region_mask].copy()
    df_region['lead_snp'] = 0 # Mark lead SNP as 1
    mask_lead_snp = df_region[args.colname_pos] == position
    df_region.loc[mask_lead_snp, 'lead_snp'] = 1
    df_region['region_index'] = count # Track number of regions
    return chromosome, position, df_region

if __name__ == "__main__":
    # Already has --output_path and --output_prefix
    args = process_args(True,
                        {'flag_name':'--gwas_result', 'type': 'str',
                         'help': 'GWAS result to be processed. Assume tab as delimiter'},
                        {'flag_name':'--delim', 'default':'tab', 'choices':['tab', 'space', ','], 'type': 'str',
                         'help': 'Delimiter in the GWAS result file. Default is tab'},
                        {'flag_name':'--pval_threshold', 'default': 5e-8, 'type': 'float',
                         'help': 'P value threshold to filter when looking for a lead SNP'},
                        {'flag_name':'--window_size', 'default': 1000000, 'type': 'int',
                         'help': 'Window size to take around a lead SNP in bp (pos +/- window_size/2). Default is 1000000=1Mb'},
                        {'flag_name':'--colname_id', 'default': 'SNP', 'type': 'str',
                         'help': 'Column name of the SNP IDs. Should not have duplicates by this column, so use chr:pos:ref:alt instead of RSIDs'},
                        {'flag_name':'--colname_chr', 'default': 'CHR', 'type': 'str',
                         'help': 'Column name of the chromosome number in the GWAS result'},
                        {'flag_name':'--colname_pval', 'default': 'P', 'type': 'str',
                         'help': 'Column name of the p value in the GWAS result'},
                        {'flag_name':'--colname_pos', 'default': 'POS', 'type': 'str',
                         'help': 'Column name of the position in the GWAS result'},
                        {'flag_name':'--extra_pval_threshold', 'default': -1, 'type': 'float',
                         'help': 'Extra p value threshold to filter SNPs in the region. Default -1=no filtering and output all SNPs'},
                        {'flag_name':'--get_zscore', 'action': 'store_true',
                         'help': 'Calculate z score as beta/se'},
                        {'flag_name':'--colname_beta', 'default': 'BETA', 'type': 'str',
                         'help': 'Column name of the beta in the GWAS result. Only used when --get_zscore is true'},
                        {'flag_name':'--colname_se', 'default': 'SE', 'type': 'str',
                         'help': 'Column name of the standard error of beta in the GWAS result. Only used when --get_zscore is true'})
    
    if args.delim == 'tab':
        args.delim = '\t'
    elif args.delim == 'space':
        args.delim = ' '
    
    logging.info('\n# Load GWAS result')
    # Use (Nullable integer data type) Int64 or Int32 instead of int incase there are missing values
    df_gwas_result = pd.read_csv(args.gwas_result, sep=args.delim,
                                 dtype={args.colname_pos: 'Int64', args.colname_pval: float})
    
    # Drop duplicate variants to avoid issue (IndexError: tuple index out of range)
    # Should not have dupliates if using chr:pos:ref:alt format
    # But still drop duplicates just in case (rsID may have duplicates)
    df_gwas_result.drop_duplicates(subset=[args.colname_pval], inplace=True)
    
    logging.info('# - shape (%s, %s)' % df_gwas_result.shape)
    
    # Filter significant SNPs (lead SNPs)
    df_lead_snps = df_gwas_result[df_gwas_result[args.colname_pval]<=args.pval_threshold].copy()
    count, lst_regions = 0, []
    
    # Iteratively find a region
    while len(df_lead_snps)>0:
        count += 1 # track number of regions
        # Take the first row in case more SNPs share the same P values
        chromosome, position, df_region = find_a_single_region(df_lead_snps, df_gwas_result, count)
        if len(df_region)>0:
            lst_regions.append(df_region)
        
        # Remove the top SNP and any other lead SNPs in the window, and move on to the next one by pvalue
        # Keep non-lead SNPs, since overlap regions are merged later
        mask_keep = (df_lead_snps[args.colname_pos]<position-args.window_size/2) | (df_lead_snps[args.colname_pos]>position+args.window_size/2)
        df_lead_snps = df_lead_snps[mask_keep].copy()
        
        print(f'\r# Process region {count}    ', end='', flush=True)
    print()

    if len(lst_regions) > 0:
        df_all_regions = pd.concat(lst_regions)
        df_all_regions.sort_values(by=[args.colname_chr, 'region_index', args.colname_pos], inplace=True)
        # output_fn = f'{args.output_path}/{args.output_prefix}.regions'
        # df_all_regions.to_csv(output_fn, sep='\t', index=False)

        # Merge regions
        logging.info('# Merge overlapping regions and relabel region index')
        output_fn = f'{args.output_path}/{args.output_prefix}.regions.overlap_merged'
        
        n_after_merged, df_regions_merged = merge_regions(df_regions=df_all_regions,
                                                          colname_id=args.colname_id,
                                                          colname_index='region_index',
                                                          colname_pos=args.colname_pval,
                                                          colname_chr=args.colname_chr)
        
        logging.info('# - N regions merged: %s' % n_after_merged)
        if args.extra_pval_threshold !=-1:
            # if need to filter SNPs by p value
            logging.info('# Filter SNPs by extra p values threshold %s' % args.extra_pval_threshold)
            df_regions_merged = df_regions_merged[df_regions_merged[args.colname_pval]<args.extra_pval_threshold].copy()

        if args.get_zscore:
            logging.info('# Calcualte z score (beta/se)')
            df_regions_merged['z_score'] = df_regions_merged[args.colname_beta]/df_regions_merged[args.colname_se]
            
        logging.info('# Save regions to file')
        df_regions_merged.to_csv(output_fn, sep='\t', index=False)
    else: # No region found
        logging.info('# No region to output')
      
    logging.info('\n# DONE')


    
    











