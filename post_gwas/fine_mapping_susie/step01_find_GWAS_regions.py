'''
Find GWAS region based on lead SNP (SNP with the smallest P value)
Only use distance to look for SNPs.
Merge regions that overlap.
Check the --clump function in plink if need to consider LD

Example call:

python /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/code/utils/01_find_GWAS_regions.py \
--output_path /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/gwas_regions \
--output_prefix test \
--gwas_result "/data100t1/home/wanying/CCHC/lipidomics/output/lip_species_GWAS_noadj_BMI_AGE2/TG-O-54:4-_[NL-18:2].fastGWA" \
--delim tab \
--pval_threshold 5e-8 \
--window_size 1000000 \
--colname_chr CHR \
--colname_pos POS \
--colname_beta BETA \
--colname_se SE \
--get_zscore
'''

import pandas as pd
import numpy as np
import logging
# from step02_merge_overlapping_regions import merge_regions

import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils')
from setup_logger_and_parse_args import process_args

# ########## Helper functions ##########
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

# Merge regions if they overlap
# Load GWAS regions detected
def find_pairs_of_regions(df_regions,
                          colname_id='SNP',
                          colname_index='region_index',
                          colname_pos='POS',
                          colname_chr='CHR'):
    '''
    Find pairs of indices of regions to be merged. Eg. (1,2), (3,4), (4,5)
    It becomes a graph, with each region as a node and a pair of indices as edge
    Param:
    - df_regions: a DataFrame contains regions and indices to each region
    - colname_id: id column to find duplicate snps. Cannot just use pos and chr due to multiallelic sites
    - colname_index: column header of the region index
    - colname_pos: column header of the position
    - colname_chr: column header of the chromosome
    Return:
    - n_merged: number of regions merged (due to overlap)
    - df_regions_merged: cleaned dataframe
    '''
    # Find duplicate SNPs
    mask = df_regions[colname_id].duplicated(keep=False)
    df_duplicate = df_regions[mask].copy()
    
    # Store pairs (or more than 2) of regions need to be merged
    regions_to_merge = set()
    for snp, df in df_duplicate.groupby(colname_id):
        inx_pair = df[colname_index].unique().tolist()
        # Sort the indices in the pair so that the smaller value comes first
        inx_pair.sort()
        regions_to_merge.add(tuple(inx_pair)) # Set only accepts hashable items
    return regions_to_merge
    
# Merge GWAS regions if they overlap
def find_connected_nodes(edges, node0):
    '''
    Given a collection of edges, find all nodes connected to node0
    (Use width first search)
    '''
    edges = edges.copy() # Deep copy to modify without affecting original list
    connceted_nodes = []
    edges_to_pop = [] # track edges to be removed from next iteration
    # Find directly connected nodes first
    for edge in edges:
        if node0 in edge:
            edges_to_pop.append(edge)
            # Find a node directly connected to node0 first
            if node0==edge[0]:
                node1 = edge[1]
            else:
                node1 = edge[0]
            connceted_nodes.append(node1)
            
    for e in edges_to_pop:
        edges.pop(edges.index(e))
    
    # Recursively find nodes directly connected to other nodes
    other_connected_nodes =[]
    for node in connceted_nodes:
        other_connected_nodes += find_connected_nodes(edges, node)
        
    # Remove duplicates and sort
    result = list(set([node0]+connceted_nodes+other_connected_nodes))
    result.sort()
    return result 

def merge_regions(df_regions,
                  colname_id='SNP',
                  colname_index='region_index',
                  colname_pos='POS',
                  colname_chr='CHR'):
    '''
    According to FINEGENE fine mapping pipeline, merge regions if they overlap.
    Process index pairs and merge regions if they overlap
    Eg: Input index pairs with (1,2), (2,3), (4,5), and output (1,2,3), (4,5).
    Label region 1, 2 and 3 with new label "1,2,3"
    
    Param:
    - df_regions: a DataFrame contains regions and indices to each region
    - colname_index: column header of the region index
    - colname_pos: column header of the position
    - colname_chr: column header of the chromosome
    Return:
    - output: Regions with overlapping regions merged
    '''
    logging.info('\n# Get regions to be merged, create new labels')
    regions_to_merge = list(find_pairs_of_regions(df_regions=df_regions,
                                                  colname_id=colname_id,
                                                  colname_index=colname_index,
                                                  colname_pos=colname_pos,
                                                  colname_chr=colname_chr))
    
    lst_nodes = df_regions[colname_index].unique().tolist() # Get nodes (ie. index of each region)
    lst_nodes.sort()
    dict_new_reion_index = dict() # New label to reach region
    for node in lst_nodes:
        dict_new_reion_index[node] = find_connected_nodes(edges=regions_to_merge, node0=node)

    logging.info('# Merge regions and relabel')
    df_regions.drop_duplicates(subset=[colname_id], inplace=True)
    df_regions['merged_region_index'] = df_regions[colname_index].apply(lambda x: ','.join([str(v) for v in dict_new_reion_index[x]]))
    
    # number of regions after merging
    n_after_merged = len(df_regions['merged_region_index'].unique())
    return n_after_merged, df_regions
    
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
                         'help': 'Window size to take around a lead SNP in bp. Default is 1000000=1Mb'},
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
    df_gwas_result = pd.read_csv(args.gwas_result, sep=args.delim,
                                 dtype={args.colname_pos: int, args.colname_pval: float})
    
    logging.info('# - shape (%s, %s)' % df_gwas_result.shape)
    
    # Filter significant SNPs
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
                                                          colname_id='SNP',
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


    
    











