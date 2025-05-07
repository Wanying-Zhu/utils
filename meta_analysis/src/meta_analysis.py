# By Wanying Zhu
# Run meta meta-analysis on multiple regression
# Only performs two-result meta-analysis for now
# Need columns of pval and beta, and sample sizes as scalars
# Modified from ChatGPT script
'''
Example call
1. Print out the help message for each flag:
python ./src/meta_analysis.py --help

2. Run meta analysis
# Basic version
result1=./example_data/result1.txt
result2=./example_data/result2.txt
python ./src/meta_analysis.py \
--input_files ${result1} ${result2} \
--output_path ./example_output \
--output_prefix python_meta_output \
--input_delimiter tab \
--pval_cols pval \
--beta_cols beta \
--marker_cols protein \
--sample_sizes 300 200 \
--overwrite

# Map marker IDs in reuslts 1 and 2 via file --id_mapping_fn, also save extra columns x, y, z in the output
# The ID mapping file is provided with the code in ../data/cchc_proteomics_metaID.txt
result1=./example_data/result1.extra_cols.txt
result2=./example_data/result2.extra_cols.txt
python ./src/meta_analysis_ongoing.py \
--input_files ${result1} ${result2} \
--output_path ./example_output \
--output_prefix python_meta_output.extra_cols \
--id_mapping_fn ./data/cchc_proteomics_metaID.txt \
--id_mapping_cols OID3072 OIDHT meta_id \
--extra_cols_to_keep meta_id OID3072 OIDHT protein \
--input_delimiter tab \
--pval_cols pval \
--beta_cols beta \
--marker_cols OID3072 OIDHT \
--sample_sizes 300 200 \
--overwrite

'''


import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse
import logging
import sys
import os

from _setup_logger_and_parse_args import process_args

# ########## Helper function ##########
def load_results(input_delimiter:list, input_files:list):
    '''
    Load regression results
    
    Params:
    - input_delimiter: delimiters of input files in a list
    - input_files: name of input files in a list

    Return:
    - Two dataframes of the results of regression 1 and regression 2
    '''
    if len(input_delimiter)==0: # Infer delimiter
        if input_files[0].endswith('.csv'):
            df_result1 = pd.read_csv(input_files[0])
        else:
            df_result1 = pd.read_csv(input_files[0], sep='\t')
        if input_files[1].endswith('.csv'):
            df_result2 = pd.read_csv(input_files[1])
        else:
            df_result2 = pd.read_csv(input_files[1], sep='\t')
    elif len(input_delimiter)==1:
        df_result1 = pd.read_csv(input_files[0], sep=input_delimiter[0])
        df_result2 = pd.read_csv(input_files[1], sep=input_delimiter[0])
    else:
        df_result1 = pd.read_csv(input_files[0], sep=input_delimiter[0])
        df_result2 = pd.read_csv(input_files[1], sep=input_delimiter[1])
    return df_result1, df_result2


def id_mapping(df_result, df_id_mapping,
               marker_id_in_result='',
               marker_id_in_mapping='',
               shared_id_col='meta_id',
               drop_cols = []):
    '''
    Map marker ID in the result file using an ID mapping dataframe

    Params:
    - df_result, df_id_mapping: Dataframes of regression result, ID mapping
    - marker_id_in_result, marker_id_in_mapping: marker ID columns to merge on in result and ID mapping files
    - shared_id_col: column name of shared IDs for meta analysis
    - drop_cols: columns to ignore in the result file.
                 If there are extra columns (other than marker ID) with the same names as in mapping file, 
                 then merge will rename those files and cause potentail issues.

    Return:
    - A dataframe with meta_id for meta analysis
    '''
    # Drop columns provided
    # Use suffixes=('_result', None) so is there are any overlapping columns
    # The one in id_mapping file will be used
    if len(drop_cols)==0:
        df_merged = df_result.merge(df_id_mapping,
                                    left_on=marker_id_in_result,
                                    right_on=marker_id_in_mapping,
                                    suffixes=('_result', None)).rename(columns={shared_id_col:'meta_id'})
    else:
        df_merged = (df_result.drop(columns=drop_cols).merge(df_id_mapping,
                                                             left_on=marker_id_in_result,
                                                             right_on=marker_id_in_mapping,
                                                             suffixes=('_result', None))
                                                      .rename(columns={shared_id_col:'meta_id'})
                    )
    return df_merged
    
def meta_single_marker(pvals, signs, sample_sizes):
    '''
    Perform two-result meta-analysis on a single marker, using Weighted Stouffer’s Z-score Meta-Analysis
    Uses p-value and direction of effect, weighted according to sample size (SCHEME SAMPLESIZE in METAL)
    Params:
    - pvals, signs, sample_sizes: pval, sign of beta, sample size of regression 1 and 2
    Return:
    - z_meta, p_meta: z score and pvalue of meta-analysis
    - weight: overall combined weight
    '''
    # Input: p-values, signs of effects, and sample sizes
    # pvals = np.array([0.01, 0.04])
    # signs = np.array([1, 1])  # +1 if beta > 0, -1 if beta < 0
    # ns = np.array([1000, 3000]) # Sample size
    pvals = np.array(pvals)
    signs = np.array(signs)
    sample_sizes = np.array(sample_sizes)
    
    # Step 1: Convert p-values to z-scores with direction
    z_scores = norm.isf(pvals / 2) * signs  # two-sided to z, with sign
    
    # Step 2: Compute weights (sqrt of sample sizes)
    weights = np.sqrt(sample_sizes)
    
    # Step 3: Weighted Z-score
    z_meta = np.sum(weights * z_scores) / np.sqrt(np.sum(weights**2))
    # # the overall combined weight
    # weight = np.sqrt(np.sum(weights**2))
    
    # Step 4: Meta-analysis p-value
    p_meta = 2 * norm.sf(abs(z_meta))
    
    # print(f"Meta-analysis Z-score: {z_meta:.4f}")
    # print(f"Meta-analysis p-value: {p_meta:.4e}")
    return z_meta, p_meta

def meta_all_markers(df_result1, df_result2, all_ids, merged_ids, output_fn):
    '''
    Run meta-analysis using all markers
    Params:
    - df_result1, df_result2: regression result 1 and 2 to be meta-analyzed
    - all_ids: an array of all markers
    - merged_ids: an array of shared markers between two results
    - output_fn: output file to save meta-analyzed data
    - 

    Return:
    - Save results to a file
    '''
    total_rows = len(all_ids)
    n_meta_analized = 0
    # Create output file
    fh_output = open(output_fn, 'w')
    fh_output.write('MarkerName\tWeight\tZscore\tP-value\tDirection\n') # Use the same column headers as METAL
    dict_direction = {-1:'-', 1:'+'} # For directions
    for i, marker in enumerate(all_ids):
        if marker in merged_ids: # Shared marker, can do meta analysis
            mask1, mask2 = df_result1[shared_col1]==marker, df_result2[shared_col2]==marker
            pvals = [df_result1[mask1][pval_col1].iloc[0], df_result2[mask2][pval_col2].iloc[0]]
            signs = [df_result1[mask1]['signs'].iloc[0], df_result2[mask2]['signs'].iloc[0]]
            direction = np.where(np.array(signs)<0, '-', '+')
            sample_sizes = [sample_size1, sample_size2]
            z_meta, p_meta = meta_single_marker(pvals=pvals, signs=signs, sample_sizes=sample_sizes)
            weight = sample_size1 + sample_size2 # Follow METAL output
            fh_output.write(f"{marker}\t{weight}\t{z_meta}\t{p_meta}\t{''.join(direction)}\n")
            n_meta_analized += 1
        else:
            # else do nothing, skip meta analysis and output directly
            if marker in df_result1[shared_col1].values:
                mask = df_result1[shared_col1]==marker
                weight = sample_size1
                pval = df_result1[mask][pval_col1].iloc[0]
                sign = df_result1[mask]['signs'].iloc[0]
                direction = [dict_direction[sign], '?']
            else:
                mask = df_result2[shared_col2]==marker
                weight = sample_size2
                pval = df_result2[mask][pval_col2].iloc[0]
                sign = df_result2[mask]['signs'].iloc[0]
                direction = ['?', dict_direction[sign]]
            # Get z from pval
            z = norm.isf(pval/2) * sign
            fh_output.write(f"{marker}\t{weight}\t{z}\t{pval}\t{''.join(direction)}\n")
        if i%50 ==0:
            print(f'\r# Total markers processed: {i}/{total_rows}    ', end='', flush=True)
    print(f'\r# Total markers processed: {i}/{total_rows}    ')
    
# ########## End of helper function ##########



if __name__=='__main__':
    # ########## Get arguments from terminal ##########
    args = process_args()
    output_fn = os.path.join(args.output_path, args.output_prefix+'.meta_output.txt')
    
    if len(args.input_files)>2:
        logging.info('# Error: This script only perfom two file meta-analysis for now')
        logging.info('# Exit')
        exit()
    
    # ########## Load result files ##########
    df_result1, df_result2 = load_results(input_delimiter=args.input_delimiter,
                                          input_files=args.input_files)
    
    # Get info from arguments
    for val in ['marker_cols', 'pval_cols', 'beta_cols', 'sample_sizes']:
        if len(eval(f'args.{val}'))==1:
            # Use [:-1] to remove 's' in each variable name
            exec(f'{val[:-1]}1 = args.{val}[0]')
            exec(f'{val[:-1]}2 = args.{val}[0]')
        else:
            exec(f'{val[:-1]}1=args.{val}[0]')
            exec(f'{val[:-1]}2 = args.{val}[1]')
    sample_size1, sample_size2 = int(sample_size1), int(sample_size2)


    # ########## Merge with ID mapping file if provided ##########
    if args.id_mapping_fn is not None:
        if args.id_mapping_fn.endswith('.csv'):
            df_id_mapping = pd.read_csv(args.id_mapping_fn)
        else:
            df_id_mapping = pd.read_csv(args.id_mapping_fn, sep='\t')
        if args.id_mapping_cols == []: # Use the column headers in reuslt files and meta_id
            args.id_mapping_cols = [marker_col1, marker_col2, 'meta_id']

        # Ignore extra columns in the mapping file
        cols_to_keep = [val for val in args.extra_cols_to_keep if val in df_id_mapping.columns] + args.id_mapping_cols
        
        # Merge with df_result1 and df_result2
        df_result1 = id_mapping(df_result=df_result1,
                                df_id_mapping=df_id_mapping[list(set(cols_to_keep))],
                                marker_id_in_result=marker_col1,
                                marker_id_in_mapping=args.id_mapping_cols[0],
                                shared_id_col=args.id_mapping_cols[2])
        df_result2 = id_mapping(df_result=df_result2,
                                df_id_mapping=df_id_mapping[list(set(cols_to_keep))],
                                marker_id_in_result=marker_col2,
                                marker_id_in_mapping=args.id_mapping_cols[1],
                                shared_id_col=args.id_mapping_cols[2])
    
    # Merge the two results by shared column, only run meta-analysis on these rows
    if args.id_mapping_fn is not None:
        # Shared column is already renamed to 'meta_id'
        shared_col = 'meta_id'
        marker_col1, marker_col2 = shared_col, shared_col
        
    mask = df_result1[marker_col1].isin(df_result2[marker_col2])
    merged_ids = list(df_result1[mask][marker_col1])
    
    # All markers to process (such as protein IDs)
    all_ids = set(list(df_result1[marker_col1]) + list(df_result2[marker_col2]))

    # Need columns: pval, se, beta and sample size
    df_result1['signs'] = np.sign(df_result1[beta_col1]) 
    df_result2['signs'] = np.sign(df_result2[beta_col2])
    
    total_rows = len(all_ids)
    n_meta_analized = 0
    # Create output file
    fh_output = open(output_fn, 'w')
    
    if len(args.extra_cols_to_keep) != 0:
        header = 'MarkerName\tWeight\tZscore\tP-value\tDirection\t'
        header += '\t'.join(args.extra_cols_to_keep) + '\n'
    else:
        header = 'MarkerName\tWeight\tZscore\tP-value\tDirection\n'
    fh_output.write(header) # Use the same column header as METAL
    dict_direction = {-1:'-', 1:'+'} # For directions
    for i, marker in enumerate(all_ids):
        if marker in merged_ids: # Shared marker, can do meta analysis
            mask1, mask2 = df_result1[marker_col1]==marker, df_result2[marker_col2]==marker
            pvals = [df_result1[mask1][pval_col1].iloc[0], df_result2[mask2][pval_col2].iloc[0]]
            signs = [df_result1[mask1]['signs'].iloc[0], df_result2[mask2]['signs'].iloc[0]]
            direction = np.where(np.array(signs)<0, '-', '+')
            sample_sizes = [sample_size1, sample_size2]
            z_meta, p_meta = meta_single_marker(pvals=pvals, signs=signs, sample_sizes=sample_sizes)
            weight = sample_size1 + sample_size2 # Follow METAL output

            if len(args.extra_cols_to_keep) !=0:
                out_line = f"{marker}\t{weight}\t{z_meta}\t{p_meta}\t{''.join(direction)}\t"
                extral_vals = []
                for col in args.extra_cols_to_keep:
                    extral_vals.append(str(df_result1[mask1][col].iloc[0]))
                out_line += '\t'.join(extral_vals) + '\n'
            else:
                out_line = f"{marker}\t{weight}\t{z_meta}\t{p_meta}\t{''.join(direction)}\n"
            fh_output.write(out_line)
            n_meta_analized += 1
        else:
            # else do nothing, skip meta analysis and output directly
            if marker in df_result1[marker_col1].values:
                mask = df_result1[marker_col1]==marker
                weight = sample_size1
                pval = df_result1[mask][pval_col1].iloc[0]
                sign = df_result1[mask]['signs'].iloc[0]
                direction = [dict_direction[sign], '?']
                if len(args.extra_cols_to_keep) !=0:
                    extral_vals = []
                    for col in args.extra_cols_to_keep:
                        try:
                            extral_vals.append(str(df_result1[mask][col].iloc[0]))
                        except:
                            extral_vals.append('nan')
            else:
                mask = df_result2[marker_col2]==marker
                weight = sample_size2
                pval = df_result2[mask][pval_col2].iloc[0]
                sign = df_result2[mask]['signs'].iloc[0]
                direction = ['?', dict_direction[sign]]
                if len(args.extra_cols_to_keep) !=0:
                    extral_vals = []
                    for col in args.extra_cols_to_keep:
                        try:
                            extral_vals.append(str(df_result2[mask][col].iloc[0]))
                        except:
                            extral_vals.append('nan')
            # Get z from pval
            z = norm.isf(pval/2) * sign

            if len(args.extra_cols_to_keep) !=0:
                out_line = f"{marker}\t{weight}\t{z}\t{pval}\t{''.join(direction)}\t"
                out_line += '\t'.join(extral_vals) + '\n'
            else:
                out_line = f"{marker}\t{weight}\t{z}\t{pval}\t{''.join(direction)}\n"
            fh_output.write(out_line)
        if i%50 ==0:
            print(f'\r# Total markers processed: {i}/{total_rows}    ', end='', flush=True)
    print(f'\r# Total markers processed: {i}/{total_rows}    ')
    logging.info('# Total number of markers meta-analyzed: %s' % n_meta_analized)
    logging.info('# Done')
    fh_output.close()       







