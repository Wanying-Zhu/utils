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
result1=./example_data/result1.txt
result2=./example_data/result2.txt
python ./src/meta_analysis.py \
--input_files ${result1} ${result2} \
--output_path ./example_output \
--output_prefix python_meta_output \
--input_delimiter tab \
--pval_cols pval \
--se_cols se \
--beta_cols beta \
--shared_cols protein \
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
def meta(pvals, signs, sample_sizes):
    '''
    Perform two-result meta-analysis, using Weighted Stouffer’s Z-score Meta-Analysis
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
# ########## End of helper function ##########



# ########## Get arguments from terminal ##########
args = process_args()
output_fn = os.path.join(args.output_path, args.output_prefix+'.meta_output.txt')

if len(args.input_files)>2:
    logging.info('# Error: This script only perfom two file meta-analysis for now')
    logging.info('# Exit')
    exit()

# ########## Load regression results ##########
if len(args.input_delimiter)==0: # Infer delimiter
    if args.input_files[0].endswith('.csv'):
        df_result1 = pd.read_csv(args.input_files[0])
    else:
        df_result1 = pd.read_csv(args.input_files[0], sep='\t')
    if args.input_files[1].endswith('.csv'):
        df_result2 = pd.read_csv(args.input_files[1])
    else:
        df_result2 = pd.read_csv(args.input_files[1], sep='\t')
elif len(args.input_delimiter)==1:
    df_result1 = pd.read_csv(args.input_files[0], sep=args.input_delimiter[0])
    df_result2 = pd.read_csv(args.input_files[1], sep=args.input_delimiter[0])
else:
    df_result1 = pd.read_csv(args.input_files[0], sep=args.input_delimiter[0])
    df_result2 = pd.read_csv(args.input_files[1], sep=args.input_delimiter[1])

# Get info from arguments
for val in ['shared_cols', 'pval_cols', 'se_cols', 'beta_cols', 'sample_sizes']:
    if len(eval(f'args.{val}'))==1:
        # Use [:-1] to remove 's' in each variable name
        exec(f'{val[:-1]}1 = args.{val}[0]')
        exec(f'{val[:-1]}2 = args.{val}[0]')
    else:
        exec(f'{val[:-1]}1=args.{val}[0]')
        exec(f'{val[:-1]}2 = args.{val}[1]')
sample_size1, sample_size2 = int(sample_size1), int(sample_size2)

# Merge the two results by shared column, only run meta-analysis on these rows
mask = df_result1[shared_col1].isin(df_result2[shared_col2])
merged_ids = list(df_result1[mask][shared_col1])

# All markers to process (such as proteins)
all_ids = set(list(df_result1[shared_col1]) + list(df_result2[shared_col2]))

# Need columns: pval, se, beta and sample size
df_result1['signs'] = np.sign(df_result1[beta_col1]) 
df_result2['signs'] = np.sign(df_result2[beta_col2])

total_rows = len(all_ids)
n_meta_analized = 0
# Create output file
fh_output = open(output_fn, 'w')
fh_output.write('MarkerName\tWeight\tZscore\tP-value\tDirection\n') # Use the same column header as METAL
dict_direction = {-1:'-', 1:'+'} # For direction
for i, marker in enumerate(all_ids):
    if marker in merged_ids: # Shared marker, can do meta analysis
        mask1, mask2 = df_result1[shared_col1]==marker, df_result2[shared_col2]==marker
        pvals = [df_result1[mask1][pval_col1].iloc[0], df_result2[mask2][pval_col2].iloc[0]]
        signs = [df_result1[mask1]['signs'].iloc[0], df_result2[mask2]['signs'].iloc[0]]
        direction = np.where(np.array(signs)<0, '-', '+')
        sample_sizes = [sample_size1, sample_size2]
        z_meta, p_meta = meta(pvals=pvals, signs=signs, sample_sizes=sample_sizes)
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
logging.info('# Total number of markers meta-analyzed: %s' % n_meta_analized)
logging.info('# Done')
fh_output.close()       







