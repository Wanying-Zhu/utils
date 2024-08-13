'''
Convert to phecode table

Example call:
python count_phecode.py --input_fn /data100t1/share/BioVU/phenos/official_release_0619/Below_PheCodes_20190531.csv \
--output_path result \
--output_prefix result \
--col_names GRID PHEWAS_CODE
'''

import pandas as pd
import numpy as np
import os
import argparse
import sys
import logging

def setup_log(fn_log, mode='w'):
    '''
    Print log message to console and write to a log file.
    Will overwrite existing log file by default
    Params:
    - fn_log: name of the log file
    - mode: writing mode. Change mode='a' for appending
    '''
    # f string is not fully compatible with logging, so use %s for string formatting
    logging.root.handlers = [] # Remove potential handler set up by others (especially in google colab)
    logging.basicConfig(level=logging.DEBUG,
                        handlers=[logging.FileHandler(filename=fn_log, mode=mode),
                                  logging.StreamHandler()], format='%(message)s')


def process_args():
    '''
    Process arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fn', help='Input file name', type=str,
                        default='/data100t1/share/BioVU/phenos/official_release_0619/Below_PheCodes_20190531.csv')
    parser.add_argument('--output_path', type=str, default='./')
    parser.add_argument('--output_prefix', type=str, default='output')
    parser.add_argument('--col_names', nargs='*', default=['GRID', 'PHEWAS_CODE'],
                        help='Column names of ID and phecode (ID column first, phecode column second)')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.output_path): # Create output folder if not exists
        print('# Create output path: ' + args.output_path)
        os.makedirs(args.output_path)

    # Record arguments used
    fn_log = os.path.join(args.output_path, args.output_prefix+'.log')
    setup_log(fn_log, mode='w')
    # logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message
        
    # Record script used
    cmd_used = 'python ' + ' '.join(sys.argv)

    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    
    logging.info('# Arguments used:')
    for arg in vars(args):
        cmd_used += f' --{arg} {getattr(args, arg)}'
        msg = f'# - {arg}: {getattr(args, arg)}'
        logging.info(msg)
    
    return args

args = process_args()

# Process phecode
logging.info('\n# Load phecode file: '+ args.input_fn)
if args.input_fn.endswith('.csv'):
    df_phecode = pd.read_csv(args.input_fn, dtype=str)
else:
    df_phecode = pd.read_csv(args.input_fn, dtype=str, sep='\t')
    
# print(df_phecode.head())
logging.info('# - File size: (%s, %s)' % df_phecode.shape)
df_phecode['CONST'] = 0 # Add a constant column for counting

df_phecode_count = df_phecode.groupby(by=args.col_names).count().reset_index()

lst_dfs, c = [], 0
grids = df_phecode_count[args.col_names[0]].unique()
logging.info('# - Count number of phecodes and reformat')
for phecode, df in df_phecode_count.groupby(args.col_names[1]):
    c += 1
    # PHEWAS_CODE_DESCRIPTION or PHEWAS_DATE contains count of that phecode
    df.rename(columns={'CONST':phecode}, inplace=True)
    df.set_index(keys=args.col_names[0], inplace=True)
    lst_dfs.append(df.reindex(index=grids, fill_value=0)[phecode])
    print(f'\r# - Process {c}    ', end='', flush=True)
df_merged = pd.concat(lst_dfs, axis=1)

output_fn = os.path.join(args.output_path, args.output_prefix)
logging.info('\n# Save to output: ' + output_fn+'.count.txt')
df_merged.sort_index(inplace=True)
df_merged.to_csv(output_fn+'.count.txt', sep='\t')

logging.info('# Save binary to output: ' + output_fn+'.binary.txt')
df_merged[:] = np.where(df_merged>0, 1, 0) # Replace with binary values
df_merged.to_csv(output_fn+'.binary.txt', sep='\t')

print('# DONE')
