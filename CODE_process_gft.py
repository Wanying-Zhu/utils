# Create cleaned file from Ensemble gtf
'''
python CODE_process_gft.py \
--input_gtf Homo_sapiens.GRCh38.112.gtf.gz \
--output_fn Homo_sapiens.GRCh38.112.reformat \
--chunksize 5000 \
--gzip
'''
import pandas as pd
import argparse
import os
import numpy as np
import subprocess
import datetime
import logging
import tqdm # Might need to install first

# ########## Helper functions ##########
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
                                  logging.StreamHandler()], format='%(name)s - %(levelname)s - %(message)s')
def parse_args():
    '''
    Parse arguments from commendline
    '''
    # ########## Load arguments from commandline ##########
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_gtf', type=str,
                        help='Input gtf file, for example http://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/')
    parser.add_argument('--output_path', type=str, default='./', help='Output directory. Default is current working directgory')
    parser.add_argument('--output_fn', type=str, default=None, help='Output (result) file name')
    parser.add_argument('--chunksize', type=int, default=None, help='Chunk size to process the gtf')
    parser.add_argument('--gzip', action='store_true', help='Gzip output file output if true. Default value is false')
    # parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output file if true. Default value is false')
    args = parser.parse_args()
    
    if not args.output_fn: # Create default output file name if not provided
        args.output_fn = '.'.join(os.path.split(args.input_gtf)[-1].split('.')[:-1]) + '.reformat.txt'
    
    if not os.path.isdir(args.output_path):
        print('# Create output path: ' + args.output_path)
        os.makedirs(args.output_path) # Create output folder if not exist
    return args

def load_gtf(in_fn):
    '''
    Load gtf file downloaded from Ensembl
    Params:
    - in_fn: input file name
    Return:
    - A dataframe of gtf
    '''
    if in_fn.endswith('gz'):
        df = pd.read_csv(in_fn, sep='\t', compression='gzip', comment='#', dtype=str)
    else:
        df = pd.read_csv(in_fn, sep='\t', comment='#', dtype=str)
    # Rename columns
    # the first 8 columns in GTF are the same as in GFF
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'extra']
    df.columns = col_names
    logging.info('# Load gtf file: row*col = %s*%s' % df.shape)
    return df

def load_gtf_by_chunk(in_fn, chunksize=5000):
    '''
    Load gtf file downloaded from Ensembl
    Params:
    - in_fn: input file name
    - chunksize: read in the GTF by chunksize
    Return:
    - TextFileReader object for iteration
    '''
    if in_fn.endswith('gz'):
        data = pd.read_csv(in_fn, sep='\t', compression='gzip', comment='#', dtype=str, chunksize=chunksize)
    else:
        data = pd.read_csv(in_fn, sep='\t', comment='#', dtype=str, chunksize=chunksize)
    logging.info('# Load gtf file by chunk')
    return data

def split_info(val):
    '''
    Split on a single value from last column of gtf and to get each field
    Example value: gene_id "ENSG00000228037"; gene_version "1"; gene_source "havana"; gene_biotype "lncRNA";
    Return a DataFrame with one row, such as:
        gene_id            gene_version    gene_source    gene_biotype
        ENSG00000228037    1               havana         lncRNA
    '''
    dict_fields = dict()
    tmp_list = val.split(';')
    for field in tmp_list:
        try:
            k, v = field.strip().split()
            v = v.strip('"')
            dict_fields[k] = [v]
        except:
            continue
    return pd.DataFrame(dict_fields)

def process_df(df, index=None):
    '''
    Process a dataframe by split fields in the last column,
    then add those field as new columns to the original dataframe
    Params:
    - df: a dataframe in GFT format
    - index: index of current chunk
    Return:
    Modified dataframe
    '''
    lst_dfs = []
    if index is not None:
        logging.info('# Process chunk %s' % index)
    else:
        logging.info('# Split fields in the last column')
    for val in tqdm.tqdm(df.iloc[:, -1]): # Split values in the last column
        lst_dfs.append(split_info(val))
    df_from_last_col = pd.concat(lst_dfs).reset_index().drop(columns='index')
    df_cleaned = pd.concat([df.iloc[:, :-1], df_from_last_col], axis=1)
    logging.info('# - Reformatted file: (%s,%s)' % df_cleaned.shape)
    return df_cleaned

def process_by_chunk(in_fn, chunksize):
    '''
    Split by chunk
    '''
    reader = load_gtf_by_chunk(in_fn, chunksize=chunksize)
    count = 0
    df_lst = [] # Store processed chunks
    for df in reader:
        count += 1
        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'extra']
        df_lst.append(process_df(df, index=count))  
    reader.close()
    df_all = pd.concat(df_lst).reset_index().drop(columns='index')
    return df_all

if __name__ == "__main__":
    args = parse_args()
    fn_log = os.path.join(args.output_path, '.'.join(args.output_fn.split('.')[:-1]))+'.log'
    setup_log(fn_log=fn_log)
    
    start = datetime.datetime.now()
    logging.info(start.strftime('%Y-%m-%d'))
    
    logging.info('# Arguments used')
    for arg in vars(args):
        msg = f'# - {arg}: {getattr(args, arg)}'
        logging.info(msg)
        
    logging.info('# Load GTF')
    if not args.chunksize:
        df = load_gtf(args.input_gtf)
        df_cleaned = process_df(df)
    else: # Process by chunksize
        df_cleaned = process_by_chunk(in_fn=args.input_gtf, chunksize=args.chunksize)
        
    
    logging.info('# Final reformatted file: (%s,%s)' % df_cleaned.shape)
    logging.info('# Save file')
    output_fn = os.path.join(args.output_path, args.output_fn)
    df_cleaned.to_csv(output_fn, sep='\t', index=False, na_rep='NA')
    
    if args.gzip:
        logging.info('# Gzip output')
        subprocess.run(['gzip', output_fn])
    logging.info('# DONE in %.2f min' % ((datetime.datetime.now()-start).total_seconds()/60))
    
    

    
    
