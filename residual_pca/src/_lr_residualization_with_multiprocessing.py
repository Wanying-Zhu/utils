# Author: Wanying Zhu
# Modified from linear_regression_with_multiprocessing.py
# Residualize a tabular data file

import pandas as pd
import sys
import argparse
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import datetime
from multiprocessing import Pool
import tqdm

print(datetime.datetime.now().strftime('%Y-%m-%d'))
start = datetime.datetime.now()

# ########## Helper functions ##########
def run_ols(df_data, phenotype, covars, verbose=False,
            indx=None, fn_model=None, fn_residual=None):
    '''
    Run a single linear regression using model: phenotype ~ covars
    Params:
    - phenotype: outcome of the model, such as lipid concentration or gene expression
    - covars: covariates
    - verbose: print summary of the model if True
    - indx=None: index number of current trait
    - fn_model: file name to save model params (such as number of observations, R2, etc.)
    - fn_residual: file name to save residuals
    Return:
    - Write residual and model to output files
    '''
    try:
        predictors = list(set(covars)) # In case there are duplicate covariates
        X = df_data[predictors].copy()
        X['const'] = 1 # Add a constant column for interception
        y = df_data[phenotype]
        model = sm.OLS(y, X, missing='drop')
        results = model.fit()
        # Get model parameters: Number of observations, number of predictors (regressors), R2, adjusted r2
        n, p, r2, r2_adj= results.nobs, results.df_model, results.rsquared, results.rsquared_adj
        with open(fn_model, 'a') as fh_model:
            fh_model.write(f'{phenotype}\t{n}\t{p}\t{r2}\t{r2_adj}\n')
        
        # Save residuals
        with open(fn_residual, 'a') as fh_resid:
            fh_resid.write(phenotype+'\t'+'\t'.join([str(val) for val in results.resid])+'\n')    
        if verbose: print(results.summary())
    except:
        logging.info('# - Ignore column: %s' % phenotype) # Ignore columns that can not be run (usually they are ID columns which do not contain numbers)
    
def merge_files(input_fns, output_fn, header=False):
    '''
    Merge and delete tmp files
    Params:
    - input_fns: an array of tmp file names.
    - output_fn: file name of the final merged result
    - desc: extra text to display before the progress bar
    - header: Whether the tmp files have a header line (default is no header line).
              Will only keep the header line from the first tmp file 
    '''
    fh_out = open(output_fn, 'a')
    for fn in tqdm.tqdm(input_fns):
        try:
            with open(fn) as fh:
                if header: fh.readline() # Skip the header line
                line = fh.readline().strip()
                while line != '':
                    fh_out.write(line + '\n')
                    line = fh.readline().strip()
            os.remove(fn)
        except:
            logging.info('#     - tmp file not found: %s' % fn)
    fh_out.close()

def run_ols_wapped(arguments):
    '''
    Wrap run_ols so it only takes one argument. Necessary for multiprocessing
    '''
    return run_ols(*arguments)

def run_ols_multithreading(df_data, lst_phenotype, covars, output_prefix, threads, verbose=False, fn_residual='', fn_model=''):
    '''
    Run linear regression with multiprocessing.
    Output result into temp files. Merge and remove tmp files once all processes are done
    Params: parameters used by run_ols()
    - df_data: a DataFrame containing phenotypes and covariates
    - lst_phenotype: a list of phenotypes to run OLS
    - output_prefix
    - threads: number of threads for multiprocessing
    - fn_residual, fn_model: file name to save the residuals and model
    '''
    arguments = [] # Create positional arguments to use in run_ols()
    # Store file names of tmp files for merging and cleaning
    tmp_residual, tmp_model = [], []
    for i, phenotype in enumerate(lst_phenotype):
        fn_tmp_residual, fn_tmp_model = f'{output_prefix}.residual.tmp{i}', f'{output_prefix}.model.tmp{i}'
        # print('#'*50, fn_tmp_result)
        tmp_residual.append(fn_tmp_residual)
        tmp_model.append(fn_tmp_model)
        
        # fn_residual_*.tmp* files are intermediate files and will be merged and deleted at the end
        # run_ols(df_data, phenotype, covars, verbose=False, indx=None, fn_model=None, fn_residual=None)
        args_single_run = [df_data, phenotype, covars, verbose, i, fn_tmp_model, fn_tmp_residual]
        arguments.append(args_single_run)
    
    if len(lst_phenotype)>10000:
        chunksize = len(lst_phenotype)//(threads*4) # Use the same chunckszie as multiprocessing.map when list is large
    else: chunksize = 1 # Also the default chuncksize of imap/imap_unordered
    
    with Pool(threads) as p:
        # p.imap_unordered(run_ols_wapped, arguments, chunksize=chunksize)
        for _ in tqdm.tqdm(p.imap_unordered(run_ols_wapped, arguments, chunksize=chunksize), total=len(lst_phenotype)): pass
    
    logging.info('# - Merge and clean tmp files')
    # Process tmp permutation files
    merge_files(tmp_residual, output_fn=fn_residual, header=None)
    merge_files(tmp_model, output_fn=fn_model, header=None)

def load_data(args):
    '''
    Load the input file (and covariate file)
    Params:
    - args: arguemnts from the commandline
    Return:
    - A dataframe of input file and covariate file merged on the ID column
    '''
    if args.threads>1: # Import multiprocessing if needed
        # from multiprocessing import Pool
        if args.threads>os.cpu_count(): # Limit the max number of threads to be the number of cpus
            args.threads = os.cpu_count()

    # ########## Load files ##########
    logging.info('')
    logging.info('# Loading files')
    if args.input_file.endswith('csv'):
        df_data = pd.read_csv(args.input_file)
    else:
        df_data = pd.read_csv(args.input_file, sep='\t') # Assume tab is the delimiter if the input is not a .csv file
    logging.info('# - Input file: (%s,%s)' % df_data.shape)

    # Drop ignored columns to avoid NA, and remove missing values
    for col in args.ignore_cols:
        if col not in args.covars + [args.id_col]:
            try:
                df_data.drop(columns=col, inplace=True)
            except:
                pass
    
    if args.covar_file:
        if args.covar_file.endswith('csv'):
            df_covar = pd.read_csv(args.covar_file)
        else:
            df_covar = pd.read_csv(args.covar_file, sep='\s+') # Assume tab is the delimiter if the input is not a .csv file
        logging.info('# - Covariate file: (%s,%s)' % df_covar.shape)

        # Drop ignored columns to avoid NA, and remove missing values
        for col in args.ignore_cols:
            if col not in args.covars + [args.id_col]:
                try:
                    df_covar.drop(columns=col, inplace=True)
                except:
                    pass
            
        # Merge data with covariates
        df_data = df_data.merge(df_covar, on=args.id_col)
        
    df_data.dropna(inplace=True)
    logging.info('# - Input file cleaned (no NAs): (%s,%s)' % df_data.shape)
    
    # Check if covariates columns exist in the input and covariate files
    for val in args.covars:
        if val not in df_data.columns:
            logging.info('# ERROR: Column not found in data or covariate file: %s' % val)
            exit()
    return df_data

def residualization(args):
    '''
    Residualize the input file
    Params:
    - args: arguemnts from the commandline
    Return:
    - Save residuals into a new file
    '''
    df_data = load_data(args)
    # ########## Run linear regression ##########
    # The model is: phenotype ~ covariates + condition
    logging.info('')
    logging.info('# Run linear regression')
    
    # Save model params (number of observations, R2, etc.)
    fn_model = f'{args.output_path}/{args.output_prefix}.model'
    with open(fn_model, 'w') as fh_model: # Write header line
        fh_model.write('phenotype\tN_observations\tN_predictors\tR2\tAdjusted_R2\n')
        
    # Save residuals to output file
    fn_residual = f'{args.output_path}/{args.output_prefix}.residual'
    with open(fn_residual, 'w') as fh_resid: # Write header line (sample IDs)
        fh_resid.write(f'{args.id_col}\t'+'\t'.join(df_data[args.id_col])+'\n')
    
    lst_phenotype = [] # Create a list of phenotypes to iterate
    for val in df_data.columns:
        if val not in args.covars + args.ignore_cols + [args.id_col]: # Exclude covariates and id columns
            lst_phenotype.append(val)
                
    if args.threads>1: # Multiprocessing
        logging.info('# Run regression with multiprocessing')
        # run_ols_multithreading(df_data, lst_phenotype, covars, output_prefix, verbose=False, fn_residual='', fn_model='')
        run_ols_multithreading(df_data=df_data, lst_phenotype=lst_phenotype, covars=args.covars, output_prefix=args.output_prefix,
                               threads=args.threads, verbose=args.verbose, fn_model=fn_model, fn_residual=fn_residual)
    else: # Regular sequential runs
        for i, phenotype in enumerate(lst_phenotype):
            print(f'\r# - Processing: {i+1}/{len(lst_phenotype)}'+' '*20, end='', flush=True)
            # run_ols(df_data, phenotype, covars, verbose=False, indx=None, fn_model=None, fn_residual=None):
            run_ols(df_data=df_data, phenotype=phenotype, covars=args.covars,
                    verbose=args.verbose, indx=i+1, fn_model=fn_model, fn_residual=fn_residual)
        print(f'\r# - Processing: {i+1}/{len(lst_phenotype)}'+' '*20)

    # Reformat residual file so that rows are samples, columns are features
    logging.info('# Reformat residual file so that rows are samples, columns are features')
    pd.read_csv(fn_residual, sep='\t').T.to_csv(fn_residual, sep='\t', header=False)
    
    # Record running time
    duration = datetime.datetime.now() - start
    duration_d = duration.days
    duration_h = duration.total_seconds()//3600
    duration_m = duration.total_seconds()//60 - duration_h
    duration_s = duration.total_seconds()%60
    msg = f'# - Linear regression finished in {duration_d} days, {duration_h} hours, {duration_m} minutes, {duration_s:.4f} seconds'
    logging.info(msg)
    logging.info('')
    


