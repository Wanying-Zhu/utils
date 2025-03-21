# Author: Wanying Zhu

# TODOs
# - multiprocessing

'''
Example runs in terminal:
To run a single trait:
python linear_regression.py \
--input_file /data100t1/home/gabi_castrodad/gabi_summer_project_files/data/lipid_class.with_covar_pcair_pc.log2.discovery.csv \
--output_path result \
--output_fn result.txt \
--covars AGE_AT_VISIT GENDER PC1 PC2 PC3 PC4 PC5 \
--condition VISIT \
--phenotype Sph

To run all phenotypes provided in the input file, and generate plots:
python linear_regression.py \
--input_file example_data/input.txt \
--output_path output \
--output_fn result.txt \
--covars AGE GENDER PC1 PC2 PC3 \
--condition CONDITION1 \
--ignore_cols ID CONDITION2 \
--permute 1000 \
--plot

* If only need one predictor, use the same one as --condition and --covars

'''

import pandas as pd
import argparse
import logging
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
import os
import datetime
# from multiprocessing import Pool

print(datetime.datetime.now().strftime('%Y-%m-%d'))
start = datetime.datetime.now()

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

def run_ols(df_data, phenotype, covars, condition, fh_output, permute=False, verbose=False,
            fn_permute_pvals='', fn_permute_betas='', fn_permute_stds='', indx=None):
    '''
    Run a single linear regression using model: phenotype ~ covars + condition
    Params:
    - phenotype: outcome of the model, such as lipid concentration or gene expression
    - covars: covariates
    - condition: such as T2D, liver/kidney traits
    - fh_output: output file handle to write result
    - permute=False: Run permutation
    - fn_permute_pvals, fn_permute_betas, fn_permute_stds: output file of permutation results if permute=True
    - verbose: print summary of the model if True
    - indx=None: index number of current trait
    Return:
    - Write pval, std, beta to output file fh_output
    '''
    try:
        predictors = list(set([args.condition] + args.covars)) # In case there are duplicate covariates
        X = df_data[predictors].copy()
        X['const'] = 1 # Constant column
        y = df_data[phenotype]
        model = sm.OLS(y, X, missing='drop')
        results = model.fit()
        pval = results.pvalues[args.condition]
        beta = results.params[args.condition]
        std = results.bse[args.condition]
    
        if not args.permute:
            fh_output.write(f'{phenotype}\t{args.condition}\t{pval}\t{beta}\t{std}\n')
        else: # Permutation!
            pvals, betas, stds = np.zeros(args.permute), np.zeros(args.permute), np.zeros(args.permute) # Store permutation results
            for i in range(args.permute):
                if not indx:
                    print(f'\r# - Permutation run {phenotype}: {i+1}/{args.permute}        ', end='', flush=True)
                else:
                    print(f'\r# - Permutation run ({indx}) {phenotype}: {i+1}/{args.permute}        ', end='', flush=True)
                    
                y_permute = np.random.permutation(df_data[phenotype])
                model_permute = sm.OLS(y_permute, X, missing='drop')
                results_permute = model_permute.fit()
                pvals[i] = results_permute.pvalues[args.condition]
                betas[i] = results_permute.params[args.condition]
                stds[i] = results_permute.bse[args.condition]
            # Calculate the p value from permutation, by comparing permutation coefficients with original coefficient
            pval_permutation = np.sum(np.abs(betas)>=np.abs(beta)) / args.permute
            # pval_permutation = np.sum(pvals<pval) / args.permute # Not recommended!
            fh_output.write(f'{phenotype}\t{args.condition}\t{pval}\t{beta}\t{std}\t{pval_permutation}\n')
            
            # Save permutation p values, betas, stds to separate files
            with open(fn_permute_pvals, 'a') as fh_permute_pvals:
                fh_permute_pvals.write(f'{phenotype}\t{args.condition}\t')
                fh_permute_pvals.write('\t'.join([str(v) for v in pvals])+'\n')
                
            with open(fn_permute_betas, 'a') as fh_permute_betas:
                fh_permute_betas.write(f'{phenotype}\t{args.condition}\t')
                fh_permute_betas.write('\t'.join([str(v) for v in betas])+'\n')
                
            with open(fn_permute_stds, 'a') as fh_permute_stds:
                fh_permute_stds.write(f'{phenotype}\t{args.condition}\t')
                fh_permute_stds.write('\t'.join([str(v) for v in stds])+'\n')         
        if verbose: print(results.summary())
    except:
        logging.info('# - Ignore column: %s' % phenotype) # Ignore columns that can not be run (usually they are ID columns which do not contain numbers)

def merge_files(input_fns, output_fn, header='infer'):
    '''
    Merge and delete tmp files
    Params:
    - input_fns: an array of tmp file names
    - output_fn: file name of the final merged result
    - header: 'infer' or None. Whether the tmp files have a header line.
              Will only keep the header line from the first tmp file 
    '''
    lst_dfs = []
    for fn in input_fns:
        lst_dfs.append(pd.read_csv(fn, sep='\t', header=header))
        os.remove(fn)
    df_all = pd.concat(lst_dfs)
    df_all.to_csv(output_fn, sep='\t', index=False)
    

# ########## Load arguments from commandline ##########
parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str,
                    help='Input file contains phenotype (with or without condition and covariates), in row x col = sample x feature format. Can be a .csv file or tab delimited file')
parser.add_argument('--output_path', type=str, default=None, help='Output directory')
parser.add_argument('--output_fn', type=str, default=None, help='Output (result) file name')
parser.add_argument('--covar_file', type=str, default=None,
                    help='(Optional) Name of the file containing covariates and/or condition if not in the input_file, in row x col = sample x feature format. Can be a .csv file or tab delimited file')
parser.add_argument('--id_col', type=str, default='LABID', help='Shared ID column between the input file and covariate file')
parser.add_argument('--ignore_cols', type=str, nargs='*', default=['RRID', 'LABID'],
                    help='Columns to be ignored from the analysis. Used if no single phenotype is specified')
parser.add_argument('--covars', type=str, nargs='*', help='Covariates of the model, such as sex, age, PCs')
parser.add_argument('--condition', type=str, help='Condition of the model, such as T2D')
parser.add_argument('--phenotype', type=str, default=None,
                    help='Outcome of the model, such as lipid concentration. Run all columns in the input file if None (except covariates and condition)')
parser.add_argument('--permute', default=None, type=int, help='Run permutation')
parser.add_argument('--verbose', action='store_true', help='Print more information. Default value is false')
parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output file if true. Default value is false')
parser.add_argument('--plot', action='store_true', help='Create volcano plot and QQ plot. Default value is false')
parser.add_argument('--threads', default=1, type=int, help='Number of threads for multiprocessing. Default is none (1)')
args = parser.parse_args()

if '.' not in args.output_fn:
    args.output_fn = args.output_fn + '.txt'

if args.plot: # Import code if plotting is needed
    import sys
    sys.path.append('/data100t1/home/wanying/lab_code/utils')
    from plot_lr_result import plot_lr_result
if args.threads>1:
    from multiprocessing import Pool
    if args.threads>os.cpu_count(): # Limit the max number threads to be number of cpus
        args.threads = os.cpu_count()

if not args.output_fn: args.output_fn = args.phenotype + '_' + '_'.join(args.covars) # Create default output file name if not provided
if not args.output_path: args.output_path = './' # Set default output fold to current folder

if not os.path.isdir(args.output_path):
    print('# Create output path: ' + args.output_path)
    os.makedirs(args.output_path) # Create output folder if not exist
    
fn_log = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.log'
setup_log(fn_log, mode='w') # Set up a log file
logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message
logging.info('# Arguments used')
for arg in vars(args):
    msg = f'# - {arg}: {getattr(args, arg)}'
    logging.info(msg)

# ########## Sanity checks ##########
if not os.path.isfile(args.input_file): # Check if input file exists
    msg = '# ERROR: Input file not found: ' + args.input_file
    logging.info(msg)
    exit()
    
if args.covar_file: # If covar_file is provided
    if not os.path.isfile(args.covar_file): # Check if covariate file exists
        msg = '# ERROR: Covarriate file not found: ' + args.covar_file
        logging.info(msg)
        exit()
        
# Check if the output file exists when overwrite=False
if os.path.isfile(f'{args.output_path}/{args.output_fn}') and not args.overwrite:
    logging.info('# Output file exists: %s' % args.output_path+'/'+args.output_fn)
    logging.info('# Skip saving and exit')
    exit()

# ########## Load files ##########
logging.info('')
logging.info('# Loading files')
if args.input_file.endswith('csv'):
    df_data = pd.read_csv(args.input_file)
else:
    df_data = pd.read_csv(args.input_file, sep='\t') # Assume tab is the delimiter if the input is not a .csv file
logging.info('# - Input file: (%s,%s)' % df_data.shape)

if args.covar_file:
    if args.covar_file.endswith('csv'):
        df_covar = pd.read_csv(args.covar_file)
    else:
        df_covar = pd.read_csv(args.covar_file, sep='\t') # Assume tab is the delimiter if the input is not a .csv file
    logging.info('# - Covariate file: (%s,%s)' % df_covar.shape)
    
    # Merge data with covariates
    df_data = df_data.merge(df_covar, on=args.id_col)

# Check if covariates, phenotype, and condition columns exist in the input and covariate files
for val in args.covars + [args.condition]:
    if val not in df_data.columns:
        logging.info('# ERROR: Column not found in data or covariate file: %s' % val)
        exit()
if args.phenotype: # If a single phenotype is provided
    if args.phenotype not in df_data.columns:
        logging.info('# ERROR: Phenotype not found in data or covariate file: %s' % args.phenotype)
        exit()
        
# ########## Run linear regression ##########
# The model is: phenotype ~ covariates + condition
# p value and beta are from the condition
logging.info('')
logging.info('# Run linear regression')
fh_output = open(f'{args.output_path}/{args.output_fn}', 'w')
if not args.permute:
    fh_output.write('phenotype\tcondition\tpval\tbeta\tstd\n') # Write header line
else:
    fh_output.write('phenotype\tcondition\tpval\tbeta\tstd\tpval_permute\n') # Write header line

# Save permutation p values, betas, stds to a separate file if needed
fn_permute_pvals = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.permute_pvals'
fn_permute_betas = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.permute_betas'
fn_permute_stds = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.permute_stds'
# Create empty output to append permutation results
if args.permute:
    for fn in [fn_permute_pvals, fn_permute_betas, fn_permute_stds]:
        with open(fn, 'w') as fh: continue
            
if args.phenotype: # Run a single phenotype
    run_ols(df_data=df_data, phenotype=args.phenotype, covars=args.covars, condition=args.condition, fh_output=fh_output, verbose=args.verbose,
            permute=args.permute, fn_permute_pvals=fn_permute_pvals, fn_permute_betas=fn_permute_betas, fn_permute_stds=fn_permute_stds)
    fh_output.close()
else:
    lst_phenotype = [] # Create a list of phenotypes to iterate
    for val in df_data.columns:
        if val not in args.covars + args.ignore_cols + [args.condition, args.id_col]: # Exclude covariates, condition and other ignored columns
            lst_phenotype.append(val)
                
    if args.threads>1: # Multiprocessing
        fh_output.close() # Not needed for now, so close file handle
        logging.info('# Run regression with multiprocessing')
        with Pool(args.threads) as p:
            # Create positional arguments to use in run_ols()
            arguments = []
            # Store file names of tmp files for merging and cleaning
            tmp_results, tmp_permute_pvals, tmp_permute_betas, tmp_permute_stds = [], [], [], []
            for i, phenotype in enumerate(lst_phenotype):
                fn_tmp_result = f'{args.output_path}/{args.output_fn}.tmp{i}'
                fn_tmp_perm_pvals = fn_permute_pvals+f'.tmp{i}'
                fn_tmp_perm_betas = fn_permute_betas+f'.tmp{i}'
                fn_tmp_perm_stds = fn_permute_stds+f'.tmp{i}'
                tmp_results.append(fn_tmp_result)
                tmp_permute_pvals.append(fn_tmp_perm_pvals)
                tmp_permute_betas.append(fn_tmp_perm_betas)
                tmp_permute_stds.append(fn_tmp_perm_stds)
                
                # fn_permute_*.tmp* files are intermediate files and will be merged and deleted at the end
                fh_output_tmp = open(fn_tmp_result, 'w')
                args_single_run = [df_data, phenotype, args.covars, args.condition, fh_output_tmp, args.permute,
                                   args.verbose, fn_tmp_perm_pvals, fn_tmp_perm_betas, fn_tmp_perm_stds, i]
                arguments.append(args_single_run)
            p.starmap(run_ols, arguments)
        logging.info('# - Merge and clean tmp files')
        merge_files(tmp_results, output_fn=f'{args.output_path}/{args.output_fn}')
        merge_files(tmp_permute_pvals, output_fn=fn_permute_pvals, header=None)
        merge_files(tmp_permute_betas, output_fn=fn_permute_betas, header=None)
        merge_files(tmp_permute_stds, output_fn=fn_permute_stds, header=None)
    
    else: # Regular sequential runs
        for i, phenotype in enumerate(lst_phenotype):
            print(f'\r# - Processing: {i+1}/{len(lst_phenotype)}'+' '*20, end='', flush=True)
            run_ols(df_data=df_data, phenotype=phenotype, covars=args.covars, condition=args.condition, fh_output=fh_output, verbose=args.verbose,
                    permute=args.permute, fn_permute_pvals=fn_permute_pvals, fn_permute_betas=fn_permute_betas, fn_permute_stds=fn_permute_stds, indx=i+1)
        print(f'\r# - Processing: {i+1}/{len(lst_phenotype)}'+' '*20)
        fh_output.close()

# Record running time
duration = datetime.datetime.now() - start
duration_d = duration.days
duration_h = duration.total_seconds()//3600
duration_m = duration.total_seconds()//60 - duration_h
duration_s = duration.total_seconds()%60
msg = f'# - Linear regression finished in {duration_d} days, {duration_h} hours, {duration_m} minutes, {duration_s:.4f} seconds'
logging.info(msg)
logging.info('')

# Perform FDR correction and plot (ignore if only a single was analyzed)
if not args.phenotype:
    logging.info('# Perform FDR correction')
    df_output = pd.read_csv(f'{args.output_path}/{args.output_fn}', sep='\t')
    df_output['pval_fdr'] = fdrcorrection(df_output['pval'])[1]
    df_output.to_csv(f'{args.output_path}/{args.output_fn}', sep='\t', index=False)

    if args.plot:
        logging.info('# Create plots')
        if args.permute:
            adjusted_pval_col=['pval_fdr', 'pval_permute']
        else:
            adjusted_pval_col=['pval_fdr']
        plot_lr_result(df_result=df_output, raw_pval_col='pval',
                       adjusted_pval_col=adjusted_pval_col,
                       beta_col='beta', output_path=args.output_path,
                       output_fn_prefix='.'.join(args.output_fn.split('.')[:-1]))

logging.info('# Done')


