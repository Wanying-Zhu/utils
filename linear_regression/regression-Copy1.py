# Author: Wanying Zhu
# Run OLS or linear mixed model

'''
Example runs in terminal:
To run a single trait:
python regression.py \
--input_file example_data/input.csv \
--output_path output \
--output_fn result_single_lipid.txt \
--covars AGE GENDER PC1 PC2 PC3 \
--condition CONDITION1 \
--phenotype LIPID4 \
--ignore_cols ID \
--id_col ID \
--overwrite \
--get_residual \
--verbose

To run all phenotypes provided in the input file, and generate plots,
python regression.py \
--input_file example_data/input.csv \
--output_path output_v2 \
--output_fn result.txt \
--covars AGE GENDER PC1 PC2 PC3 \
--plot \
--condition CONDITION1 \
--model mixed \
--group_col RRID \
--ignore_cols ID CONDITION2 \
--id_col ID \
--get_residual \
--overwrite \
--permute 100

* If only need one predictor, use the same one as --condition and --covars

'''
# TODO: allow additional covar file (such as list of sample IDs, hidden covariates)
# parser.add_argument('--additional_covar_fn', type=str, help='File name of additional covariates')

import pandas as pd
import sys
import argparse
import logging
import os
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
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
                                  logging.StreamHandler()], format='%(message)s')

def process_args():
    '''
    Process arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str,
                        help='Input file contains phenotype (with or without condition and covariates), in row x col = sample x feature format. Can be a .csv file or tab delimited file')
    parser.add_argument('--output_path', type=str, default='./', help='Output directory. Default is current working directgory')
    parser.add_argument('--output_fn', type=str, default=None, help='Output (result) file name')
    parser.add_argument('--covar_file', type=str, default=None,
                        help='(Optional) Name of the file containing covariates and/or condition if not in the input_file, in row x col = sample x feature format. Can be a .csv file or tab delimited file')
    parser.add_argument('--id_col', type=str, default='LABID',
                        help='Shared ID column between the input file and covariate file. Not used if run a single phenotype')
    parser.add_argument('--ignore_cols', type=str, nargs='*', default=['RRID', 'LABID'],
                        help='Columns to be ignored from the analysis. Used if no single phenotype is specified')
    parser.add_argument('--covars', type=str, nargs='*', help='Covariates of the model, such as sex, age, PCs')
    parser.add_argument('--condition', type=str, help='Condition of the model, such as T2D')
    parser.add_argument('--phenotype', type=str, default=None,
                        help='Outcome of the model, such as lipid concentration. Run all columns in the input file if None (except covariates and condition)')
    parser.add_argument('--model', default='ols', type=str, choices=['ols', 'mixed'],
                        help='ols:linear regression, mixed:linear mixed model')
    parser.add_argument('--permute', default=None, type=int, help='Run permutation')
    
    parser.add_argument('--group_col', default=None, type=str,
                        help='Column name of groups for linear mixed model')
    parser.add_argument('--verbose', action='store_true', help='Print more information. Default value is false')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output file if true. Default value is false')
    parser.add_argument('--plot', action='store_true', help='Create volcano plot and QQ plot. Default value is false')
    parser.add_argument('--get_residual', action='store_true', help='Save redisual to output')
    args = parser.parse_args()
    
    if '.' not in args.output_fn:
        args.output_fn = args.output_fn + '.txt'
        
    if not os.path.isdir(args.output_path):
        print('# Create output path: ' + args.output_path)
        os.makedirs(args.output_path) # Create output folder if not exist
    
    # Set up logger file
    fn_log = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.log'
    setup_log(fn_log, mode='w') # Set up a log file
    logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message

    # Record command line call used
    cmd_used = 'python ' + ' '.join(sys.argv)
    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    
    logging.info('# Arguments used:')
    for arg in vars(args):
        msg = f'# - {arg}: {getattr(args, arg)}'
        logging.info(msg)   
    return args


def run_single_regression(df_data, phenotype, covars, condition, fn_output, permute=False, verbose=False,
                   fn_permute_pvals='output.permute_pvals', fn_permute_betas='output.permute_betas',
                   fn_permute_stds='output.permute_stds', indx=None, fn_model='output.model',
                   get_residual=False, fn_residual='', model_type='ols'):
    '''
    Run a single linear regression using model: phenotype ~ covars + condition
    Or linear mixed model with random effect of intercept
    Params:
    - phenotype: outcome of the model, such as lipid concentration or gene expression
    - covars: covariates
    - condition: such as T2D, liver/kidney traits
    - fn_output: output file name to write the result
    - permute=False: Run permutation
    - fn_permute_pvals, fn_permute_betas, fn_permute_stds: output file of permutation results if permute=True
    - verbose: print summary of the model if True
    - indx=None: index number of current trait
    - fn_model: file name to save model params (such as number of observations, R2, etc.)
    - get_residual: True: save residuals to output file
    - fn_residual: file name to save residuals
    - model: "ols" for linear regression or "mixed" for linear mixed model
    Return:
    - Write pval, std, beta to output file fh_output
    '''
    # try:
    predictors = list(set([args.condition] + args.covars)) # In case there are duplicate covariates
    X = df_data[predictors].copy()
    X['const'] = 1 # Add a constant column for interception
    y = df_data[phenotype]
    if model_type=='ols':
        model = sm.OLS(y, X, missing='drop')
    elif model_type=='mixed':
        # exog_re=X['const'] can be skipped (random intercept is default setting)
        model = sm.MixedLM(y, X, groups=df_data[args.group_col], exog_re=X['const'], missing='drop')
    results = model.fit()
    pval = results.pvalues[args.condition]
    beta = results.params[args.condition]
    std = results.bse[args.condition]
    # Get model parameters: Number of observations, number of predictors (regressors), R2, adjusted r2
    n, p, r2, r2_adj= results.nobs, results.df_model, results.rsquared, results.rsquared_adj
    with open(fn_model, 'a') as fh_model:
        fh_model.write(f'{phenotype}\t{n}\t{p}\t{r2}\t{r2_adj}\n')
                    
    # Get info of the model
    if get_residual: # Save residuals if needed 
        with open(fn_residual, 'a') as fh_resid:
            fh_resid.write(phenotype+'\t'+'\t'.join([str(val) for val in results.resid])+'\n')
    
    with open(fn_output, 'a') as fh_output:
        if not args.permute:
            fh_output.write(f'{phenotype}\t{args.condition}\t{pval}\t{beta}\t{std}\n')
        else: # Permutation!
            pvals, betas, stds = np.zeros(args.permute), np.zeros(args.permute), np.zeros(args.permute) # Store permutation results
            for i in range(args.permute):
                if not indx:
                    if verbose: print(f'\r# - Permutation run {phenotype}: {i+1}/{args.permute}        ', end='', flush=True)
                else:
                    if verbose: print(f'\r# - Permutation run ({indx}) {phenotype}: {i+1}/{args.permute}        ', end='', flush=True)
                    
                y_permute = np.random.permutation(df_data[phenotype])
                if model_type=='ols':
                    model_permute = sm.OLS(y_permute, X, missing='drop')
                elif model_type=='mixed':
                    model_permute = sm.MixedLM(y_permute, X, missing='drop', groups=df_data[args.group_col], exog_re=X['const'])
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
    # except:
    #     logging.info('# - Ignore column: %s' % phenotype) # Ignore columns that can not be run (usually they are ID columns which do not contain numbers)


# ########## End of helper functions ##########

if __name__=='__main__':
    # Load arguments from commandline
    args = process_args()
    
    if args.plot: # Import code if plotting is needed
        import sys
        sys.path.append('/data100t1/home/wanying/lab_code/utils')
        from plot_lr_result import plot_lr_result
    
    if not args.output_fn: args.output_fn = args.phenotype + '_' + '_'.join(args.covars) # Create default output file name if not provided
    
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
    # Check if group column exists if doing mixed effect model
    if args.model == 'mixed':
        if args.group_col not in df_data.columns:
            logging.info('# ERROR: Group column not found for mixed effect model: %s' % args.group_col)
            exit()
            
    # ########## Run linear regression ##########
    # The model is: phenotype ~ covariates + condition
    # p value and beta are from the condition
    logging.info('')
    logging.info('# Run regression')
    fn_output = f'{args.output_path}/{args.output_fn}'
    
    # Save model params (number of observations, R2, etc.)
    fn_model = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.model'
    with open(fn_model, 'w') as fh_model: # Write header line
        fh_model.write('phenotype\tN_observations\tN_predictors\tR2\tAdjusted_R2\n')
        
    # Save residual if needed
    fn_residual = f'{args.output_path}/{".".join(args.output_fn.split(".")[:-1])}.residual'
    if args.get_residual: # Save residuals if needed
        with open(fn_residual, 'w') as fh_resid: # Write header line (sample IDs)
            fh_resid.write('phenotype\t'+'\t'.join(df_data[args.id_col])+'\n')
    
    with open(fn_output, 'w') as fh_output:
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
        run_single_regression(df_data=df_data, phenotype=args.phenotype, covars=args.covars,
                              condition=args.condition, fn_output=fn_output, verbose=args.verbose,
                              permute=args.permute, model_type=args.model, fn_permute_pvals=fn_permute_pvals,
                              fn_permute_betas=fn_permute_betas, fn_permute_stds=fn_permute_stds,
                              fn_model=fn_model, get_residual=args.get_residual, fn_residual=fn_residual)
    else:
        lst_phenotype = [] # Create a list of phenotypes to iterate
        for val in df_data.columns:
            if val not in args.covars + args.ignore_cols + [args.condition, args.id_col]: # Exclude covariates and condition columns
                lst_phenotype.append(val)
                    
        # Regular sequential runs
        for i, phenotype in enumerate(lst_phenotype):
            print(f'\r# - Processing: {i+1}/{len(lst_phenotype)}'+' '*20, end='', flush=True)
            run_single_regression(df_data=df_data, phenotype=phenotype, covars=args.covars,
                                  condition=args.condition, fn_output=fn_output, verbose=args.verbose,
                                  permute=args.permute, model_type=args.model, fn_permute_pvals=fn_permute_pvals,
                                  fn_permute_betas=fn_permute_betas, fn_permute_stds=fn_permute_stds,
                                  indx=i+1, fn_model=fn_model, get_residual=args.get_residual, fn_residual=fn_residual)
        print(f'\r# - Processing: {i+1}/{len(lst_phenotype)}'+' '*20)
    
    # Record running time
    duration = datetime.datetime.now() - start
    duration_d = duration.days
    duration_h = duration.total_seconds()//3600
    duration_m = duration.total_seconds()//60 - duration_h
    duration_s = duration.total_seconds()%60
    msg = f'# - Regression finished in {duration_d} days, {duration_h} hours, {duration_m} minutes, {duration_s:.4f} seconds'
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


