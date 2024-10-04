import logging
import argparse
import os
import sys

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
    Process arguments and set up a log file
    '''
    parser = argparse.ArgumentParser()
    
    # Add some default arguments
    parser.add_argument('--input_file', type=str,
                        help='Input file name, must in row x col = sample x feature format (eg. rows are samples, columns are gene expressions)')
    parser.add_argument('--output_path', type=str, default='./')
    parser.add_argument('--output_prefix', type=str, default='')
    parser.add_argument('--covar_file', type=str, default=None,
                        help='(Optional) Name of the file containing covariates if not in the input_file, in row x col = sample x feature format. Can be a .csv file or tab-delimited file. (Or a list of samples to use in the input)')
    parser.add_argument('--id_col', type=str, default='LABID', help='Shared ID column between the input file and covariate file')
    parser.add_argument('--ignore_cols', type=str, nargs='*', default=['RRID', 'LABID'],
                        help='Columns to be ignored from the analysis. Used if no single phenotype is specified')
    parser.add_argument('--covars', type=str, nargs='*', help='Covariates of the model, such as sex, age, PCs')
    parser.add_argument('--phenotype', type=str, default=None,
                        help='Outcome of the model, such as lipid concentration. Run all columns in the input file if None (except covariates and condition)')
    parser.add_argument('--scale', action='store_false', help='(True) Scale (unit variation) the data before PCA. Centering is implemented by PCA already')
    parser.add_argument('--verbose', action='store_true', help='Print more information. Default value is false')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output file if true. Default value is false')
    parser.add_argument('--threads', default=1, type=int, help='Number of threads for multiprocessing. Default is none (1)')
    
    args = parser.parse_args()

    # ########## Sanity checks ##########
    if not os.path.isdir(args.output_path): # Create output folder if not exists
        print('# Create output path: ' + args.output_path)
        os.makedirs(args.output_path)
    if args.output_prefix == '':
        args.output_prefix = os.path.split(args.input_file)[-1]
        if '.' in args.output_prefix:
            args.output_prefix = '.'.join(args.output_prefix.split('.')[:-1])
    # Record arguments used
    fn_log = os.path.join(args.output_path, args.output_prefix+'.log')
    setup_log(fn_log, mode='w')
    logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message
    
    if args.covar_file: # If covar_file is provided, check if covariate file exists
        if not os.path.isfile(args.covar_file):
            msg = '# ERROR: Covarriate file not found: ' + args.covar_file
            logging.info(msg)
            exit()
    # Check if the output file exists when overwrite=False
    if os.path.isfile(f'{args.output_path}/{args.output_prefix}.residual') and not args.overwrite:
        logging.info('# Output file exists: %s' % args.output_path+'/'+args.output_prefix)
        logging.info('# Skip saving and exit')
        exit()
        
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
    
