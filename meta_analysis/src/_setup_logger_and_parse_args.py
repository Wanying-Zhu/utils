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
    parser.add_argument('--input_files', nargs='+',
                        help='Input file names of regression results, separated by space. Files must have column headers, columns of pvalue and beta. Provide delimiter(s) or infer from the suffix')
    parser.add_argument('--input_delimiter', nargs='*', default=[],
                        help="Delimiters of input file. Or provide one value if delimiters are the same in all input files. Use ',', 'tab', 'space' or other values")
    parser.add_argument('--output_path', type=str, default='./',
                        help='Output path')
    parser.add_argument('--output_prefix', type=str, default='meta_output',
                        help='Output prefix')

    parser.add_argument('--pval_cols', nargs='+',
                        help='''Column names of pvalue in each input file, separated by space.
                        Or provide one value if column names are the same in all input files
                        ''')
    parser.add_argument('--se_cols', nargs='+',
                        help='''Column names of standard error in each input file, separated by space.
                        Or provide one value if column names are the same in all input files
                        ''')
    parser.add_argument('--beta_cols', nargs='+',
                        help='''Column names of beta in each input file, separated by space.
                        Or provide one value if column names are the same in all input files
                        ''')
    parser.add_argument('--sample_sizes', nargs='+',
                        help='''Sample sizes of each regression (integer value).
                        Or provide one value if column names are the same in all input files
                        ''')
    
    parser.add_argument('--shared_cols', nargs='+',
                        help='''Names of shared ID column to merge the regression results, separated by space.
                        Or provide one value if column names are the same in all input files
                        ''')
    
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output file if True. Default value is False')
    args = parser.parse_args()
    # Convert delimiter(s) to valid string(s)
    dict_delimiter = {'tab': '\t', 'space': ' '}
    if len(args.input_delimiter) != 0:
        for i in range(len(args.input_delimiter)):
            args.input_delimiter[i] = dict_delimiter[args.input_delimiter[i]]
            
    # ########## Sanity checks ##########
    if not os.path.isdir(args.output_path): # Create output folder if not exists
        print('# Create output path: ' + args.output_path)
        os.makedirs(args.output_path)
    
    # Record arguments used
    fn_log = os.path.join(args.output_path, args.output_prefix+'.log')
    setup_log(fn_log, mode='w')
    # logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message
    
    # Check if the output file exists when overwrite=False
    output_fn = os.path.join(args.output_path, args.output_prefix+'.meta_output.txt')
    if os.path.isfile(output_fn) and not args.overwrite:
        logging.info('# Output file exists: %s' % output_fn)
        logging.info('# Change output prefix or use --overwrite. Skip saving and exit')
        exit()

    # Number of input files must >=2
    if len(args.input_files)<2:
        logging.info('# Error: Must have at least 2 input files\n')
        logging.info('# Exit')
        exit()
        
    # Check if number of column names is valid
    for val in ['pval_cols', 'se_cols', 'beta_cols', 'shared_cols']:
        if len(eval(f'args.{val}'))!=1 and len(eval(f'args.{val}'))!=len(args.input_files):
            logging.info('# Error: %s' % '--'+val)
            logging.info('# - Number of column names must match number of input files')
            logging.info('# - Or provide one value if column names are the same in all input files\n')
            logging.info('# Exit')
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
    
