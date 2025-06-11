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
                        help='''Input file names of regression results, separated by space. 
                        Files must have column headers, columns of pvalue and beta. 
                        Provide delimiter(s) or infer from the suffix
                        ''')
    parser.add_argument('--input_delimiter', nargs='*', default=[],
                        help="(Optional) Delimiters of input file. Or provide one value if delimiters are the same in all input files. Use ',', 'tab', 'space' or other values")
    parser.add_argument('--output_path', type=str, default='./',
                        help='Output path')
    parser.add_argument('--output_prefix', type=str, default='meta_output',
                        help='Output prefix')
    
    parser.add_argument('--id_mapping_fn', nargs='?', default=None, type=str,
                        help='''(Optional) Needed when marker IDs do not match in results 1 and 2. 
                                A file with an ID mapping scheme. Expects 3 columns:
                                (1) Marker ID column in result 1, 
                                (2) Marker ID column in result 2, 
                                (3) Shared marker IDs. 
                                Column names of this file are provided by --id_mapping_cols
                             ''')
    parser.add_argument('--id_mapping_cols', nargs='+', default=[],
                        help='''(Optional) Column names in file --id_mapping_fn. 
                        The 1st, 2nd and 3rd values refer to column names corresponding to these columns in --id_mapping_fn:  
                        (1) ID column map to result 1,
                        (2) ID column map to result 2, and 
                        (3) shared ID column. 
                        If omitted, then assume the ID columns are the same as --marker_cols, 
                        and use meata_id for the shared ID column 
                        ''')
    
    parser.add_argument('--pval_cols', nargs='+',
                        help='''Column names of pvalue in each input file, separated by space.
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
    parser.add_argument('--marker_cols', nargs='+',
                        help='''Names of shared marker ID column to merge the regression results, separated by space.
                        Or provide one value if column names are the same in all input files.
                        If the ID columns do not match, then must supply a file via --id_mapping_fn.
                        This file should contain at least 3 columns:
                        (1) An ID column of result 1: with the same column header as the marker ID column in result 1; 
                        (2) An ID column of result 2: with the same column header as the marker column in result 2; 
                        ''')
    parser.add_argument('--extra_cols_to_keep', nargs='*', type=str, default=[],
                        help='Extra columns in the individual result (and ID mapping file if provided) to keep in the meta output')
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
        
    # Check if the number of column names is valid
    for val in ['pval_cols', 'beta_cols', 'marker_cols']:
        if len(eval(f'args.{val}'))!=1 and len(eval(f'args.{val}'))!=len(args.input_files):
            logging.info('# Error: %s' % '--'+val)
            logging.info('# - Number of column names must match number of input files')
            logging.info('# - Or provide one value if column names are the same in all input files\n')
            logging.info('# Exit')
            exit()

    # If an ID mapping file is provided, check if column names of this file is provided
    if args.id_mapping_fn is not None:
        if len(args.id_mapping_cols)!=0 and len(args.id_mapping_cols)!=len(args.input_files)+1:
            # The column names can be omitted, or must be the same number as imput files + 1 (extra shared ID column)
            logging.info('# Error: %s' % '--id_mapping_fn')
            logging.info('# - Number of column names must match number of input files+1')
            logging.info('# - Or omit and use the column names are the same as --marker_cols and meta_id for shared ID column\n')
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
    
