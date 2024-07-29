import logging
import argparse

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
    # If only need to write to log file, use below
    # logging.basicConfig(level=logging.DEBUG, filename='app.log', filemode='a', format='%(name)s - %(levelname)s - %(message)s')

    # logging.debug('This is a debug message')
    # logging.info('This is an info message')
    # logging.warning('This is a warning message')
    # logging.error('This is an error message')
    # logging.critical('This is a critical message')
    '''
    # Once arguments are parsed, record script used
    cmd_used = 'python ' + os.path.basename(__file__)
    logging.info('# Arguments used:')

    for arg in vars(args):
        cmd_used += f' --{arg} {getattr(args, arg)}'
        msg = f'# - {arg}: {getattr(args, arg)}'
        logging.info(msg)
    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    '''

def process_args(log_args, *args):
    '''
    Process arguments
    Example call
    process_args(True,
                 {'flag_name':'--flag1', 'default': 'flag1 default value', 'type': 'str', 'help': 'flag1 help message'}
                 {'flag_name':'--flag2', 'default': 'flag2 default value', 'type': 'int', 'help': 'flag2 help message'} )
    Params:
    - log_args: If true, save arguments into log file
    - *args: any other arguments need to be parsed, keep details in a dictionary.
             Must have a key of 'flag_name'. Other key/value pairs are passed to parser.add_argument().
             For example: {'flag_name':'--flag1', 'default': 'flag1 default value', 'type': 'str', 'help': 'flag1 help message'}
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fn', help='Input file name', type=str)
    parser.add_argument('--output_path', type=str, default='./')
    parser.add_argument('--output_prefix', type=str, default='output')

    # Parse other arguments
    for a in args:
        # create a string of arguments to pass to parser.add_argument()
        args_to_add = f"'{a['flag_name']}'"
        for k, v in a.items():
            if k != 'flag_name':
                if k == 'type': # type is not a string
                    args_to_add += f",{k}={v}"
                else:
                    args_to_add += f",{k}='{v}'"
            
        cmd = f"parser.add_argument({args_to_add})"
        eval(cmd)
    
    terminal_args = parser.parse_args()
    
    if not os.path.isdir(terminal_args.output_path): # Create output folder if not exists
        os.makedirs(terminal_args.output_path)

    # Record arguments used
    if log_args:
        fn_log = os.path.join(terminal_args.output_path, terminal_args.output_prefix+'.log')
        setup_log(fn_log, mode='w')
        
    # Record script used
    cmd_used = 'python ' + os.path.basename(__file__)
    logging.info('# Arguments used:')
    for arg in vars(terminal_args):
        cmd_used += f' --{arg} {getattr(terminal_args, arg)}'
        msg = f'# - {arg}: {getattr(terminal_args, arg)}'
        logging.info(msg)
    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    return terminal_args
    