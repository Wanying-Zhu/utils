# Process arguments and sanity checks
import os
import time
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-o', '--output_prefix', type=str,
                               help='Prefix of output file ')
    parser.add_argument('--output_dir', type=str, help='Output directory. Default is current directory', default='./')
    parser.add_argument('--input', type=str, help='Input file name')
    parser.add_argument('--overwrite', action='store_true',
                        help='On/off flag. If true (on), will rerun and overwrite existing files')
    args = parser.parse_args()

    # Sanity checks
    sanity_checks(args)

    # logging.info('# Arguments used:')
    # for arg in vars(args):
    #     logging.info('# - %s: %s' % (arg, getattr(args, arg)))
    
    return args

def sanity_checks(args):
    # # Make all directory to absolute path
    # args.output_dir = os.path.expanduser(args.output_dir)

    # Check if files exist
    if not os.path.isfile(args.input):
        print('# ERROR: file not found:', args.input, end='')
        exit()
    else:
        print(' PASS')

    print('# Check output directory: ', end='')
    if not os.path.isdir(args.output_dir):
        print('\n#\t - Output directory does not exist. Create one at:', args.output_dir)
        os.makedirs(args.output_dir)
    else:
        print(' PASS')

    # if os.path.isfile(os.path.join(args.output_dir, args.output_prefix)) and not args.overwrite:
    #     logging.info('# Output file exists. Set --overwrite True if really want to rerun')
    #     exit()