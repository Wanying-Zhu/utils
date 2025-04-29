# By Wanying Zhu
# Run meta meta-analysis on multiple regression 
'''
Example call
1. Print out the help message for each flag:
python /data100t1/home/wanying/lab_code/utils/meta_analysis/src/meta_analysis.py --help

2. Run meta analysis
python /data100t1/home/wanying/lab_code/utils/meta_analysis/src/meta_analysis.py \
--input_files xxx xxx \
--output_path ./ \
--output_prefix output \
--pval_cols P \
--se_cols SE \
--beta_cols BETA \
--shared_cols protein_id \
--overwrite

'''


import pandas as pd
import numpy as np
import argparse
import logging
import sys
import os

from _setup_logger_and_parse_args import process_args

# ########## Get arguments from terminal ##########
args = process_args()
output_fn = os.path.join(args.output_path, args.output_prefix+'.meta_output.txt')

# ########## Load regression results ##########
# if 
