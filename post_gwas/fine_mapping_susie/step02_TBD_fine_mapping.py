# Necessary files for susie fine-mapping
# 1. GWAS summary stats, with BETA, SE, number of samples
# 2. Phenotype file (or variance of phenotype)
# 3. Genotype file or LD matrix

import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils')
from setup_logger_and_parse_args import process_args

# Prepare files for finemapping

# Calcualte LD matrix
'''
# Pairwise LD measures for multiple SNPs (genome-wide)
# https://zzz.bwh.harvard.edu/plink/ld.shtml
plink1.9 --file mydata --r square --keep-allele-oder --out
'''