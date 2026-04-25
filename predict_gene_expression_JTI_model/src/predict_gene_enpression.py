# Author: Wanying Zhu
# Date: 2026-03-26

# Predict gene expression based on genotype VCF
# (must be tabix indexed with .tbi or .csi file in the same directory)
# Use the predixcan-like .db file to load weights

# !!![Important]!!!
# This code looks up variants via chromosome number and position to gain speed,
# So the original predixcan or JTI models (.db) will not work
# Need to use the .db file with extra columns snp_id_GRch37 and/or snp_id_GRch37 in the weight table
# The SNP ids in these columns are in format: chr:pos:ref:alt (eg. chr4:77356246:T:C or 4:77356246:T:C)
# (Your genotype VCF can still use rsid in the ID column, as long as the .db file contains the extra columns)
# What can you do:
# Solution 1: Create your own .db models (for example, update from predixcan or JTI model)
# Solution 2: Use Wanying's updated versions of JTI models here:
# /data100t1/home/wanying/BioVU/202505_hypophosphatasia/data/harmonized_predixcan_db/*.with_snp_id.db

'''
Example call
1. To run a few genes from a list
database=/data100t1/home/wanying/BioVU/202505_hypophosphatasia/data/harmonized_predixcan_db/JTI_Whole_Blood.with_snp_id.db

python /data100t1/home/wanying/lab_code/utils/predict_gene_expression_JTI_model/src/predict_gene_enpression.py \
--model_db_path ${database} \
--model_db_snp_key snp_id_GRch38 \
--vcf_genotypes /data100t1/share/BioVU/agd_250k/vcf-converted/agd250k_chr1.primary_pass.vcf.gz \
--output_path /data100t1/home/wanying/BioVU/20250724_predixan_AGD_250k/output \
--output_prefix output.selected_gene_from_list \
--only_entries ENSG00000001629 ENSG00000162551 ENSG00000001460 \
--overwrite \
--save_vcf \
--verbose

# Add --chr_in_vcf if chr is included in the VCF #CHR column
# To run one or few genes: --only_entries ENSG00000162551

2. To run a few genes provided in a file
database=/data100t1/home/wanying/BioVU/202505_hypophosphatasia/data/harmonized_predixcan_db/JTI_Whole_Blood.with_snp_id.db

python /data100t1/home/wanying/lab_code/utils/predict_gene_expression_JTI_model/src/predict_gene_enpression.py \
--model_db_path ${database} \
--model_db_snp_key snp_id_GRch38 \
--vcf_genotypes /data100t1/share/BioVU/agd_250k/vcf-converted/agd250k_chr1.primary_pass.vcf.gz \
--output_path /data100t1/home/wanying/BioVU/20250724_predixan_AGD_250k/output \
--output_prefix output.selected_genes_from_file \
--only_entries_fn /data100t1/home/wanying/BioVU/20250724_predixan_AGD_250k/code/utils/example_file/selected_gene_list.txt \
--overwrite


3. To run all genes of a given chromosome (eg. chr1)
database=/data100t1/home/wanying/BioVU/202505_hypophosphatasia/data/harmonized_predixcan_db/JTI_Whole_Blood.with_snp_id.db
ref_file=/data100t1/home/wanying/BioVU/20250724_predixan_AGD_250k/code/utils/data/gene_chr_reference.chr1.txt
python /data100t1/home/wanying/lab_code/utils/predict_gene_expression_JTI_model/src/predict_gene_enpression.py \
--model_db_path ${database} \
--model_db_snp_key snp_id_GRch38 \
--vcf_genotypes /data100t1/share/BioVU/agd_250k/vcf-converted/agd250k_chr1.primary_pass.vcf.gz \
--output_path /data100t1/home/wanying/BioVU/20250724_predixan_AGD_250k/output \
--output_prefix output.selected_genes_from_file \
--only_entries_fn <(awk 'NR>1 {print $1}' ${ref_file}) \
--overwrite
'''

import pandas as pd
import numpy as np
import logging
import argparse
import os
import sys
import time
import sqlite3
import subprocess
import gzip
import functools

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

def record_args(args):
    '''
    Record commandline call and arguments used
    '''
    # Record script used
    cmd_used = 'python ' + ' '.join(sys.argv)
    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    logging.info('# Arguments used:')
    for arg in vars(args):
        msg = f'# - {arg}: {getattr(args, arg)}'
        logging.info(msg)
    logging.info('')
        
def add_arguments():
    '''
    Take arguments from commendline
    '''
    parser = argparse.ArgumentParser(description='Predict transcriptome from a model.')
    
    parser.add_argument("--model_db_path", type=str,
                        help="Name of model db in data folder. "
                             "PrediXcan will filter input GWAS snps that are not present in the model")
    parser.add_argument("--model_db_snp_key", default='rsid', type=str,
                        help="Deafult: rsid. Specify a key to use as snp_id")
    parser.add_argument("--vcf_genotypes", type=str,
                        help="File name of the genotypes vcf to use. Must be tabix indexed with .tbi or .csi file in the same directory")
    parser.add_argument("--build", default=38, choices=[37,38],
                        help='GRCh build of the genotype VCF file, either 37 or 38')
    # Eg. tabix (system tabix does not recognize .csi index), /usr/bin/tabix, /data100t1/home/wanying/miniforge3/envs/jupyter_env/bin/tabix
    parser.add_argument("--tabix", type=str, default='/data100t1/home/wanying/miniforge3/envs/jupyter_env/bin/tabix', 
                        help="The path to tabix")
    parser.add_argument("--chr_in_vcf", action='store_true',
                        help="If true, assuming the #CHROM (1st) column of the VCF contains string 'chr' (ie. chr1 rather than 1)")
    parser.add_argument("--output_path", default='./', type=str,
                        help="Output directory. Will create the folders if not exist")
    parser.add_argument("--output_prefix", default='output', type=str,
                        help="Prefix of the file to put results in")
    parser.add_argument("--only_entries", default=[], nargs="*",
                        help="Compute only these entries in the models. Provide as genes ensembl IDs (ENSG00000139618), separate by space")
    parser.add_argument("--only_entries_fn", default=None, type=str,
                        help="Compute only these entries in the models. Provide in a file, with one genes ensembl ID per row and no header line")
    parser.add_argument("--overwrite", action='store_true',
                        help="Overwrite existing files if True")
    parser.add_argument("--save_vcf", action='store_true',
                        help="Save vcf of loaded SNPs (model_snps.loaded.vcf.gz), and a list of missing SNPs (model_snps.missing) (Do not recommand if predicting all or many genes)")                 
    parser.add_argument("--verbose", action='store_true',
                        help='Print out prediction status of each gene')
    args = parser.parse_args()
    # args, unknown = parser.parse_known_args()
    
    # Sanity checks
    if not os.path.isdir(args.output_path):
        os.makedirs(args.output_path)
    log_fn = os.path.join(args.output_path, args.output_prefix+'.log')
    setup_log(log_fn)
    logging.info('# Load arguments and setup log file')
    
    exit_flag = False
    if not args.model_db_path:
        logging.info('# Error: required argument missing: --model_db_path')
        exit_flag = True
    if not args.vcf_genotypes:
        logging.info('# Error: required argument missing: --vcf_genotypes')
        exit_flag = True

    output_fn = os.path.join(args.output_path, args.output_prefix+'.predicted_expression')
    output_summary_fn = os.path.join(args.output_path, args.output_prefix+'.summary.txt')
    if os.path.isfile(log_fn) and os.path.isfile(output_fn):
        if not args.overwrite:
            logging.info('# Error: output files exist. Remove the files or use --overwrite ')
            exit_flag
        else:
            # Remove existing file: (exclude the log file)
            # - *.model_snps.missing
            # - *.model_snps.loaded.vcf.gz (This is the on that will hold gzip, so it must be removed. Other files can be overwritten easily)
            # - *.model_summary
            for suffix in ['model_snps.missing', 'model_snps.loaded.vcf.gz']:
                subprocess.run(f"rm {os.path.join(args.output_path, args.output_prefix+'.*.{suffix}')}", shell=True)
            
    if exit_flag:
        logging.info('# Exit')
        exit()
        
    record_args(args)
    return args

def log_running_time(start_time):
    '''
    Log time elapsed
    Params:
    - start: start time (by time.time())
    '''
    end_time = time.time()
    elapsed_time = end_time - start_time
    # Print the elapsed time
    if elapsed_time>3600:
        logging.info('Elapsed time: %.2f hours' % (elapsed_time/3600))
    elif elapsed_time>60:
        logging.info('Elapsed time: %.2f minutes' % (elapsed_time/60))
    else:
        logging.info('Elapsed time: %.2f seconds' % (elapsed_time/1))

def dosage2number(val, na_rep=['.','NA'], flip=False):
    '''
    Convert genotype dosage string to number
    Assume missing values are '.' or 'NA' be default
    Params:
    - na_rep: a array of strings to label missing values
    - flip: flip dosage if True (flipped dosage = 2-original dosage)
    - val: string value to convert
    Return:
    - numeric value of dosage
    '''
    if val not in [na_rep]:
        if not flip:
            return float(val)
        else:
            return 2-float(val)
    else: # Return 0 for missing values, this is OK for predictions
        return 0

def genotype2number(val, na_rep='NA', flip=False):
    '''
    Convert genotype dosage string to number
    Change missing values as needed
    Prams:
    - na_rep: a array of strings to label missing values
    - flip: flip genotype if True (flipped gt = 2-original gt)
    - val: string value to convert
    Return:
    - numeric value of genotype
    '''
    # A dictionaty to map genotype string to numbers
    dict_gt_conversion = {'0/0':0, '0/1':1, '1/0':1, '1/1':2, './.':0,
                          '0|0':0, '0|1':1, '1|0':1, '1|1':2, '.|.':0,
                          '.':0, '..':0, na_rep:0}
    if not flip:
        return dict_gt_conversion[val]
    else:
        return 2-dict_gt_conversion[val] 

def get_genotypes(snps, snp_chr_pos, ref_alleles, alt_alleles, vcf,
                  model_snps_output_prefix, save_vcf=False, verbose=False):
    '''
    Get the genotypes for the variants needed to predict the expression of a single gene.
    Tabix is used to complete the query, so variant's chromosome and position are needed
    ! Variant's chr and position are included ('snp_id_GRch37', 'snp_id_GRch38') in Wanying's updated .db file,
    but not available in the original JTI models.
    The code will compare ref and alt alleles, and change genotype (or dosage) to 2-origianl_value if need to flip
    Params:
    - snps: an array of SNPs to look up (actual query will be based on chr and pos)
    - snp_chr_pos: array of chr:pos-pos to query vcf
    - ref_alleles, alt_alleles: list of alt or ref allele (used to match SNPs)
    - vcf: gzip and tabix-indexed vcf. Assuming the values only contains genotype such as '0/1', '0|0'
           TODO: Need additional functionalities to parse values with multiple fields from VCF
    - model_snps_output_prefix: save dosage or genotype of loaded SNPs, and list of missing SNPs for future reference
    - save_vcf: If true, save a vcf of loaded snps and a file of unfound snps
    - verbose: print out snp loading status
    Return:
    - Save loaded genotypes to a vcf file and gzip
    - n_snps_used: number of SNPs loaded
    - genotype: a dataframe with columns: snp_id, genotypes for each individual in the vcf
    '''
    # use tabix to query vcf file by position
    n_total_snps = len(snp_chr_pos)
    n_snps_found, n_snps_not_found = 0, 0 # Number of SNPs found and not found
    lst_gt = [] # Genotype or dosage of all SNPs
    lst_snps = [] # Track which SNPs are loaded

    if save_vcf:
        # Output header lines to loaded_snp_vcf
        loaded_snp_vcf = f'{model_snps_output_prefix}.model_snps.loaded.vcf'
        cmd = f"{TABIX} -H {vcf} > {loaded_snp_vcf}"
        cmd_run= subprocess.run(cmd, shell=True)
        loaded_snp_fh = open(loaded_snp_vcf, 'a') # Keep adding more lines

        # Also save SNPs that are not found
        missing_snp_fh = open(f'{model_snps_output_prefix}.model_snps.missing', 'w')

    for i, chr_pos in enumerate(snp_chr_pos):
        '''
        Something strange on the server vgipiper06 (I do not have a good explaination for it):
        If using tabix in this python code, then it points to and old version (Version: 0.2.5 (r1005)) /data100t1/gapps/tabix.
        This version does not recognize .csi index and throws our error:
         - [tabix] the index file either does not exist or is older than the vcf file. Please reindex.

        But when using tabix in the terminal, it actually points to a newever version of tabix (Version: 1.20),
        although "which tabix" gives the same pointer to /data100t1/gapps/tabix.

        When using "whereis tabix", I found multiple locations of tabix (both works for this script):
        - /usr/bin/tabix
        - /belowshare/vumcshare/data100t1/home/wanying/miniforge3/envs/jupyter_env/bin/tabix
        - /belowshare/vumcshare/data100t1/gapps/tabix (old version, does not work)
        
        Therefore, /usr/bin/tabix is specified here to avoid confusions
        '''
        if chr_pos==-1: # SNP not found: SNP does not have valid chr:pos (None value) in the .db file
            n_snps_not_found += 1
            if save_vcf:
                missing_snp_fh.write(snps[i]+'\t'+chr_pos+'\n')
            continue
                
        ref_allele, alt_allele = ref_alleles[i], alt_alleles[i]
        cmd = f"{TABIX} {vcf} {chr_pos}"
        cmd_run= subprocess.run(cmd, shell=True, text=True, capture_output=True)
        result = cmd_run.stdout
        err_msg = cmd_run.stderr # Capture any error message from terminal

        if len(err_msg)>0: # SNP not found: Log error
            n_snps_not_found += 1
            logging.info('# Error loading genotype: %s' % err_msg)
            if save_vcf:
                missing_snp_fh.write(snps[i]+'\t'+chr_pos+'\n')
        else: # SNP found: Porcess genotype into numbers
            lst_results = result.split('\n')
            if len(lst_results)==1: # No SNPs found
                n_snps_not_found += 1
                if save_vcf:
                    missing_snp_fh.write(snps[i]+'\t'+chr_pos+'\n')
                continue

            snp_loaded = False # Track if current SNP is loaded
            for line in lst_results:
                if line!='':
                    # Convert to numbers
                    # Check if ref and alt alleles match in the VCF, otherwise flip the dosage
                    tmp = line.split()
                    ref, alt = tmp[3], tmp[4]

                    # Check header lines of the VCF:
                    # ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                    # ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
                    if tmp[8]=='GT': # Genotype to numbers
                        if ref_allele==alt and alt_allele==ref: # Need to flip
                            values = [genotype2number(val, flip=True) for val in tmp[9:]]
                            n_snps_found += 1
                            snp_loaded = True
                        elif ref_allele==ref and alt_allele==alt: # Do not flip
                            values = [genotype2number(val, flip=False) for val in tmp[9:]]
                            n_snps_found += 1
                            snp_loaded = True
                        else: # Mismatch, not this SNP
                            continue
                    elif tmp[8]=='DS': # Dosages to numbers
                        if ref_allele==alt and alt_allele==ref: # Need to flip
                            values = [dosage2number(val, flip=True) for val in tmp[9:]]
                            n_snps_found += 1
                            snp_loaded = True
                        elif ref_allele==ref and alt_allele==alt: # Do not flip
                            values = [dosage2number(val, flip=False) for val in tmp[9:]]
                            n_snps_found += 1
                            snp_loaded = True
                        else: # Mismatch, not this SNP
                            continue
                    else:
                        logging.info('# ERROR: Unrecognized fileds in the values. For now the code is designed to run on VCFs with only genotypes or dosages in the values')
                        logging.info('# Exit')
                        exit()
                else:
                    # Reach end of the line, so nothing
                    continue
                
                if not snp_loaded: # If no match found
                    n_snps_not_found += 1
                    if save_vcf:
                        missing_snp_fh.write(snps[i]+'\t'+chr_pos+'\n')
                else:
                    lst_snps.append(snps[i])
                    lst_gt.append(values)
                    if save_vcf:
                     loaded_snp_fh.write(line+'\n')
            if verbose:
                print(f'\r# Load genotypes for variants in the model: {n_snps_found}/{n_total_snps}; {n_snps_not_found} SNPs not found',
                      flush=True, end='')
    if verbose:
        print()
        logging.info('# Load genotypes for variants in the model: %s/%s; %s SNPs not found' % (n_snps_found, n_total_snps, n_snps_not_found))
    if save_vcf:
        loaded_snp_fh.close()
        missing_snp_fh.close()

    # Get column headers to label loaded dataframe of genotypes
    with gzip.open(vcf, 'rt') as fh:
        for line in fh:
            if line[:2] != '##':
                break
    headers = line.strip().split()[9:]
    # Compress loaded SNPs file, load into dataframe and return
    if save_vcf:
        subprocess.run(f'gzip {loaded_snp_vcf}', shell=True)
    # df = pd.read_csv(f'{loaded_snp_vcf}.gz', compression='gzip', sep='\t', comment='#', header=None)
    # df.columns = headers
    # # Columns are #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG001_PrecFDA_H
    # return df.drop(columns=['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    df_gt = pd.DataFrame(lst_gt, columns=headers, index=lst_snps)
    return n_snps_found, df_gt
    
def get_weights(gene, conn):
    '''
    Get weights for the variants needed to predict the expression of a single gene
    Params:
    - gene: gene of interest
    - conn: connection to the database
    Return:
    - weights: a dataframe with at least these columns: snp_id (or user defined snp_id column), ref_allele, eff_allele, weights
      (All columns in the updated db are: 'rsid', 'gene', 'weight', 'ref_allele', 'eff_allele', 'snp_id_GRch37', 'snp_id_GRch38')
    '''
    # snps = ['chr4:76435093:T:C', 'chr4:163447245:T:C', 'chr4:163447245:T:C', 'chr4:186837848:C:T']
    # cmd = f"SELECT * FROM weights WHERE snp_id_GRch38 IN ({','.join('?'*len(snps))})"
    cmd = f"SELECT * FROM weights WHERE gene='{gene}'"
    df_weights = pd.read_sql_query(sql=cmd, con=conn)
    return df_weights
    
def get_snp_chr_pos(snp_id_col, chr_in_vcf, df_weights):
    '''
    Get the chromosomes and positions for SNPs to query VCF file.
    Assuming there are columns in the weight table with SNP id in chr:pos:ref:alt format
    Params:
    - snp_id_col: column name for the SNP ids (valid values are: 'rsid', 'snp_id_GRch37', 'snp_id_GRch38')
    - chr_in_vcf: If true, assuming 'chr' is included in the chromosome column of the VCF
    - df_weights: a dataframe of snp weights from the .db model for a given gene
    Return:
    - An array of SNP info as chr:position-position to query VCF
    '''
    # Column header of SNP_id are limited to these 3 names
    if snp_id_col not in ['rsid', 'snp_id_GRch37', 'snp_id_GRch38']:
        logging.info('# ERROR: This code expects SNP id column to be rsid, snp_id_GRch37 or snp_id_GRch38')
        logging.info('# Modify your model if using customized .db file with different name for the SNP id column')
        logging.info('# Exit')
    if snp_id_col=='rsid':
        # Load SNP chromosome and position to query vcf
        try:
            if build==37:
                snp_ids = df_weights['snp_id_GRch37'].values
            elif build==38:
                snp_ids = df_weights['snp_id_GRch38'].values
            else:
                logging.info('# Error: Wrong GRCh build. Valid values are 37 or 38')
                logging.info('# Exit')
                exit()
        except:
            logging.info('# Error: this code look up variants via chromosome number and position.')
            logging.info('# However not columns named snp_id_GRch37 or snp_id_GRch38 was found in the weight table')
            logging.info('# Hint: Did you use the original JTI or predixcan .db file which only contains rsids?')
            logging.info('# If yes, ask Wanying for updated JTI .db file with snp_id_GRch37 and snp_id_GRch38.')
            logging.info('# Exit')
            exit()
    else:
        snp_ids = df_weights[snp_id_col].values

    # Some rsids do not have valid chr:pos in the new database, they need special care (label as -1)
    snp_ids = [snp if snp is not None else -1 for snp in snp_ids]
    
    # Check if string 'str' is in the SNP id and vcf
    if len(snp_ids)==0: # If no SNPs in the model (should not happen)
        logging.info('# Warning: no SNPs in the model from .db file (which should not happen)')
    elif 'chr' in snp_ids[0]:
        # logging.info("# String 'chr' is detected in the chromosome number. Make sure the format matches in the VCF file")
        if not chr_in_vcf:
            # Remove 'chr' if the vcf file does not have chr in the chromosome column
            # Some rsids do not have valid chr:pos in the new database, they need special care (label as -1)
            snp_ids = [snp.split('chr')[-1] if snp!=-1 else -1 for snp in snp_ids]
    else:
        # Add chr to snp list if needed (May not be necessary for GRCh37 VCFs)
        logging.info("# String 'chr' is not in the chromosome number. Make sure the format matches in the VCF file")
        if chr_in_vcf:
            snp_ids = ['chr'+snp for snp in snp_ids]

    # Return chr:pos-pos to query vcf
    # Original ID format is chr1:20633561:C:T
    return [f"{snp.split(':')[0]}:{snp.split(':')[1]}-{snp.split(':')[1]}" if snp!=-1 else -1 for snp in snp_ids]
        
            
def predict_a_single_gene(gene, snp_id_col, chr_in_vcf, vcf,
                          model_snps_output_prefix, build, conn,
                          save_vcf=False, verbose=False):
    '''
    Get the predicted expression of a single gene
    Params:
    - gene: gene id as stored in the .db file
    - snp_id_col: column name of the snp id column
    - chr_in_vcf: (True/False) If ture, assume string chr is included in the chromosome column of the VCF
    - vcf: tabixed vcf file to query genotype
    - model_snps_output_prefix: Prefix of file names. Save load SNPs and missing SNPs for reference
    - build: GRCh build (37 or 38)
    - conn: connection to the database
    - save_vcf: If true, save tabix quaried results of loaded snps into a vcf and gzip.
                Also keep a list of missing snps.
                DO NOT turn this flag on if predicting many genes, otherwise may take up a lot of storage space.
    - verbose: print out status of each gene
    Return:
    - n_snps_used: number of SNPs loaded
    - expression: predicted expression of the given gene
    '''
    # Load weights
    df_weights = get_weights(gene=gene, conn=conn)
    if len(df_weights)==0:
        # Gene not found
        return 0, None

    # Need these columns from the weight dataframe
    cols = [snp_id_col, 'ref_allele', 'eff_allele', 'weight']
    snps = df_weights[snp_id_col].values # Get genotype of these SNPs from VCF
    lst_ref = df_weights['ref_allele'].values
    lst_alt = df_weights['eff_allele'].values
    # Get an array of SNP chromosome and position (chr:pos-pos) to query VCF for genotype
    snp_chr_pos = get_snp_chr_pos(snp_id_col=snp_id_col, chr_in_vcf=chr_in_vcf, df_weights=df_weights)

    n_snps_used, genotypes = get_genotypes(snps=snps, snp_chr_pos=snp_chr_pos,
                                           ref_alleles=lst_ref, alt_alleles=lst_alt,
                                           vcf=vcf,
                                           model_snps_output_prefix=f'{model_snps_output_prefix}.{gene}',
                                           save_vcf=save_vcf,
                                           verbose=verbose)
   
    if len(genotypes)==0: # If no SNPs found at all
        return 0, None
    weights = df_weights.set_index(keys=snp_id_col).reindex(genotypes.index)['weight'].values
    expressions = np.matmul(weights, genotypes.values)
    sample_ids = genotypes.columns
    return n_snps_used, pd.DataFrame({'sample_id':sample_ids, gene:expressions})
    
if __name__=='__main__':
    # Load arguments from the terminal
    args = add_arguments()
    TABIX = args.tabix # Set the global tabix (system tabix does not recognize .csi index file)
    outout_fn = os.path.join(args.output_path, args.output_prefix+'.predicted_expression')
    model_summary_fn = os.path.join(args.output_path, args.output_prefix+'.model_summary') # Save model info (same as predixcan)
    fh_model_summary = open(model_summary_fn, 'w')
    fh_model_summary.write('gene\tgene_name\tn_snps_in_model\tn_snps_used\tpred_perf_r2\tpred_perf_pval\n')
    # Record the start time
    start_time = time.time()
    
    # Connect to database
    conn = sqlite3.connect(args.model_db_path)
    # Load genes (including names and model details) from table 'extra'
    cur = conn.cursor()
    cmd = f"SELECT * FROM extra"
    df_genes = pd.read_sql_query(sql=cmd, con=conn)
    # print(df_genes.head(10))
    # Load selected genes from file or list provided
    if args.only_entries_fn is not None:
        lst_selected_genes = []
        with open(args.only_entries_fn) as fh:
            for line in fh:
                if line.strip() != '':
                    lst_selected_genes.append(line.strip())
        df_genes = df_genes[df_genes['gene'].isin(lst_selected_genes)]
    elif args.only_entries != []: # Only predict selected genes
        df_genes = df_genes[df_genes['gene'].isin(args.only_entries)]
    
    # Skip rest if no gene left for prediction
    if len(df_genes)==0:
        logging.info('# No gene left for prediction')
        logging.info('# Total running time:')
        log_running_time(start_time=start_time)
        logging.info('# Exit')
        exit()
    
    df_all_expressions = ''
    total_n_genes, count = len(df_genes), 0
    for i, (_, row) in enumerate(df_genes.iterrows()):
        cur_start_fime = time.time()
        
        # To match the output of predixcan
        gene = row['gene']
        gene_name = row['genename']
        n_snps_in_model = row['n.snps.in.model']
        pred_perf_r2 = row['pred.perf.R2']
        pred_perf_pval = row['pred.perf.pval']
        if args.verbose:
            logging.info(f'\n# Predict {gene} ({gene_name})')
        # Calculate predicted expression
        n_snps_used, df_expressions = predict_a_single_gene(gene=gene,
                                                            snp_id_col=args.model_db_snp_key,
                                                            chr_in_vcf=args.chr_in_vcf,
                                                            vcf=args.vcf_genotypes,
                                                            model_snps_output_prefix = os.path.join(args.output_path, args.output_prefix),
                                                            build=args.build,
                                                            conn=conn,
                                                            save_vcf=args.save_vcf,
                                                            verbose=args.verbose)

        if df_expressions is None:
            # Move to next if no SNPs for this gene is found
            print(f'\r# Processed genes {i+1}/{total_n_genes}; {count} successful predictions', flush=True, end='')
            continue

        fh_model_summary.write(f'{gene}\t{gene_name}\t{n_snps_in_model}\t{n_snps_used}\t{pred_perf_r2}\t{pred_perf_pval}\n')
        fh_model_summary.flush()
        # Combine with other genes
        if len(df_all_expressions)==0:
            df_all_expressions = df_expressions.copy()
        else:
            df_all_expressions[gene] = df_expressions[gene]
        count += 1
        if args.verbose:
            log_running_time(cur_start_fime)
        else: # Turn it off when running verbose mode
            print(f'\r# Processed genes {i+1}/{total_n_genes}; {count} successful predictions', flush=True, end='')
    print()
    logging.info('# Processed genes %s/%s; %s successful predictions' % (i+1, total_n_genes, count))
    
    # Write to output file
    logging.info('# Save predictions to output file')
    if len(df_all_expressions) != 0:
        df_all_expressions.to_csv(outout_fn, sep='\t', index=False)
    fh_model_summary.close()
    cur.close()
    
    # Calculate the time elapsed
    logging.info('# Total running time:')
    log_running_time(start_time=start_time)




