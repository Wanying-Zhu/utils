'''
(Similar to 02_check_flip_dosage.py)
Takes a list of SNPs with ref, alt, chr, pos, effect allele and prs weights,
Query VCF file by chr and position.
Flip dosage as 2-old_dosage if effect allele is not alt allele
Save processed flipped_dosage*weight in output file

Example call:
python /data100t1/home/wanying/BioVU/202405_PAGE_PRS/code/utils/get_weighted_dosage.py \
--lipid HDL \
--method PT \
--chromosome 1 \
--ancestry EUR \
--vcf example_data/AFR.chr1.vcf.gz \
--snp_list example_data/chr1_snp.list \
--output_fn output/chr1.weighted_dosage \
--impute \
--mean_af_col_anc MeanAF_African_American \
--mean_af_col MeanAF \
--overwrite
'''

import pandas as pd
import numpy as np
import sys
import subprocess
import os
import gzip
import argparse
import logging

# #################### Helper funcitons ####################

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
    parser.add_argument('--lipid', choices=['HDL', 'nonHDL', 'LDL', 'logTG', 'TC'])
    parser.add_argument('--method', choices=['PRSCS', 'PT'])
    parser.add_argument('--chromosome')
    parser.add_argument('--ancestry', default='EUR', type=str, choices=['AFR', 'EAS', 'EUR', 'HIS', 'SAS'])
    parser.add_argument('--vcf', type=str, default='/data100t1/share/BioVU/TOPMed_imputed',
                        help='VCF file to look up SNPs')
    parser.add_argument('--snp_list', type=str,
                        default='/data100t1/home/wanying/BioVU/202405_PAGE_PRS/supporting_files/redo/HDL/PT/chr22_snp.list',
                        help='File name of SNP list with pos, alt, ref, effect allele and PRS weight')
    parser.add_argument('--output_fn', type=str)
    parser.add_argument('--impute', action='store_true',
                        help='Impute missing SNPs using mean allele count (= mean allele frequency * 2)')
    parser.add_argument('--mean_af_col_anc', type=str, default=None,
                        help='Column name of ancestry specific mean allele frequency (for imputation). Such as MeanAF_European.')
    parser.add_argument('--mean_af_col', type=str, default='MeanAF',
                        help='Column name of mean allele frequency (for imputation). Use this one if ancestry-specific AF is missing')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing BGzipped file')
    args = parser.parse_args()
    output_path = os.path.split(args.output_fn)[0]
    if not os.path.isdir(output_path): # Create output folder if not exists
        os.makedirs(output_path)
    return args

def process_dosage_of_single_snp(line, ref, alt, effect_allele):
    '''
    Read in ref, alt and effect allele in the snp list, 
    Extract dosages from a single line of vcf and flip according the information.
    Params:
    - line: a line from VFC query
    - ref, alt, effect_allele: ref, alt, effect allele info from the SNP list
    '''
    
    # Get dosage
    lst_dosage = [] # Need to handle .|. values. Set np.nan for missing values
    
    # Get ref and alt alleles from vcf to compare with SNP list
    _,_,_, ref_vcf, alt_vcf ,_,_,_,_, dosages = line.split(maxsplit=9)
    for val in dosages.split():
        try: lst_dosage.append(val.split(':')[1])
        except: lst_dosage.append(np.nan)
    
    # In case there are values as .|.:.:.:. that slip through the missing value check above
    for i in range(len(lst_dosage)):
        try:
            lst_dosage[i] = float(lst_dosage[i])
        except:
            # Store none for missing dosage
            # lst_dosage[i] = np.nan

            # If any missing value is detected, mark this SNP as missing,
            # so that all samples use the same set of SNPs in PRS calculation
            return None
    # Convert to numpy array for fast calculation
    array_dosage = np.array(lst_dosage)

    flip = False
    if ref==ref_vcf and alt==alt_vcf:
        if effect_allele==alt_vcf: # No flip needed
            pass
        else:
            flip=True
    elif alt==ref_vcf and ref==alt_vcf:
        if effect_allele==alt_vcf: # No flip needed
            pass
        else:
            flip=True
    else:
        # A SNP at the given position is found, but not the one in the list
        # Label the SNP as not found
        return None
        
    if flip: array_dosage = 2-array_dosage
    return array_dosage

def process_multiallelic(lines, pos, ref, alt):
    '''
    Take in tabix search result of multiallelic site in a list, return the true one based on SNP id searched
    Params:
    - lines: result from the tabix in a list. One entry per line
    - pos, ref, alt: pos, ref, alt of the SNP to be looked for. Such as 14903694, G, A
    Return:
    The line of snp looking for. Or empty string if no match is found
    '''
    for line in lines:
        tmp_lst = line.split(maxsplit=9)
        cur_pos, cur_ref, cur_alt = tmp_lst[1], tmp_lst[3], tmp_lst[4]
        if cur_pos==str(pos):
            if (cur_ref==ref and cur_alt==alt) or (cur_ref==alt and cur_alt==ref):
                return line
    return ''

def calcualte_weighted_dosage(vcf_fn, df_snp_list, output_fn, impute=False,
                              mean_af_col_anc=None, mean_af_col=None):
    '''
    Takes in a snp list with flag indicates if the dosage needs to be flipped
    Save weight*dosage to an output file
    Params:
    - dosage_fn: Dosage vcf to look up snps
    - df_snp_list: a dataframe of snps to look for, with snp ids, positions, ref, alt, effect allele, weights
    - output_fn: output file to save processes dosage*weight
    - impute: impute missing SNPs using mean allele frequency if true
    - mean_af_col_anc: mean allele frequency in current ancestral group samples (for imputation)
    - mean_af_col: mean allele frequency in all samples (for imputation)
    '''
    not_found_fn = output_fn + '.not_found' # Save any snps that are missing from the dosage vcf
    out_fh = open(output_fn, 'w')
    not_found_fh = open(not_found_fn, 'w')
    not_found_fh.write('ID\timputed_weighted_dosage\n')
    
    # Get header of dosage file
    with gzip.open(vcf_fn, 'rt') as fh:
        line = fh.readline().strip()
        while line[:2] == '##':
            line = fh.readline().strip()
        header = line.split(maxsplit=9)[-1]
    out_fh.write('CHROM\tPOS\tID\tREF\tALT\t' + header + '\n')
    
    c_total, count_multiallelic, count_not_found = 0, 0, 0
    for i in range(len(df_snp_list)):
        c_total += 1
        snp = df_snp_list['GRCh38'].iloc[i]
        chr_num = df_snp_list['chr'].iloc[i]
        pos = df_snp_list['pos'].iloc[i]
        ref = df_snp_list['ref'].iloc[i]
        alt = df_snp_list['alt'].iloc[i]
        effect_allele = df_snp_list['PRS_effect_allele'].iloc[i]
        weight = df_snp_list['PRS_weight'].iloc[i]

        if impute:
            mean_af = df_snp_list[mean_af_col_anc].iloc[i] # Mean allele frequency of current ancestral group
            if np.isnan(mean_af): # If ancestry specific allele frequency is missing, use the mean AF instead
                mean_af = df_snp_list[mean_af_col].iloc[i] # Mean allele frequency in all samples
            w_dosage = mean_af*2 * weight # use mean_af*2 as imputed dosage
        
        cmd = f'tabix {vcf_fn} {chr_num}:{pos}-{pos}'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.rstrip()
        if result == '':
            count_not_found += 1
            if impute:
                not_found_fh.write(f'{snp}\t{w_dosage}\n')
            else:
                not_found_fh.write(f'{snp}\tNA\n')
            continue # If not found

        lines = result.split('\n')
        
        # Sanity checks
        if len(lines)>1: # If multiallelic site
            count_multiallelic += 1
            line = process_multiallelic(lines, pos, ref, alt)
            if line != '':
                # Get dosage
                lst_dosage = process_dosage_of_single_snp(line=line, ref=ref, alt=alt,
                                                          effect_allele=effect_allele)
            else: # No match
                count_not_found += 1
                # not_found_fh.write(snp+'\n')

                if impute:
                    not_found_fh.write(f'{snp}\t{w_dosage}\n')
                else:
                    not_found_fh.write(f'{snp}\tNA\n')
                continue
        else:
            # Get dosage
            lst_dosage = process_dosage_of_single_snp(line=lines[0], ref=ref, alt=alt,
                                                      effect_allele=effect_allele)
        if lst_dosage is not None:
            out_fh.write(f'{chr_num}\t{pos}\t{snp}\t{ref}\t{alt}\t' + '\t'.join([str(val) for val in weight*lst_dosage]) + '\n')
        else:
            # If lst_dosage is none, exclude the SNP from PRS and label it as missing
            count_not_found += 1
            if impute:
                not_found_fh.write(f'{snp}\t{w_dosage}\n')
            else:
                not_found_fh.write(f'{snp}\tNA\n')
            
        if c_total>10000 and c_total%100==0:
            print(f'\r# Processed: {c_total}/{len(df_snp_list)}    ', end='', flush=True)
        else:
            print(f'\r# Processed: {c_total}/{len(df_snp_list)}    ', end='', flush=True)
    msg = f'\r# Processed: {c_total}/{len(df_snp_list)}; N_multiallelic: {count_multiallelic}; N_not_found: {count_not_found} '
    logging.info(msg)
    out_fh.close()
    not_found_fh.close()

# #################### End of helper functions ####################    
args = process_args()

fn_log = '.'.join(args.output_fn.split('.')[:-1]) + '.log'
setup_log(fn_log, mode='w')

# Record script used
cmd_used = 'python ' + os.path.basename(__file__)
logging.info('# Arguments used:')
for arg in vars(args):
    cmd_used += f' --{arg} {getattr(args, arg)}'
    msg = f'# - {arg}: {getattr(args, arg)}'
    logging.info(msg)

logging.info('\n# Call used:')
logging.info('python ' + ' '.join(sys.argv))

if os.path.isfile(args.output_fn+'.gz'):
    if not args.overwrite:
        logging.info('# BGzipped output file already exists. Skip')
        exit()
    else:
        # Remove existing BGzipped file and tabix file
        subprocess.run(f'rm {args.output_fn}.gz', shell=True)
        subprocess.run(f'rm {args.output_fn}.gz.tbi', shell=True)

# ##################### Load files ####################
df_snp_list = pd.read_csv(args.snp_list, sep='\t', dtype={'pos':int})

if args.impute: # Ignore SNPs without mean allele frequencies (cannot be imputed)
    # Save SNPs that cannot be imputed
    fn_cannot_impute = args.output_fn + '.cannot_impute'
    df_snp_list[df_snp_list['MeanAF'].isna()].to_csv(fn_cannot_impute, index=False)
    df_snp_list = df_snp_list[df_snp_list['MeanAF'].notna()].copy()

calcualte_weighted_dosage(args.vcf,
                          df_snp_list.sort_values(by='pos'),
                          args.output_fn, impute=args.impute,
                          mean_af_col_anc=args.mean_af_col_anc,
                          mean_af_col=args.mean_af_col)

logging.info('# bgzip flipped dosage')
subprocess.run(f'bgzip {args.output_fn}', shell=True)

logging.info('# tabix flipped dosage')
subprocess.run(f'tabix -b 2 -e 2 {args.output_fn}.gz', shell=True)

logging.info('\n# DONE\n')
