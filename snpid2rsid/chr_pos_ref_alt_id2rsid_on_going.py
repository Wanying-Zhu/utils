# Find rsID based on chr:position:ref:alt ID
# Download VCF reference file fromhttps://ftp.ncbi.nih.gov/snp/latest_release/VCF/
'''
Example call
python /data100t1/home/wanying/lab_code/utils/snpid2rsid/chr_pos_ref_alt_id2rsid.py \
--input_fn example_data/snp_list.txt \
--delimiter tab \
--snp_col_indx 0 \
--snp_id_delimiter : \
--input_header \
--reference /data100t1/home/wanying/shared_data_files/ncbi/GRCh38_vcf/GCF_000001405.40.gz \
--build 38 \
--output_path example_output \
--output_prefix output \
--inplace

python /data100t1/home/wanying/lab_code/utils/snpid2rsid/chr_pos_ref_alt_id2rsid.py \
--input_fn example_data/snp_list_v2.txt \
--delimiter tab \
--snp_col_indx 0 \
--snp_id_delimiter _ \
--reference /data100t1/home/wanying/shared_data_files/ncbi/GRCh38_vcf/GCF_000001405.40.gz \
--build 38 \
--output_path example_output \
--output_prefix output \
--inplace
'''

import logging
import argparse
import os
import sys
import subprocess

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
    Process and return arguments, set up logger file
    '''
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_fn', type=str,
                        help='''Input file name with SNP chr:pos:ref:alt IDs in a column.
                        Other columns are not used.
                        Examples of SNP IDs: 1:1000:A:G, chr1:1000:A:G, CHR1:1000:A:G.
                        Format of chromosome must match that in the reference file''')
    parser.add_argument('--delimiter', type=str, default='tab',
                        choices=['tab', ',', 'space'],
                        help='Delimiter of the input file')
    parser.add_argument('--snp_col_indx', type=int, help='Column index (zero-based) of the SNP id column')
    parser.add_argument('--snp_id_delimiter', type=str, default=':',
                        help='Delimiter of the SNP id. Eg. snp_id_delimiter of chr1:123:A:T is :')
    parser.add_argument('--input_header', action='store_true', help='True = the input file has a header line')
    parser.add_argument('--reference', type=str,
                        default='/data100t1/home/wanying/shared_data_files/ncbi/GRCh38_vcf/GCF_000001405.40.gz',
                        help='Bgzipped VCF reference file. Assuming the 3rd column (ID column) contains rsIDs. Must contain the corresponding tabix index file in the same directory')
    parser.add_argument('--build', type=int, choices=[37, 38], default='38',
                        help='GRCh assembly number, 37 (GRCh37) or 38(GRch38)')
    parser.add_argument('--output_path', type=str, default='./')
    parser.add_argument('--output_prefix', type=str, default='output')
    parser.add_argument('--inplace', action='store_true',
                        help='If Ture, add an SNP ID column ("rsID") to the input file with new rsIDs. Else only output old SNP IDs and rsIDs into a new file')
    
    terminal_args = parser.parse_args()
    
    if not os.path.isdir(terminal_args.output_path): # Create output folder if not exists
        print('# Create output path: ' + terminal_args.output_path)
        os.makedirs(terminal_args.output_path)

    # Record arguments used
    fn_log = os.path.join(terminal_args.output_path, terminal_args.output_prefix+'.log')
    setup_log(fn_log, mode='w')
    # logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message
        
    # Record script used
    cmd_used = 'python ' + ' '.join(sys.argv)

    logging.info('\n# Call used:')
    logging.info(cmd_used+'\n')
    
    logging.info('# Arguments used:')
    for arg in vars(terminal_args):
        cmd_used += f' --{arg} {getattr(terminal_args, arg)}'
        msg = f'# - {arg}: {getattr(terminal_args, arg)}'
        logging.info(msg)
    
    return terminal_args

def if_ref_alt_allele_match(queried_ref, queried_alt, reference_ref, reference_alt):
    '''
    Check if ref and alt alleles of two SNPs match.
    Necessary for multiallelic sites. A single postion may map multiple rsIDs.
    
    There are few situations to consider when determine the flag:
    1. Ref and alt are the same: 1:1000:A:C and 1:1000:A:C
    2. Flipped allele: 1:1000:A:C and 1:1000:C:A
    3. Flipped strand: 1:1000:A:C and 1:1000:T:C
    4. Different ways to represent inserstion and deletion: 1:1000:AAC:C, 1:1000:I:D !!!!!!!!!!!!! 1:1000:AA:- !!!!!!!!!!!!!!!!!!!!!!
    
    Params:
    - queried_ref, queried_alt: ref and alt alleles of the queried SNP
    - reference_ref, reference_alt: ref and alt alleles of a SNP at the same locus in the reference file
    Return:
    - True (match) or False (different)
    '''
    # Make sure everything is in uppercase
    queried_ref, queried_alt = queried_ref.upper(), queried_alt.upper()
    reference_ref, reference_alt = reference_ref.upper(), reference_alt.upper()
    
    dict_flip_strand = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    flag = False
    # ########## 1. Ref and alt are the same: 1:1000:A:C and 1:1000:A:C ##########
    if queried_ref==reference_ref and queried_alt==reference_alt: # match
        flag = True
    
    # ########## 2. Flipped allele: 1:1000:A:C and 1:1000:C:A ##########
    elif queried_ref==reference_alt and queried_alt==reference_ref: # flipped match
        flag = True
    elif ',' in reference_ref or ',' in reference_alt:
        # Handle multiallelic site
        # Eg: chr1:75631396A:G in the reference file is: rs211737, ref=A, alt=C,G,T
        flag = if_multiallelic_match(queried_ref, queried_alt, reference_ref, reference_alt)
    else:
        # ########## 3. Flipped strand: 1:1000:A:C and 1:1000:T:C ##########
        # Try to flip strand
        queried_ref_flip_strand = dict_flip_strand.get(queried_ref)
        queried_alt_flip_strand = dict_flip_strand.get(queried_alt)
        
        if queried_ref_flip_strand==reference_ref and queried_alt==reference_alt:
            # Ref strand flipped
            flag = True
        elif queried_ref_flip_strand==reference_alt and queried_alt==reference_ref:
            # Ref strand flipped, ref/alt allele flipped
            flag = True
        elif queried_ref==reference_ref and queried_alt_flip_strand==reference_ref:
            # Alt strand flipped
            flag = True
        elif queried_ref==reference_alt and queried_alt_flip_strand==reference_ref:
            # Alt strand flipped, ref/alt allele flipped
            flag = True
        elif queried_ref_flip_strand==reference_ref and queried_alt_flip_strand==reference_alt:
            # Both strand flipped
            flag = True
        elif queried_ref_flip_strand==reference_alt and queried_alt_flip_strand==reference_ref:
            # Both strand flipped, ref/alt allele flipped
            flag = True

        # ########## 4. Different ways to represent inserstion and deletion: 1:1000:AAC:C and 1:1000:I:D ##########
        # 'D' and 'I' refer to deletion and insertion
        if queried_ref in ['D', 'I'] and queried_alt in ['D', 'I']:
            if len(reference_ref) > len(reference_alt): # Deletion in reference file
                if queried_ref=='I' and queried_alt=='D':
                    flag = True
            elif len(reference_ref) < len(reference_alt): # Insertion in reference file
                if queried_ref=='D' and queried_alt=='I':
                    flag = True
    return flag

def if_multiallelic_match(queried_ref, queried_alt, reference_ref, reference_alt):
    '''
    Further map a multiallelic site
    Example: chr1:75631396A:G in the reference is: rs211737, ref=A, alt=C,G,T
            In this example, A>C, A>G, A>T all have the same rsid
    Params:
    - queried_ref, queried_alt: ref and alt alleles of the queried SNP
    - reference_ref, reference_alt: ref and alt alleles of a SNP at the same locus in the reference file
    Return:
    - True (match) or False (different)
    '''
    # Make sure everything is in uppercase
    queried_ref, queried_alt = queried_ref.upper(), queried_alt.upper()
    reference_ref, reference_alt = reference_ref.upper(), reference_alt.upper()
    
    dict_flip_strand = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    flag = False
    if ',' in reference_ref:
        alleles = reference_ref.split(',')
        for ref_ref in alleles:
            # Yeah I used recursion by calling if_ref_alt_allele_match,
            # but there there at most one layer of recursion, so this should not cause any issue
            flag = if_ref_alt_allele_match(queried_ref, queried_alt,
                                           reference_ref=ref_ref,
                                           reference_alt=reference_alt)
            if flag: # Stop if find a match
                break
    elif ',' in reference_alt:
        alleles = reference_alt.split(',')
        for ref_alt in alleles:
            flag = if_ref_alt_allele_match(queried_ref, queried_alt,
                                           reference_ref=reference_ref,
                                           reference_alt=ref_alt)
            if flag: # Stop if find a match
                break
    return flag
    
    
def find_single_snp(snp_id, vcf_reference):
    '''
    Given a SNP id in format of chr:position:ref:alt, return the rsID using local VCF reference
    Params:
    - snp_id: chr:position:ref:alt.
              Style of chr must match that in the reference file, ie. both chr22 or 22
    - vcf_Reference: Reference file
    Return:
    - rsID, or 'NA's if rsID is not found
    '''
    # Eg. tabix GCF_000001405.40.gz NC_000001.11:10001-10002]
    chr_num, pos, ref, alt = snp_id.split(':')
    arg = f'tabix {vcf_reference} {chr_num}:{pos}-{pos}'
    result = subprocess.run(arg, shell=True, text=True, capture_output=True).stdout
    if result == '':
        rsid = snp_id
    else:
        # Need to compare ref and alt allele if this is a multiallelic site
        result = result.strip().split('\n')
        if len(result) == 1:
            # Assuming the 3rd column (ID column) contains rsIDs in reference
            _, _, rsid, reference_ref, reference_alt, _ = result[0].split('\t', maxsplit=5)
            return rsid, reference_ref, reference_alt
        else:
            # Multiallelic site found            
            for val in result:
                _, _, rsid, reference_ref, reference_alt, _ = val.split('\t', maxsplit=5)
                if if_ref_alt_allele_match(ref, alt, reference_ref, reference_alt):
                    return rsid, reference_ref, reference_alt # Find a match
    return 'NA', 'NA', 'NA' # No match was found
        
if __name__=='__main__':
    # Terminal mode
    args = process_args()

    # Format of chromosome number is different in NCBI
    # Reference
    # - (GRCh37): https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.25/
    # - (GRCh38): https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
    if args.build==37: # GRCh37
        dict_chr = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11',
                    '4': 'NC_000004.11', '5': 'NC_000005.9', '6': 'NC_000006.11',
                    '7': 'NC_000007.13', '8': 'NC_000008.10', '9': 'NC_000009.11',
                    '10': 'NC_000010.10', '11': 'NC_000011.9', '12': 'NC_000012.11',
                    '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9',
                    '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9',
                    '19': 'NC_000019.9', '20': 'NC_000020.10', '21': 'NC_000021.8',
                    '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9',
                    'MT': 'NC_012920.1'}
    else: # GRCh38
        dict_chr = {'1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12',
                    '4': 'NC_000004.12', '5': 'NC_000005.10', '6': 'NC_000006.12',
                    '7': 'NC_000007.14', '8': 'NC_000008.11', '9': 'NC_000009.12',
                    '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
                    '13': 'NC_000013.11', '14': 'NC_000014.9', '15': 'NC_000015.10',
                    '16': 'NC_000016.10', '17': 'NC_000017.11', '18': 'NC_000018.10',
                    '19': 'NC_000019.10', '20': 'NC_000020.11', '21': 'NC_000021.9',
                    '22': 'NC_000022.11', 'X': 'NC_000023.11', 'Y': 'NC_000024.10',
                    'MT': 'NC_012920.1'}
    
    dict_delimiter = {'tab':'\t', ',':',', 'space':' '}
    args.delimiter = dict_delimiter[args.delimiter]

    # Read input file line by line, find rsID based on SNP ID chr:pos:ref:alt
    fn_output = os.path.join(args.output_path, args.output_prefix+'.converted')
    fh_output = open(fn_output, 'w')
    if not args.inplace:
        fh_output.write('queried_snp_id\trsID\tref_allele_in_reference\talt_allele_in_reference\n')
    with open(args.input_fn) as fh:
        if args.input_header:
            headers = fh.readline() # Skip header line
            if args.inplace:
                # Re-use the headers as the original input file if args.inplace is true
                # Cannot use ',' when write inplace, since the alleles might contain '',
                fh_output.write('\t'.join(headers.strip().split(args.delimiter))+'\trsID\tref_allele_in_reference\talt_allele_in_reference\n')
        count = 0
        for line in fh:
            snp_id = line.strip().split(args.delimiter)[args.snp_col_indx]
            chr_num, pos, ref, alt = snp_id.split(args.snp_id_delimiter)
            if 'chr' not in chr_num.lower(): # 1:1000:A:G
                chr_num = dict_chr[chr_num]
            elif 'CHR' in chr_num: # CHR1:1000:A:G
                chr_num = dict_chr[chr_num.split('CHR')[-1]]
            else: # chr1:1000:A:G
                chr_num = dict_chr[chr_num.split('chr')[-1]]
            query_snp_id = ':'.join([chr_num, pos, ref, alt])
            returned_id, returned_ref, returned_alt = find_single_snp(query_snp_id, args.reference)
            if not args.inplace:
                fh_output.write(f'{snp_id}\t{returned_id}\t{returned_ref}\t{returned_alt}\n')
            else:
                fh_output.write('\t'.join(line.strip().split(args.delimiter))+f'\t{returned_id}\t{returned_ref}\t{returned_alt}\n')
            count += 1
            if count%10==0:
                print(f'\r# N SNPs processed: {count}    ', end='', flush=True)
                fh_output.flush()
    fh_output.close()
    print(f'\n# Done    ')


