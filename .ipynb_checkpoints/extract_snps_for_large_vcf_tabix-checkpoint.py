# If VFC is tabix indexed, this should be the preferred way
# In fact it might be better to tabix the file whenever possible, save time in SNP lookups

# This code takes in a list of variants, with chromosome numbers and positions,
# then search in the gzipped (.gz) genotype data files and extract these snps
# Results are outputted in a file with user supplied file name

import subprocess

def find_variants(lst_pos, output_fn, input_fn, mode='w', keep_multiallelic=True, bgzip=True, flush=False):
    '''
    This function takes in positions of a list of snp (allow duplicates),
    then checks the input genotype file and output found SNP into an output file
    Parameters:
        - lst_pos: positions (int) of snps. Must be list-like for iteration.
                   For example: chr1:1234
        - output_fn: path and file name to output found SNPs
        - input_fn: path and file name of input genotype file. Must be a single chromosome
                    Assume all variants are sorted base on position in the input file
        - mode: output mode. 'w'=(over)write, 'a'=append
        - keep_multiallelic: Default is to output all SNPs of multiallelic sites.
                            Will eclude multiallelic SNPs if set to false
        - bgzip: bgzip filtered vcf if True
        - flush: actively write to output file (good long procedure, will still have simething when the server crashes)
    Returns:
        - Output found SNPs and genotypes into output file.
        - Save any missing SNPs into output_fn+'.not_found_pos'
    '''
    
    print('\n#', '#' * 25, 'Start SNP lookups', '#' * 25)
    lst_pos = list(lst_pos) # Convert lst_pos to a list, so that .sort() and .remove() will work
    output_fh = open(output_fn, mode)
    not_found_fh = open(output_fn+'.not_found_pos', 'w') # Write a position into this file if not found in genotype file
    if not keep_multiallelic:  # Skip multiallelic sites but save them for future references
        multiallelic_fh = open(output_fn + '.multiallelic_skipped', 'w')
    total, c = len(lst_pos), 0 # Track progress
    headers = subprocess.run(f'tabix -H {input_fn}', shell=True, text=True, capture_output=True).stdout # Get header lines
    output_fh.write(headers)
    
    c_multiallelic, c_missing = 0, 0 # count number of multiallelic sites, number of SNPs not found
    if not keep_multiallelic:
        print('# - Skip multiallelic sites ')
        
    for pos in lst_pos:
        pos = pos + '-' + pos.split(':')[1] # Reformat position to search by tabix
        cmd = f'tabix {input_fn} {pos}'
        result = subprocess.run(cmd, shell=True, text=True, capture_output=True).stdout
        if not keep_multiallelic: # ignore multiallelic sites
            if len(result.strip().split('\n'))>1: # if returned result contains more than one line
                c_multiallelic += 1
                c += 1
                multiallelic_fh.write(pos.split('-')[0]+'\n') # Save skipped sites for future reference
                if flush: multiallelic_fh.flush()
                continue
        if len(result) == 0: # If SNP not found
            not_found_fh.write(pos.split('-')[0]+'\n')
            if flush: not_found_fh.flush()
            c_missing += 1
        output_fh.write(result)
        if flush: output_fh.flush()
        c += 1
        if (total>=10000 and c%100 == 0) or (total<10000):
            print(f'\r# - Number of SNPs checked: {c}/{total}; Not found: {c_missing}; Multiallelic SNPs skipped: {c_multiallelic}    ',
                  flush=True, end='')
    print(f'\r# - Number of SNPs checked: {c}/{total}; Not found: {c_missing}; Multiallelic SNPs skipped: {c_multiallelic}    ')
    output_fh.close()
    not_found_fh.close()
    if not keep_multiallelic: multiallelic_fh.close()
    
    if bgzip:
        print('\n# - bgzip output file')
        subprocess.run(f'bgzip {output_fn}', shell=True)
    print('\n#', '#'*25, 'DONE', '#'*25)


# Example call
# find_variants(['chr10:10537', 'chr10:10574'],
#               input_fn='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs/max_unrelated_set_chr10.vcf.gz',
#               output_fn='/data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_sample_vcfs/202312_redo_training_vcfs/out.txt',
#               mode='w',
#               keep_multiallelic=False)