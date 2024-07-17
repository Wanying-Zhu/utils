# This code takes in a list of variants, with chromosome number and positions,
# then search in the gzipped (.gz) genotype data files and extract these snps
# Results are outputted in a file with user supplied file name
# This code does not use pandas so should be more appropriate for very large dataset (like BioVU)
# ############################################# Caveat(s) #############################################
# - Multiallelic site: search by position only may not get the correct SNP for multiallelic sites.
#         Need to check ref and alt alleles as well (if known) using find_variants_multiallelic().
# - Currently does not support searching for dupicate positions
# #####################################################################################################

from xopen import xopen # xopen is faster than gzip
    
def find_variants(lst_pos, output_fn, input_fn, input_col_name='POS', threads=1, verbose=False):
    '''
    This function takes in positions of a list of snp (no duplication),
    then checks the input genotype file and output found SNP into a output file
    Parameters:
        - lst_pos: positions (int) of snps. Must be array-like for iteration
        - output_fn: path and file name to output found SNPs
        - input_fn: path and file name of input genotype file (VCF). Must be a single chromosome
                    Assume all variants are sorted base on position in the input file
        - input_col_name='POS': column name of position in input genotype file. Usually 'POS', could be 'current_pos'
        - threads: multi-threading
        - verbose: whether to print number of snps processed
    Returns:
        - Output found SNPs and genotyeps into output file.
        - Print out any missing SNPs into console
    '''
    lst_pos = list(lst_pos) # Convert lst_pos to a list, so that .sort() and .remove() will work
    lst_pos.sort() # Sort positions to be found
    input_fh = xopen(input_fn, threads=threads)
    output_fh = open(output_fn, 'w')
    # not_found_fh = open(output_fn+'.not_found_pos', 'w') # Write a position into this file if not found in genotype file
    line = input_fh.readline().strip()

    while line[0:2] == '##': # Read trhough header lines
        line = input_fh.readline().strip()

    output_fh.write(line) # Now reach the column header line, write it into output file
    output_fh.write('\n')

    line = input_fh.readline().strip()
    count = 0
    while len(lst_pos) != 0 and line != '':
        current_pos = int(line.split(maxsplit=2)[1]) # get position of current row

        if current_pos<=lst_pos[0]:
            if current_pos==lst_pos[0]: # If found a variant, write into output file
                output_fh.write(line)
                output_fh.write('\n')
                lst_pos.remove(current_pos) # Update list of positions

            line = input_fh.readline().strip() # Put code below here, since if current_pos > lst_pos[0], does not need to read a new line
            if verbose:
                count += 1
                if count % 50000 == 0:
                    print('.', end='', flush=True)
                if count % 2500000 == 0:
                    print(count, 'rows processed, now at position', current_pos)

        elif current_pos > lst_pos[0]: # If this variant is not in the genotype file, skip it and look for the next one
            # not_found_fh.write(str(lst_pos[0])+'\n')
            lst_pos.remove(lst_pos[0])

    input_fh.close()
    output_fh.close()
    # not_found_fh.close()


def find_variants_multiallelic(lst_pos_ref_alt, output_fn, input_fn,
                               input_col_names=['POS', 'REF', 'ALT'],
                               threads=1, verbose=False):
    '''
    This function takes in positions of a list of snp (no duplication),
    then checks the input genotype file and output found SNP into a output file
    Parameters:
        - lst_pos_ref_alt: Tuples of (position (int), ref_allele (str), alt_allele (str)) of snps.
                           Must be array-like for iteration.
                           Must be ordered by position, then ref allele and alt allele.
        - output_fn: path and file name to output found SNPs
        - input_fn: path and file name of input genotype file (VCF). Must be a single chromosome
                    Assume all variants are sorted base on position in the input file
        - input_col_names=['POS', 'REF', 'ALT']': column names of position, ref allele and alt allele in input genotype file.
        - threads: multi-threading
        - verbose: whether to print number of snps processed
    Returns:
        - Output found SNPs and genotyeps into output file.
        - Print out any missing SNPs into console
    '''
    # Convert input to a list for pop remove
    lst_snps = list(lst_pos_ref_alt)
    input_fh = xopen(input_fn, threads=threads)
    output_fh = open(output_fn, 'w')
    
    not_found_fh = open(output_fn+'.not_found_pos', 'w') # Write a position into this file if not found in genotype file
    line = input_fh.readline().strip()

    while line[0:2] == '##': # Read through header lines
        line = input_fh.readline().strip()

    output_fh.write(line + '\n') # Now reach the column header line, write it into output file
    header_lst = line.split(maxsplit=9) # Get indices of position, ref allele and alt allele
    indx_pos, indx_ref, indx_alt = header_lst.index(input_col_names[0]), header_lst.index(input_col_names[1]), header_lst.index(input_col_names[2])
    line = input_fh.readline().strip()
    count, count_not_find = 0, 0
    while len(lst_snps) != 0 and line != '':
        tmp = line.split(maxsplit=9)
        current_pos, ref_allele, alt_allele = int(tmp[indx_pos]), tmp[indx_ref], tmp[indx_alt]

        # Position, ref allele and alt allele to look for in the input list
        pos_to_find, ref_allele_to_find, alt_allele_to_find = lst_snps.pop(0)
        find = False # Track if a SNP is found
        # line = input_fh.readline().strip()
        while current_pos<=pos_to_find and not find and line != '':
            tmp = line.split(maxsplit=9)
            current_pos, ref_allele, alt_allele = int(tmp[indx_pos]), tmp[indx_ref], tmp[indx_alt]

            if verbose:
                count += 1
                if count % 50000 == 0:
                    print('.', end='', flush=True)
                if count % 2500000 == 0:
                    print(count, 'rows processed, now at position', current_pos)

            # If found a variant, write it into output file
            if current_pos==pos_to_find and ref_allele==ref_allele_to_find and alt_allele==alt_allele_to_find:
                find = True
                output_fh.write(line + '\n')
                line = input_fh.readline().strip()
                #
                # # If the SNP is found, no need to loop anymore. Move on to the next input SNP
                break
            line = input_fh.readline().strip()


        # If exit the while loop then (1) SNP at pos_to_find is already found
        #   or (2) current_pos>pos_to_find
        #   or (3) Reach to end of the file
        # Keep looking for the next one until position list is empty
        if not find:
            count_not_find += 1
            not_found_fh.write(f'{pos_to_find}\t{ref_allele_to_find}\t{alt_allele_to_find}\n')


    input_fh.close()
    output_fh.close()
    if verbose: print(f'\n#Number of SNPs not found: {count_not_find}\n#DONE\n')
    not_found_fh.close()