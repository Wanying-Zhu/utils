import gzip
def get_dosage(vcf_fn = '/data100t1/home/wanying/ImageVU/PRS/biovu_prs/merged_biovu_eur_grs1_all_snp.vcf',
               index_format_col = 8, dosage_field = 'DS', cols_to_save='all', output_fn = None):
    '''
    Parameters:
    - vcf_fn: name of input vcf files, can be either compressed (vcf.gz) or uncompressed (.vcf)
    - index_format_col: default is 8. Zero-based index of FORMAT column
    - dosage_field: default is 'DS'. Name of dosage field stated in FORMAT
    - cols_to_save: Columns index to include in the output, default keep all columns. Eg. cols_to_save = [0, 1, 3, 4]
    - output_fn: will use [vcf_fn].dosage if not provided
    
    Return:
    - Output dosage, save in output_fn
    
    '''
    if output_fn is None:
        output_fn = f'{vcf_fn}.dosage'
    # Extract dosage from vcf adn save to output file
    # vcf_fn = '/data100t1/home/wanying/ImageVU/PRS/biovu_prs/merged_biovu_eur_grs1_all_snp.vcf'
    # index_format_col = 8 # Index of FORMAT column
    # dosage_field = 'DS' # Name of dosage field
    # output_fn = f'{vcf_fn}.dosage'
    output_fh = open(output_fn, 'w')
    # fh = open(vcf_fn)
    # line = fh.readline()
    try:
        fh = open(vcf_fn)
        line = fh.readline()
    except:
        fh = gzip.open(vcf_fn, 'rt')
        line = fh.readline()
    
    header = '' # Save header line to write
    while line[0] == '#': # Skip header line
        header = line
        line = fh.readline()
    line = line.strip()

    header = header.strip()
    tmp = header.split(maxsplit=index_format_col+1)
    # Write column header to output file
    columns = tmp[:-1]
    id_columns= tmp[-1] # Sample ids
    if cols_to_save != 'all':
        columns = [columns[x] for x in cols_to_save]
    headerline = '\t'.join(columns) + '\t' + id_columns
    output_fh.write(headerline)

    # Find index of Dosage field in FORMAT (usually DS field)
    tmp = line.split(maxsplit=index_format_col+1)
    index_dosage = tmp[index_format_col].split(':').index(dosage_field)
    while line != '':
        tmp = line.split(maxsplit=index_format_col+1)
        columns = tmp[:-1] # Other columns
        if cols_to_save != 'all':
            columns = [columns[x] for x in cols_to_save]

        genotype_data = tmp[-1]
        tmp_genotype = genotype_data.split()
        dosage_vals = [] # Store dosage values
        for val in tmp_genotype:
            dosage_vals.append(val.split(':')[index_dosage])
        output_line = '\n'+'\t'.join(columns) + '\t' + '\t'.join(dosage_vals)
        output_fh.write(output_line)
        line = fh.readline().strip()
        
    fh.close()
    output_fh.close()
