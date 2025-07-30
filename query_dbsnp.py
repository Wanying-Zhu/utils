'''
! Recommend only use this code to lookup chr:pos of rsID
! Use code in snpid2rsid to map SNP ID to rsID

Example usage:
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils')
from query_dbsnp import query_rsids, query_positions
query_rsids(['rs1234', 'rs4567'], output='result.txt', not_found='not_found.txt')

query_positions(['5:11893', '5:11963'], output='result.txt', not_found='not_found.txt')
query_positions(['8:75707664', '8:75707678'], output='result.txt', not_found='not_found.txt')
'''


import requests
from bs4 import BeautifulSoup

def query_rsids(rsids, output='result.txt', not_found='not_found.txt', flush=False, mode='w'):
    '''
    Query dbsnp with  rsID and get B37 and B38 positions
    Params
    - rsids: an array of rsIDs to query
    - output: file name of the outputs
    - not_found: a file to store rsIDs that were not found
    - flush: if ture flush result to output file immediately
    - mode: 'w'=write, 'a'=append to (existing) output file
    '''
    fh_output = open(output, mode=mode)
    fh_not_found = open(not_found, mode=mode)
    fh_output.write('rsID\tGRCh37\tGRCh38\n')
    fh_not_found.write('rsID\n')
    for i, rsid in enumerate(rsids):
        # Define the URL you want to query
        url = 'https://www.ncbi.nlm.nih.gov/snp/?term=' + rsid
        
        # Send a GET request to the URL
        response = requests.get(url)
        
        try:
            txt = None # Reset txt
            # Parse the HTML content
            soup = BeautifulSoup(response.text, 'html.parser')
            for tag in soup.find_all('dd'):
                if 'GRCh' in tag.text:
                    txt = tag.text
                    break
            b38, b37, _ = txt.split(')')
            b38 = b38.split('\n(')[0]
            b37 = b37.split('\n(')[0]
            line = f'{rsid}\t{b37}\t{b38}'
            fh_output.write(line+'\n')
            if flush:
                fh_output.flush()
        except:
            fh_not_found.write(rsid+'\n')
            if flush:
                fh_not_found.flush()
        print(f'\r# Processing: {i+1}/{len(rsids)}        ', end='', flush=True)
    print('\n# Done')
    fh_output.close()

def query_positions(positions, output='result.txt', not_found='not_found.txt', flush=False, mode='w'):
    '''
    Query dbsnp with chr_number:position and get rsID, B37 and B38 positions
    Params
    - positions: an array of positions to query. Such as [1:1234, 2:3456]
    - output: file name of the outputs
    - not_found: a file to store positions that were not found
    - flush: if ture flush result to output file immediately
    - mode: 'w'=write, 'a'=append to (existing) output file
    '''
    fh_output = open(output, mode=mode)
    fh_not_found = open(not_found, mode=mode)
    fh_output.write('ID\trsID\tGRCh37\tGRCh38\n')
    fh_not_found.write('ID\n')
    for i, snp in enumerate(positions):
    # Define the URL you want to query
        chr_num, pos = snp.split(':')
        url = 'https://www.ncbi.nlm.nih.gov/snp/?term=' + chr_num + '%3A' + pos
        # Send a GET request to the URL
        response = requests.get(url)
        
        try:
            # Parse the HTML content
            txt, rsid = None, 'NA' # Reset txt and rsid
            soup = BeautifulSoup(response.text, 'html.parser')
            for tag in soup.find_all('dd'):
                if 'GRCh' in tag.text:
                    txt = tag.text
                    break
            b38, b37, _ = txt.split(')')
            b38 = b38.split('\n(')[0]
            b37 = b37.split('\n(')[0]
    
            for tag in soup.find_all('a'): # Get rsID
                if 'snp/rs' in tag.get('href'):
                    rsid = tag.text
                    break
            line = f'{snp}\t{rsid}\t{b37}\t{b38}'
            fh_output.write(line+'\n')
            if flush:
                fh_output.flush()
        except:
            fh_not_found.write(snp+'\n')
        print(f'\r# Processing: {i+1}/{len(positions)}        ', end='', flush=True)
    print('\n# Done')
    fh_output.close() 

if __name__ == '__main__': # If run as a script
    import argparse
    import sys
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fn', help='Input file name. One SNP per line', type=str)
    parser.add_argument('--output_path', type=str, default='./')
    parser.add_argument('--output_prefix', type=str, default='output')
    parser.add_argument('--snp_indx', help='Index of the SNP id column. Default is the first column (index=0)', type=int, default=0)
    parser.add_argument('--delimiter', help='Delimiter of the file. Default is whitespace (tab or space)', type=str, default='whitespace')
    parser.add_argument('--batch_size', help='Process SNPs by batch size (if need to convert a lot of SNPs). Default (-1) is to load all at once',
                        type=int, default=-1)

    args = parser.parse_args()
    if args.delimiter == 'whitespace': # None is to split on whitespace in string.split
        args.delimiter = None

    # Record script used
    cmd_used = 'python ' + ' '.join(sys.argv)
    for arg in vars(terminal_args):
        cmd_used += f' --{arg} {getattr(terminal_args, arg)}'
        print('# Call used:')
        print(cmd_used)
    
    print('\n# Query dbSNP:')
    fh_snp = open(args.input_fn)
    lst_snps = []
    count = 0
    for line in fh_snp:
        if line.strip() != '':
            snp = line.strip().split(sep=args.delimiter)[args.snp_indx]
            lst_snps.append(snp)
            count += 1
            if agrs.batch_size != -1:
                # Process by bath if needed
                pass
    print('# SNPs loaded:', len(lst_snps))
    
    fh_snp.close()
    
