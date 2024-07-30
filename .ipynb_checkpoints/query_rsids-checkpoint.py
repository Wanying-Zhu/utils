'''
Example usage:

sys.path.append('/data100t1/home/wanying/lab_code/utils')
from query_rsids import query_rsids
query_rsids(['rs1234', 'rs4567'], output='result.txt', not_found='not_found.txt')
'''


import requests
from bs4 import BeautifulSoup

def query_rsids(rsids, output='result.txt', not_found='not_found.txt'):
    '''
    Query dbsnp with  rsID and get B37 and B38 positions
    Params
    - rsids: an array of rsIDs to query
    - output: file name of the outputs
    - not_found: a file to store rsIDs that were not found
    '''
    fh_output = open(output, 'w')
    fh_not_found = open(not_found, 'w')
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
            b38, b37, _ = txt.split(')')
            b38 = b38.split('\n(')[0]
            b37 = b37.split('\n(')[0]
            line = f'{rsid}\t{b37}\t{b38}'
            fh_output.write(line+'\n')
        except:
            fh_not_found.write(rsid+'\n')
        print(f'\r# Processing: {i}/{len(rsids)}        ', end='', flush=True)
    print('\n# Done')
    fh_output.close()
    
        
    
