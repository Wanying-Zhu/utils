def find_ensgid_and_pos(gene_name='NT5C2', build='hg38'):
    '''
    Search Genecard website for Ensmbl ID and position of a given gene
    Parameters
    - gene_name: Name of the gene to be searched
    - build: Defines build of the returned position, default is 'hg38' (GRCh38/hg38). All other velues will return position in build GRCh37/hg19
    Return
    - ensg: Ensembl ID of the gene, such as ENSG00000076685
    - pos: position of the gene, such as chr10:103087185-103277605
    '''
    from bs4 import BeautifulSoup
    import re
    import urllib.request

    genecard_url = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='
    req = urllib.request.Request(genecard_url+gene_name, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as response:
        html_bytes = response.read()
    html = html_bytes.decode("utf-8") # Decode the bytes to a string using UTF-8
    html_text = BeautifulSoup(html, 'html.parser').get_text() # Convert returned html to plain text format
    
    m_ensgid = re.search('Ensembl:.*', html_text) # Find Ensembl ID
    
    if build == 'hg38':
        # Use DOTALL to make . match newline character as well
        m_pos = re.search('Latest Assembly.*\(GRCh38/hg38\)', html_text, flags=re.DOTALL) # Find position
    else:
        m_pos = re.search('Previous Assembly.*\(GRCh37/hg19\)', html_text, flags=re.DOTALL) # Find position
        
    try:
        ensg = m_ensgid.group(0).split()[-1] # Eg: ['Ensembl:', 'ENSG00000076685']
        pos = m_pos.group(0).split()[-2] # Eg: ['Latest', 'Assembly', 'chr10:103,087,185-103,277,605', '(GRCh38/hg38)']
        chr_num = pos.split(':')[0] # Get chromosome number
        start = pos.split(':')[-1].split('-')[0].replace(',', '') # Remove ',' in numbers
        end = pos.split(':')[-1].split('-')[-1].replace(',','')
        pos = f'{chr_num}:{start}-{end}'
        return ensg, pos # Return matched string
    except:
        print('\tWARNING: Ensembl ID or position is not find')
        print('\t- Search result of ensembl ID:', m_ensgid)
        print('\t- Search result of position:', m_ensgid)
