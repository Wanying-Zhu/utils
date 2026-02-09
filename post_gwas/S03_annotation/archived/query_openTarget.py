# Code by Hannah
# Modified by Wanying
# Query OpenTarget genes for variant (V2G)


'''
Example usage:
from query_openTarget import v2g
v2g(variant_ids=['11_61800281_C_A', 'rs61902204'], output_prefix='output')

'''
# To use the code get CHR_POS_REF_ALT in b38 (?) format first, or use rsid
### This script takes a list of variant ids in rsid format (one variant per line) and ###
# 1. Finds the gene related to that rsid depending on Open Targets assigning variants to genes pipeline (aka "V2G")
#       See: https://genetics-docs.opentargets.org/our-approach/data-pipeline
# 2. A) writes a file that captures only the top n (default 100) genes identified from provided variant list
#    B) writes a file listing any variants not found in Open Targets portal (not every variant is annotated in Open Targets)


#Import libraries
import requests #HTTP request
import json #json formatting
import os

def query(variant_id, outfile, notfoundfile, n_genes):
    # outfile = open(os.path.join(output_path, output_prefix+f'.top{n_genes}genelist_OpenTarget.txt'), mode=mode)
    # notfoundfile = open(os.path.join(output_path, output_prefix+f'.variantnotfound_OpenTarget.txt'), mode=mode)
    keep_gene_dict = {}
    #build query string
    query_string = """
    query useSearchToConvertRSIDIntoIDFormat($myVariantID: String!) {
    search(queryString:$myVariantID){
        totalVariants
        variants{
        id
        }
    }
    }
    """

    #Set variable(s) object of arguments to be passed to endpoint
    variable = {"myVariantID": variant_id}

    #Set base URL of Genetics Portal GraphQL API endpoint
    base_url = "https://api.genetics.opentargets.org/graphql"

    #Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variable})
    #print(r.status_code)

    #Transform API response into JSON 
    api_response_as_json = json.loads(r.text)

    #Print first element of JSON response data
    #print(api_response_as_json["data"]["search"]["variants"][0]["id"])

    #Store variant ID as new variable
    refid = api_response_as_json["data"]["search"]["variants"][0]["id"]

    #Now, we want to query our new variant ID (refid) and find its associated functional gene

    #build query string

    query_string = """
    query myQuery ($myrefID: String!) {
    genesForVariant(variantId: $myrefID){
        overallScore
        gene {
        symbol
        }
    }
    }
    """

    #Set object of argument to be passed to endpoint
    variant = {"myrefID": refid}

    #Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variant})
    #print(r.status_code)

    #Transform API response into JSON 
    api_response_as_json = json.loads(r.text)

    #Print JSON response data
    #print(api_response_as_json["data"]["genesForVariant"])
    json_lst = (api_response_as_json["data"]["genesForVariant"])

    curr_score = 0

    for item in json_lst:
        temp_dict = item

        gene_name = (temp_dict["gene"]["symbol"])
        functional_score = (temp_dict["overallScore"])
        if functional_score > curr_score:
            curr_score = functional_score
            keep_gene_name = gene_name

    if len(keep_gene_dict) == n_genes:
        #stop script, we have our list of n_genes (unique) genes !
        notfoundfile.close()
        outfile.close()
        print('# DONE')
        return
    elif keep_gene_name not in keep_gene_dict and len(keep_gene_dict) < n_genes:
            keep_gene_dict[keep_gene_name] = keep_gene_name
            outfile.write(keep_gene_name + '\n')


def v2g(variant_ids:list, output_prefix:str, output_path:str='./',
        n_genes:int=100, mode:str='w'):
    '''
    Query variant to gene (V2G)
    Params:
    - variant_id: SNP in the format of CHR_POS_REF_ALT in B38
    - output_prefix: prefix of the output file (and SNPs not found)
    - output_path: path of the ouput file
    - n_genes: Max number of genes to output
    - mode: mode to write to the output file 'w'=write or 'a'=append
    Return
    Save results to files
    '''
    ## IMPORTANT: variant file must have snps listed from more significant to least significant (low to high pval) ##
    
    #Save functional genes in dictionary
    # keep_gene_dict = {}
    outfile = open(os.path.join(output_path, output_prefix+f'.top{n_genes}genelist_OpenTarget.txt'), mode=mode)
    notfoundfile = open(os.path.join(output_path, output_prefix+f'.variantnotfound_OpenTarget.txt'), mode=mode)
    for i, variant_id in enumerate(variant_ids):
        if 'RS' not in variant_id.upper(): # Also test allele flipped version
            ch, pos, ref, alt = variant_id.split('_')
            variant_id_flipped = f'{ch}_{pos}_{alt}_{ref}'
            
        try:
            try:
                query(variant_id, outfile, notfoundfile, n_genes)
            except: # Try flipped id if original one returned nothing
                try:
                    print('\n# - Test allele flipped ID')
                    query(variant_id_flipped, outfile, notfoundfile, n_genes)
                except:
                    print('\n# - Still no result')
                    raise Exception("No result")
                    
        except:
            notinopentarget = variant_id
            notfoundfile.write(notinopentarget + '\n') #record variant as not found
        print(f'\r# Processed {i+1}/{len(variant_ids)}    ', end='', flush=True)
        
    print(f'\r# Processed {i+1}/{len(variant_ids)}    ')
    notfoundfile.close()
    outfile.close()
    print('# DONE')
