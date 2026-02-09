'''
Query OpenTargets to get mapped gene for a list of variants
Use function query_OpenTarget()

Example usage:

import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils/post_gwas/S03_annotation')
from query_OpenTarget import query_OpenTarget

query_OpenTarget(variants=['1_123456_T_C', '2_123456_T_C'],
                 output_prefix='output_file_prefix', 
                 output_path='./')


'''



# Import relevant libraries to make HTTP requests and parse JSON response
import requests
import json

def single_variant_query(variant_id:str):
    '''
    Use GraphQL API to query OpenTarget for a single variant
    Parameters
    ----------
    - variant : str
        variant in the format of 1_123456_T_C. The function will try fliped ref/alt alleles if no result is found with the original input

    Returns
    -------
    - nearest_likely_gene: sort the result by min_distance_from_footprint, the first entry is returned as the nearest likely gene
    - min_distance_from_footprint
    - final_predicted_consequence: such as missense_variant
    - final_biotype: such as protein_coding
    - final_distance_from_tss: distance from transcription start site
    
    Example
    -------
    quesy_OpenTarget('1_123456_T_C')    
    '''

    # Build query string to get general information and genetic constraint and tractability assessments 
    query_string = """
    query VariantEffectPredictorQuery($variantId: String!) {
        variant(variantId: $variantId) {
        id
        transcriptConsequences {
            variantConsequences {
            id
            label
            }
            aminoAcidChange
            uniprotAccessions
            codons
            distanceFromFootprint
            distanceFromTss
            target {
            id
            approvedSymbol
            biotype
            }
            impact
            consequenceScore
            transcriptIndex
            transcriptId
            lofteePrediction
            siftPrediction
            polyphenPrediction
        }
        referenceAllele
        alternateAllele
        }
    }
    """

    # Set variables object of arguments to be passed to endpoint
    variables = {"variantId": variant_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    # print(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)
    # print(api_response)
    
    # Process the API response
    if not api_response['data']['variant']:
        # print('# No result found')
        return None
    else:
        for i, v in enumerate(api_response['data']['variant']['transcriptConsequences']):
            gene_name = v['target']['approvedSymbol']
            biotype = v['target']['biotype']
            distance_from_tss = v['distanceFromTss']
            distance_from_footprint = v['distanceFromFootprint']
            predicted_consequence = v['variantConsequences'][0]['label']

            if i==0:
                min_distance_from_footprint = distance_from_footprint
                nearest_likely_gene = gene_name
                final_predicted_consequence = predicted_consequence
                final_biotype = biotype
                final_distance_from_tss = distance_from_tss
            else:
                if distance_from_footprint < min_distance_from_footprint:
                    min_distance_from_footprint = distance_from_footprint
                    nearest_likely_gene = gene_name
                    final_predicted_consequence = predicted_consequence
                    final_biotype = biotype
                    final_distance_from_tss = distance_from_tss
        return(nearest_likely_gene, final_predicted_consequence, final_biotype, min_distance_from_footprint, final_distance_from_tss)


def query_OpenTarget(variants:list, output_prefix:str, output_path:str):
    '''
    Use GraphQL API to query OpenTarget for a list of variants
    Parameters
    ----------
    - variant : list
        list of variants in the format of [1_123456_T_C]
    - output_prefix: prefix of output file. The code will create two files, one for mapped gene, the other for SNPs without mapped gene
    - output_path: path to save the output files

    Returns
    -------
    Return a list of tuples, each tuple contains below information for the given variant:
    - nearest_likely_gene: sort the result by min_distance_from_footprint, the first entry is returned as the nearest likely gene
    - min_distance_from_footprint
    - final_predicted_consequence: such as missense_variant
    - final_biotype: such as protein_coding
    - final_distance_from_tss: distance from transcription start site
    
    Example
    -------
    query_OpenTarget(['1_123456_T_C', '2_123456_T_C'])    
    '''
    fh_output = open(f'{output_path}/{output_prefix}.OpenTarget.tsv', 'w')
    fh_not_found = open(f'{output_path}/{output_prefix}.OpenTarget.not_found.tsv', 'w')
    fh_output.write('variant_id\tnearest_likely_gene\tpredicted_consequence\tbiotype\tdistance_from_footprint\tdistance_from_tss\n')
    fh_not_found.write('variant_id\n')
    count_not_found, count_found = 0, 0
    for i, variant_id in enumerate(variants):
        res = single_variant_query(variant_id)
        if not res: # If no result is found, try fliped ref/alt alleles
            chr_num, pos, ref, alt = variant_id.split('_')
            res = single_variant_query(f'{chr_num}_{pos}_{alt}_{ref}')
        if not res:
            print(f'\n# - No result found for variant {variant_id}')
            fh_not_found.write(f'{variant_id}\n')
            count_not_found += 1
            continue
        
        # Save returned results to outpout file
        fh_output.write(f'{variant_id}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\t{res[4]}\n')
        count_found += 1
        if len(variants) > 100 and i % 100 == 0:
            print(f'\r# Processing variant {i+1}/{len(variants)} (found/not found: {count_found}/{count_not_found})', end='', flush=True)
        else:
            print(f'\r# Processing variant {i+1}/{len(variants)} (found/not found: {count_found}/{count_not_found})', end='', flush=True)
    fh_output.close()
    fh_not_found.close()