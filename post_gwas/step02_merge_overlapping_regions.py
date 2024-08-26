# merge regions if they overlap
# Load GWAS region output from 01_find_GWAS_regions.py

def find_pairs_of_regions(df_regions,
                          colname_id='SNP',
                          colname_index='region_index',
                          colname_pos='POS',
                          colname_chr='CHR'):
    '''
    Find pairs of indices of regions to be merged. Eg. (1,2), (3,4), (4,5)
    It becomes a graph, with each region as a node and a pair of indices as edge
    Param:
    - df_regions: a DataFrame contains regions and indices to each region
    - colname_id: id column to find duplicate snps. Cannot just use pos and chr due to multiallelic sites
    - colname_index: column header of the region index
    - colname_pos: column header of the position
    - colname_chr: column header of the chromosome
    Return:
    - n_merged: number of regions merged (due to overlap)
    - df_regions_merged: cleaned dataframe
    '''
    # Find duplicate SNPs
    mask = df_regions[colname_id].duplicated(keep=False)
    df_duplicate = df_regions[mask].copy()
    
    # Store pairs (or more than 2) of regions need to be merged
    regions_to_merge = set()
    for snp, df in df_duplicate.groupby(colname_id):
        inx_pair = df[colname_index].unique().tolist()
        # Sort the indices in the pair so that the smaller value comes first
        inx_pair.sort()
        regions_to_merge.add(tuple(inx_pair)) # Set only accepts hashable items
    return regions_to_merge
    

def find_connected_nodes(edges, node0):
    '''
    Given a collection of edges, find all nodes connected to node0
    (Use width first search)
    '''
    edges = edges.copy() # Deep copy to modify without affecting original list
    connceted_nodes = []
    edges_to_pop = [] # track edges to be removed from next iteration
    # Find directly connected nodes first
    for edge in edges:
        if node0 in edge:
            edges_to_pop.append(edge)
            # Find a node directly connected to node0 first
            if node0==edge[0]:
                node1 = edge[1]
            else:
                node1 = edge[0]
            connceted_nodes.append(node1)
            
    for e in edges_to_pop:
        edges.pop(edges.index(e))
    
    # Recursively find nodes directly connected to other nodes
    other_connected_nodes =[]
    for node in connceted_nodes:
        other_connected_nodes += find_connected_nodes(edges, node)
        
    # Remove duplicates and sort
    result = list(set([node0]+connceted_nodes+other_connected_nodes))
    result.sort()
    return result 

def merge_regions(df_regions,
                  colname_id='SNP',
                  colname_index='region_index',
                  colname_pos='POS',
                  colname_chr='CHR'):
    '''
    According to FINEGENE fine mapping pipeline, merge regions if they overlap.
    Process index pairs and merge regions if they overlap
    Eg: Input index pairs with (1,2), (2,3), (4,5), and output (1,2,3), (4,5).
    Label region 1, 2 and 3 with new label "1,2,3"
    
    Param:
    - df_regions: a DataFrame contains regions and indices to each region
    - colname_index: column header of the region index
    - colname_pos: column header of the position
    - colname_chr: column header of the chromosome
    Return:
    - output: Regions with overlapping regions merged
    '''
    print('# Get regions to be merged, create new labels')
    regions_to_merge = list(find_pairs_of_regions(df_regions=df_regions,
                                                  colname_id=colname_id,
                                                  colname_index=colname_index,
                                                  colname_pos=colname_pos,
                                                  colname_chr=colname_chr))
    
    lst_nodes = df_regions[colname_index].unique().tolist() # Get nodes (ie. index of each region)
    lst_nodes.sort()
    dict_new_reion_index = dict() # New label to reach region
    for node in lst_nodes:
        dict_new_reion_index[node] = find_connected_nodes(edges=regions_to_merge, node0=node)

    print('# Merge regions and relabel')
    df_regions.drop_duplicates(subset=[colname_id], inplace=True)
    df_regions['merged_region_index'] = df_regions[colname_index].apply(lambda x: ','.join([str(v) for v in dict_new_reion_index[x]]))
    
    # number of regions after merging
    n_after_merged = len(df_regions['merged_region_index'].unique())
    return n_after_merged, df_regions



