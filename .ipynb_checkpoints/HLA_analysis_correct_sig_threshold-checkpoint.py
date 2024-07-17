# This script is used to correct significance threshold of
# amino acid analysis reulst
# Need AA_OR.txt and AA_chisq.txt

import pandas as pd

chiq_fn = 'AA_chisq.txt'
OR_fn = "AA_OR.txt"
df_chiq = pd.read_csv(chiq_fn, sep='\t')
df_or = pd.read_csv(OR_fn, sep='\t')

df_chiq = df_chiq[['Locus', 'Position', 'df']]
df_chiq['merge'] = df_chiq['Locus']+'_'+df_chiq['Position']
df_or['merge'] = df_or['Locus']+'_'+df_or['Position']
df_merge = df_or.merge(df_chiq, left_on='merge', right_on='merge')

df_merge['adj_sig'] = 0.05/(df_merge['df']+1)
df_merge['filter'] = (df_merge['p.value']<df_merge['adj_sig'])

tmp_df = df_merge.loc[df_merge['filter']==1]
tmp_df[['Locus_x', 'Position_x', 'Residue', 'OR', 'CI.lower', 'CI.upper','p.value', 'sig']].to_csv('out.txt', sep='\t', index=False)

