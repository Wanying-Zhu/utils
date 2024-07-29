import pandas as pd
import numpy as np

# Process phecode
fn_phecode = '/data100t1/share/BioVU/phenos/official_release_0619/Below_PheCodes_20190531.csv'
print('# Load phencode file:', fn_phecode)
df_phecode = pd.read_csv(fn_phecode, dtype=str)
# print(df_phecode.head())
print('# File size:', df_phecode.shape)

df_phencode_count = df_phecode.groupby(by=['GRID', 'PHEWAS_CODE']).count().reset_index()
lst_dfs, c = [], 0
grids = df_phencode_count['GRID'].unique()
print('# Count number of phecodes and reformat')
for phecode, df in df_phencode_count.groupby('PHEWAS_CODE'):
    c += 1
    # PHEWAS_CODE_DESCRIPTION or PHEWAS_DATE contains count of that phecode
    df.rename(columns={'PHEWAS_CODE_DESCRIPTION':phecode}, inplace=True)
    df.set_index(keys='GRID', inplace=True)
    lst_dfs.append(df.reindex(index=grids, fill_value=0)[phecode])
    print(f'\r# - Process {c}    ', end='', flush=True)
df_merged = pd.concat(lst_dfs, axis=1)

output_fn = '20240729_phecode_counts_BioVU'
print('\n# Save to output:', output_fn+'.txt')
df_merged.sort_index(inplace=True)
df_merged.to_csv(output_fn+'.txt', sep='\t')

print('\n# Save binary to output:', output_fn+'.binary.txt')
df_merged[:] = np.where(df_merged>0, 1, 0) # Replace with binary values
df_merged.to_csv(output_fn+'.binary.txt', sep='\t')
