
import os
import pandas as pd

# This file merge individual files (with OR and SE from single and meta analysis) together
dir = '/data100t1/home/wanying/lab_code/utils/output/preprocessing_meta_and_single/'
lst_fn = []


# get file names of preprocessed outputs
for root, dirs, files in os.walk(dir):
    lst_fn = files

lst_pheno = ['hdl', 'ldl', 't2d', 'tc', 'tg']   # Lauren did analysis on these 5 phenotypes
pheno = lst_pheno[4]

flag_first_df = True    # If current dataframe is the first one, since only need to keep meta result once
df_merged = pd.DataFrame()  # Empty dataframe for final output

# Read through all files and merge OR and SE of single analyses
for fn in lst_fn:

    if fn.split('_')[0] == pheno:   # Process by phenotypes
        fn = dir + fn
        df = pd.read_csv(fn, sep='\t', dtype='str')
        if flag_first_df:   # If this is the first df, then copy everything
            df_merged = df.copy()
            flag_first_df = False
        else:
            df_merged = df_merged.merge(df.iloc[:, [0, 3, 4]],
                                        left_on='MarkerName', right_on='MarkerName')

output_fn = '/data100t1/home/wanying/lab_code/utils/output/merged_OR_SE/' + pheno + '_OR_SE.txt'
df_merged.to_csv(output_fn, sep='\t', index=False, na_rep='NA')

# Add two columns to the end, ie. OR and Se from meta-analysis




