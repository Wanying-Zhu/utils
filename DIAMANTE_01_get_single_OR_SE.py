# This code is used to extract OR values of DIAMANTE output files
# (in order to make forest plot)

import pandas as pd

# ---------------- Helper functions ---------------

# Get individual study names and file names went in to meta analysis from each study
def get_file_name(fn, dir_single_analysis_individual_file = '/vgipiper05/DIAMANTE/'):
    """
    Help on method get_file_name:
    get_file_name(fn, dir_single_analysis_individual_file = '/vgipiper05/DIAMANTE/')

    This function gets individual study names and file names that went in to meta analysis from each study

    Parameters
    ----------
    - fn: such as '/vgipiper05/DIAMANTE/mrmega_hdl_final_0820.in'
         A file name of the file that contains names of single analysis result files that went into meta analysis
    - dir_single_analysis_individual_file: location where folders of each single analysis are stored
                                           Default is set to '/vgipiper05/DIAMANTE/'

    Returns
    -------
    - lst_study_names: a list of study names (abbreviations)
    - dict_study_dir_fn: a dictionary of study_name:file_name (file name of single analysis result files that went into meta analysis)
    """
    lst_study_names = []  # A list to store study names
    dict_study_dir_fn = dict()  # A dictionary to store single analysis file names of each study

    # There might be more than one files from the same study
    # If so, they will be renamed as: study, study_1, study_2, study_3, etc.
    last_item = ''
    duplication_count = 1
    with open(fn, 'r') as fh:
        line = fh.readline().strip()
        while line != '':
            print(line)
            study_name = line.split('/')[1]
            if study_name not in lst_study_names:
                lst_study_names.append(study_name)
            else:
                if study_name == last_item:
                    duplication_count = duplication_count + 1
                else:
                    duplication_count = 1
                last_item = study_name
                study_name = study_name+'_'+str(duplication_count)
                lst_study_names.append(study_name)
            dict_study_dir_fn[study_name] = dir_single_analysis_individual_file + line[3:]
            line = fh.readline().strip()
    return lst_study_names, dict_study_dir_fn


# UNFINISHED !!!
# Get OR values of a variant from a single analysis result
def getOR(variant, fn_single_analysis):
    """
    Help on method getOR:
    getOR(variant, fn_single_analysis)

    This function takes in a variate rsID, and search all single analysis files for OR

    Parameters
    ----------
    - variant: rsID of the variant
    - fn_single_analysis: file name of the single analysis went into meta-analysis,
                          with absolute path to the file

    Returns
    -------
    OR of the variant
    """
    with open(fn_single_analysis, 'r') as fh:
        line = fh.readline().strip()
        while line != '':
            lst_tmp = line.split()
            line = fh.readline().strip()


# ---------------- End of helper functions ---------------


dir_meta_analysis = '/vgipiper05/DIAMANTE/mrmega/'  # Directory to meta-analysis results
dir_single_analysis_input = '/vgipiper05/DIAMANTE/mrmega/mrmega_input_files/' # Directory to individual analysis results

# File names of top hit variants in meta-analysis
lst_fn_meta_analysis = ['hdl_0820.tophits',
                        'ldl_0820.tophits',
                        't2d_0820.tophits',
                        'tc_0820.tophits',
                        'tg_0820.tophits']

# File names of names of single analysis results that went into meta-analysis
lst_fn_single_analysis = ['mrmega_hdl_final_0820.in',
                          'mrmega_ldl_final_0820.in',
                          'mrmega_t2d_final_0820.in',
                          'mrmega_tc_final_0820.in',
                          'mrmega_tg_final_0820.in']

column_title = ['MarkerName', 'Chromosome', 'Position', 'EA', 'NEA',
                'EAF', 'Nsample', 'Ncohort', 'Effects', 'beta_0',
                'se_0', 'beta_1', 'se_1', 'chisq_association', 'ndf_association',
                'P-value_association', 'chisq_ancestry_het', 'ndf_ancestry_het',
                'P-value_ancestry_het', 'chisq_residual_het',
                'ndf_residual_het', 'P-value_residual_het', 'lnBF', 'Comments']

i = 0   # 0 is hdl
pheno = lst_fn_meta_analysis[i].split('_')[0]
print('Processing '+pheno)   # Track what is going on
df_meta = pd.read_csv(dir_meta_analysis + lst_fn_meta_analysis[i], sep='\t', header=None, dtype='str')
df_meta.columns = column_title
# print(df_meta.shape)
# print(df_meta.head())
# print('-----------------------------\n')

# df_single = pd.read_csv(dir_single_analysis_input+lst_fn_single_analysis[i], sep='\t', header=None)
# print(df_single.shape)
# print(df_single.head())
# print('-----------------------------\n')


# Location fo where folders of single analysis are stored
dir_single_analysis_individual_file = '/vgipiper05/DIAMANTE/'
fn = dir_single_analysis_input + lst_fn_single_analysis[i]

# Get list of study names, and file names of each single analysis that went into meta-analysis
lst_study_names, dict_study_dir_fn = get_file_name(fn, dir_single_analysis_individual_file)
print('\n------------------\nBelow individual files went into meta-analysis:')
for k, v in dict_study_dir_fn.items():
    print(k, ':\t', v)
print('\n------------------\n')
"""
df_single_result = pd.read_csv(dict_study_dir_fn[lst_study_names[0]], sep='\t', compression='gzip')

# print(df_single_result.shape, '\n--------------')
# print(df_single_result.head())
df_single_result = df_single_result[['SNP', 'BETA', 'SE']]
# !!!!!!!!???? Choose Beta_0 from meta-analysis results? Double check with Lauren !!!!!!!!????
df_meta = df_meta[['MarkerName', 'beta_0', 'se_0']]
df_merged = df_meta.merge(df_single_result, right_on='SNP', left_on='MarkerName', how='left')
df_merged = df_merged[['MarkerName', 'beta_0', 'se_0', 'BETA', 'SE']]
column_title = ['MarkerName', 'meta_OR', 'meta_SE', 'single_OR', 'single_SE']
df_merged.columns = column_title
print(df_merged.head())
df_merged.to_csv('test.txt', sep='\t', index=False, na_rep='NA')

"""

# !!!!!!!!???? Choose Beta_0 from meta-analysis results? Double check with Lauren !!!!!!!!????
df_meta = df_meta[['MarkerName', 'beta_0', 'se_0']]

for i in range(len(lst_study_names)): # Iterate studyies for file with file names of single analysis results
    study_name = lst_study_names[i]
    # Read in single analysis result
    df_single_result = pd.read_csv(dict_study_dir_fn[study_name], sep='\t', compression='gzip', dtype='str')
    df_merged = df_meta.merge(df_single_result, left_on='MarkerName', right_on='SNP', how='left')
    # merge meta analysis OR, SE with single analysis data OR, SE
    df_merged = df_merged[['MarkerName', 'beta_0', 'se_0', 'BETA', 'SE']]
    column_title = ['MarkerName', 'meta_OR', 'meta_SE', study_name+'_OR', study_name+'_SE']
    df_merged.columns = column_title
    output_fn = '/data100t1/home/wanying/lab_code/utils/output/'+pheno+'_'+study_name+'_merged.txt'
    df_merged.to_csv(output_fn, sep='\t', index=False, na_rep='NA')

    print(pheno+'_'+study_name+' is completed,\nfile is saved at:'+output_fn+'\n---------------')

print('All done!')