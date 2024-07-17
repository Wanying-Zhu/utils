# This code calculate 95% CI, and perform plotting with R
"""
import subprocess

mean = ''
lower = ''
upper = ''
label_text = ''
fig_title = ''
subprocess.Popen(['Rscript','Rcode_forest_plot.R', mean,lower, upper, fig_title])

"""


import pandas as pd
import math
import subprocess
import numpy as np

# ---------------- Helper functions ----------------
# Calculate 95% CI for regression coefficient beta
# beta = ln(OR)
# Returns lower and upper values of 95% confidence intervals
def calculate_95_CI(beta, SE):
    # Formula for OR CI calculation is commented out as below
    # upper = math.e**(math.log(OR) + 1.96 * SE * math.log(OR))
    # lower = math.e**(math.log(OR) - 1.96 * SE * math.log(OR))
    upper = beta + 1.96 * SE * beta
    lower = beta - 1.96 * SE * beta
    return lower, upper
# ---------------- End of helper functions ----------------


dir_input = '/data100t1/home/wanying/lab_code/utils/output/merged_OR_SE/'
lst_fn = ['hdl_OR_SE.txt', 'ldl_OR_SE.txt', 't2d_OR_SE.txt', 'tc_OR_SE.txt', 'tg_OR_SE.txt']

i = 4
fn = lst_fn[i]

df = pd.read_csv(dir_input+fn, sep='\t')
# print(df.head())

# Get SNP rsIDs
lst_rsID = df['MarkerName'].values
# Get study names from column info, this list will be plot label text
# The first on is always meta-analysis
columns = df.columns
lst_study_name = []
for name in columns:
    if name[-2:] == 'OR':
        lst_study_name.append(name.split('_OR')[0])

df.set_index('MarkerName', inplace=True)
# print('\n-------------')
# print(df.head())

# process each variant based on rsID

for j in range(len(lst_rsID)):

    rsID = lst_rsID[j]
    print('Plotting:', rsID)

    lst_mean = []   # A list to store OR
    lst_lower = []  # A list to store 95 confidence interval (lower)
    lst_upper = []
    lst_label_text = lst_study_name.copy()  # Must do deep copy here
    fig_title = rsID
    output_dir = '/data100t1/home/wanying/lab_code/utils/output/forest_plots/'
    fig_name = output_dir + fn.split('_')[0]+'_' + rsID+'.jpeg'

    # Calculate CI for each study, and meta-analysis
    # Notice here OR was really regression coefficient beta!!
    # Need to convert back to OR: Beta = ln(OR)
    # Or just plot beta, and label x-axis as log(OR)
    for study_name in lst_study_name:
        beta = df.loc[rsID, study_name+'_OR'] # Odds ratio
        SE = df.loc[rsID, study_name+'_SE'] # Standard error
        lower_CI, upper_CI = calculate_95_CI(beta, SE)
        lst_mean.append(beta)
        lst_lower.append(lower_CI)
        lst_upper.append(upper_CI)

    # Move meta-analysis result to the last
    lst_mean.append(lst_mean[0])
    lst_mean = lst_mean[1:]
    lst_lower.append(lst_lower[0])
    lst_lower = lst_lower[1:]
    lst_upper.append(lst_upper[0])
    lst_upper = lst_upper[1:]
    lst_label_text.append(lst_label_text[0])
    lst_label_text = lst_label_text[1:]
    lst_label_text[-1] = 'Summary'  # Change the label of meta-analysis to summary

    # Convert lists to the correct argument formant for R script
    # Need to be sth like: 1.1,0.3,0.5 (No space, only separate by ",")
    mean = str(lst_mean[0])
    lower = str(lst_lower[0])
    upper = str(lst_upper[0])
    label_text = lst_label_text[0]
    for i in range(1, len(lst_label_text)):
        if np.isnan(lst_mean[i]): # If the data ids mossing, then don't plot this study
            # do nothing
            pass
        else:
            mean = mean + ',' + str(lst_mean[i])
            lower = lower + ',' + str(lst_lower[i])
            upper = upper + ',' + str(lst_upper[i])
            label_text = label_text + ',' + lst_label_text[i]

    print(len(lst_label_text),'\n', lst_label_text)
    print(len(lst_mean), '\n', mean, '\n=================\n')
    # print(len(lst_lower), '\n', lower,'\n------------')
    # print(len(lst_upper), '\n', upper,'\n------------')
    # print(fig_title,'\n------------')
    # print(fig_name, '\n------------')

    subprocess.Popen(['Rscript','Rcode_forest_plot.R',
                      mean, lower, upper, label_text, fig_title, fig_name])
