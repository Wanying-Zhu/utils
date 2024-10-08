'''
Example call
python /data100t1/home/wanying/CCHC/lipidomics/20240914_defferential_expresson_in_lipidomics/code/utils/get_residual_pca.py \
--input_file /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_measures/lipid_species_INVed_covar.txt \
--output_path /data100t1/home/wanying/CCHC/lipidomics/20240914_defferential_expresson_in_lipidomics/input/lipid_species_residual_PCA \
--output_prefix lipid_species_INVed_residual_scale \
--covars AGE_AT_VISIT GENDER PC1 PC2 PC3 PC4 PC5 \
--id_col LABID \
--ignore_cols RRID VISIT INTERVIEW_DATE CHOL1 trig hdlc ldlcalc ADA2010_Cat ADA2010_DM BMI1 MED1 MED2 MED3 MED4 MED5 MED6 MED7 MED8 MED9 MED10 genotype_ID MEDS ON_STATIN "MS Label" \
--overwrite \
--threads 8

python /data100t1/home/wanying/CCHC/lipidomics/20240914_defferential_expresson_in_lipidomics/code/utils/get_residual_pca.py \
--input_file /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_measures/lipid_species_INVed_covar.txt \
--output_path /data100t1/home/wanying/CCHC/lipidomics/20240914_defferential_expresson_in_lipidomics/input/lipid_species_residual_PCA \
--output_prefix test_result \
--covars AGE_AT_VISIT GENDER PC1 PC2 PC3 PC4 PC5 \
--covar_file xxx \
--id_col LABID \
--ignore_cols RRID VISIT INTERVIEW_DATE CHOL1 trig hdlc ldlcalc ADA2010_Cat ADA2010_DM BMI1 MED1 MED2 MED3 MED4 MED5 MED6 MED7 MED8 MED9 MED10 genotype_ID MEDS ON_STATIN "MS Label" \
--overwrite \
--threads 8 \
--phenotype

'''
import pandas as pd
import numpy as np
import argparse
import logging
import sys
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler # Scale the data
import matplotlib.pyplot as plt

from _setup_logger_and_parse_args import process_args
from _lr_residualization_with_multiprocessing import residualization
from _elbow_method_find_number_of_pcs import get_n_pcs_by_elbow


# ########## Get arguments and residualization ##########
args = process_args()
residualization(args)

# PCA on residuals
logging.info('# Run PCA on residuals (PCA transformed residuals)')
fn_residual = f'{args.output_path}/{args.output_prefix}.residual'
df_residual = pd.read_csv(fn_residual, sep='\t')
X = df_residual.iloc[:, 1:]  # The first column is ID column
if args.scale: # Scale the data
    logging.info('# - Scale the residuals')
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

logging.info('# - Perform PCA')
pca = PCA().fit(X)
df_transformed = pd.DataFrame(pca.transform(X)) # Get transformed data as PCs
df_transformed.columns = [f'PC{x+1}' for x in range(df_transformed.shape[-1])]
df_transformed = pd.concat([df_residual[args.id_col], df_transformed], axis=1)

var_explained_col = 'explained_variance_ratio'
df_pca_var_explained = pd.DataFrame(data = pca.explained_variance_ratio_, columns=[var_explained_col])
df_pca_var_explained['PC'] = [x+1 for x in range(len(df_pca_var_explained))]
df_pca_var_explained.to_csv(f'{args.output_path}/{args.output_prefix}.explained_variance_ratio', sep='\t', index=False)

# Find number of PCs by elbow method
k, _ = get_n_pcs_by_elbow(df_pca_var_explained=df_pca_var_explained, explained_variance_col=var_explained_col)
logging.info('# - Number of PCs to use by elbow method: %s' % k)

# Find number of PCs of x% variance explained in total
n_pcs_by_variance = np.nonzero(np.cumsum(df_pca_var_explained[var_explained_col].values)>args.total_var_explained)[0][0]
logging.info('# - Number of PCs to use to capture '+str(args.total_var_explained*100)+'% of the variance: '+str(n_pcs_by_variance))

logging.info('# - Save PCs')
out_fn = f'{args.output_path}/{args.output_prefix}.residual.elbow_{k}.pca'
df_transformed.to_csv(out_fn, sep='\t', index=False)

# Plot scree plot
plot_fn = f'{args.output_path}/{args.output_prefix}.residual.pca.scree_plot.png'
    
# df_pca_var_explained
fig, ax = plt.subplots()
ax.plot(df_pca_var_explained['PC'],df_pca_var_explained[var_explained_col], ls='--', marker='.')
plt.axvline(x=k, color='r', label=f'Elbow method: N={k}', lw=1)
plt.axvline(x=n_pcs_by_variance, color='b', lw=1, label=f'To capture {args.total_var_explained*100}% of the variance: N={n_pcs_by_variance}')
ax.legend()
ax.set_xlabel('PC')
ax.set_ylabel(var_explained_col)
ax.set_title(f'Elbow = PC{k}')
fig.savefig(plot_fn)

logging.info('# DONE')

