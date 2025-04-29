Based on decision of the CCHC group, use the residual PCA method to find hidden covariates.
Similar idea to PEER and SVA, but might have better performance.

# Steps
1. Regress out any covariates (except the phenotype of interest) that would be used in differential abundance analysis.
	* eg. RNAseq ~ sex + age + genetic_PCs
2. Perform PCA on the residuals, save transformed data as PCs
3. Determine the number of PCs to use by the elbow method
	* Similar to the PC selection using a scree plot
4. Use the PCs as covariates in differential abundance analysis
	* eg. RNAseq ~ T2D + sex + age + genetic_PCs + residual_PCs

# Usage
Check example calls at the beginning of the ```./src/get_residual_pca.py``` file.

Get help of options by ```python ./src/get_residual_pca.py --help```

```
# Example
python ./src/get_residual_pca.py \
--input_file /data100t1/home/wanying/CCHC/lipidomics/input_docs/lipidomic_measures/lipid_species_INVed_covar.txt \
--output_path /data100t1/home/wanying/CCHC/lipidomics/20240914_defferential_expresson_in_lipidomics/input/lipid_species_residual_PCA \
--output_prefix lipid_species_INVed_residual_scale \
--covars AGE_AT_VISIT GENDER PC1 PC2 PC3 PC4 PC5 \
--id_col LABID \
--ignore_cols RRID VISIT INTERVIEW_DATE CHOL1 trig hdlc ldlcalc ADA2010_Cat ADA2010_DM BMI1 MED1 MED2 MED3 MED4 MED5 MED6 MED7 MED8 MED9 MED10 genotype_ID MEDS ON_STATIN "MS Label" \
--create_new_covar_file \
--overwrite
```

# Flags
1. ```--input_file```: Input file name, must be in row x col = sample x feature format (eg. rows are samples, columns are gene expressions)
2. ```--output_path, --output_prefix```: Output path and prefix
3. ```--covar_file```: (Optional) Name of the file containing covariates if not in the ```--input_file```, in row x col = sample x feature format. Can be a .csv file or tab-delimited file. (Or a list of samples to use in the input)
4. ```--id_col```: Shared ID column between the input file and covariate file
5. ```--ignore_cols```: Columns to be ignored from the analysis, separated by space
6. ```--covars```: Covariates of the model, such as sex, age, genetic PCs, separated by space
7. ```--scale```: (Default True) Scale (unit variation) the data before PCA. Centering is implemented by PCA already
8. ```--verbose```: (Default False) Print more information.
9. ```--overwrite```: (Default False) Overwrite existing output file if true.
10. ```--total_var_explained```: (Default=0.9) Get number of PCs by this cumulative total variance explained (as a comparison to elbow method)
11. ```--create_new_covar_file```: (Default False) Merge of PCs selected by elbow method with other covariates, and save to a new file for other analysis
12. ```--cols_to_save```: (Optional) If ```--create_new_covar_file``` is True, save the id column, number of PCs by elbow method and these columns to the new file

# Outputs
1. ```prefix.log```: log file
2. ```prefix.model```: model parameters
3. ```prefix.residual```: Residualized input file
4. ```prefix.residual.elbow_*.pca```: (use this file as hidden covaraites, or ```prefix.result_with_covar_cols.residual_pca_with_covar``` if ```--create_new_covar_file``` is True)
	* PCA transfromed residuals, elbow is the number of PCs to use from elbow method
5. ```prefix.residual.pca.scree_plot.png```: PCA scree plot with elbow labeled
6. ```prefix.result_with_covar_cols.residual_pca_with_covar```: when use ```--create_new_covar_file```, save desired columne with number of PCs determined by elbow method

# Reference paper
Zhou, H.J., Li, L., Li, Y. et al. PCA outperforms popular hidden variable inference methods for molecular QTL mapping. Genome Biol 23, 210 (2022). https://doi.org/10.1186/s13059-022-02761-4
