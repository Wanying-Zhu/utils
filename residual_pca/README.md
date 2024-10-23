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
--overwrite \
--threads 8
```

# Outputs
1. ```prefix.log```: log file
2. ```prefix.model```: model parameters
3. ```prefix.residual```: Residualized input file
4. ```prefix.residual.elbow_*.pca```: (use this file as hidden covaraites)
	* PCA transfromed residuals, elbow is the number of PCs to use from elbow method
5. ```prefix.residual.pca.scree_plot.png```: PCA scree plot with elbow labeled

# Reference paper
Zhou, H.J., Li, L., Li, Y. et al. PCA outperforms popular hidden variable inference methods for molecular QTL mapping. Genome Biol 23, 210 (2022). https://doi.org/10.1186/s13059-022-02761-4 
