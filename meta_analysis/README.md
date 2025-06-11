Meta-analysis on two input files, sample size weighted. Need p-value and beta.

Check ```./src/equivalent_metal_script/metal_script.sh``` for equivalent METAL settings.

# Usage
Check example calls at the beginning of the ./src/meta_analysis.py file.

In CCHC proteomics: ```batch1=3072```, ```batch2=HT```

Get help for options by ```python ./src/meta_analysis.py --help```

```
# Example to run meta-analysis in terminal
result1=./example_data/result1.txt
result2=./example_data/result2.txt
python ./src/meta_analysis.py \
--input_files ${result1} ${result2} \
--output_path ./example_output \
--output_prefix python_meta_output \
--input_delimiter tab \
--pval_cols pval \
--beta_cols beta \
--marker_cols protein \
--sample_sizes 300 200 \
--overwrite

# Example in kids project
# Map marker IDs in reuslts 1 and 2 via file --id_mapping_fn, also save extra columns x, y, z in the output
# The ID mapping file is provided with the code in ./data/cchc_proteomics_bathc1_2_metaID_5511proteins.txt
result1=./example_data/result1.extra_cols.txt
result2=./example_data/result2.extra_cols.txt
python ./src/meta_analysis.py \
--input_files ${result1} ${result2} \
--output_path ./example_output \
--output_prefix python_meta_output.extra_cols \
--id_mapping_fn ./data/cchc_proteomics_bathc1_2_metaID_5511proteins.txt \
--id_mapping_cols OID3072 OIDHT meta_id \
--extra_cols_to_keep meta_id OID3072 OIDHT protein \
--input_delimiter tab \
--pval_cols pval \
--beta_cols beta \
--marker_cols OID3072 OIDHT \
--sample_sizes 300 200 \
--overwrite

```

# Flags
1. ```--input_files```: Input file names of regression results, separate by space. Files must have column headers, columns of pvalue and beta. Provide delimiter(s) or infer from the suffix
2. ```--input_delimiter```: (optional) Delimiters of input file. Or provide one value if delimiters are the same in all input files. Use ',', 'tab', 'space' or other values
3. ```--output_path```, ```--output_prefix```: (optional) Output path and prefix
4. ```--id_mapping_fn```: (Optional) Needed when marker IDs do not match in results 1 and 2. A file with an ID mapping scheme. Expects 3 columns (Column names of this file are provided by ```--id_mapping_cols```):
    * Marker ID column in result 1
    * Marker ID column in result 2
    * A shared marker IDs.
6. ```--id_mapping_cols```: (Optional) Column names in file --id_mapping_fn. The 1st, 2nd and 3rd values refer to column names corresponding to these columns in --id_mapping_fn:
    * ID column map to result 1
    * ID column map to result 2
    * shared ID column 
    If omitted, then assume the ID columns are the same as ```--marker_cols```, and use meata_id for the shared ID column 
7. ```--pval_cols```: Column names of pvalue in each input file, separated by space. Or provide one value if column names are the same in all input files
8. ```--beta_cols```: Column names of beta in each input file, separated by space. Or provide one value if column names are the same in all input files
9. ```--sample_sizes```: Sample sizes of each regression (integer value). Or provide one value if column names are the same in all input files
10. ```--marker_cols```: Names of shared ID column to merge the regression results, separated by space. Or provide one value if column names are the same in all input files
11. ```--extra_cols_to_keep```: Extra columns in the individual result (and ID mapping file if provided) to keep in the meta output
12. ```--overwrite```: (optional) Overwrite existing output file if True. Default value is False


# Output
Available columns in ```*.meta_output.txt```:

1. **MarkerName**: marker used to merge the two results. Eg. protein id
2. **Weight**: If the marker is meta-analyzed, weight=sample size1 + sample size2. Otherwise, output the sample size of the corresponding analysis
3. **Zscore**: Meta-analyzed z score, or z score of the corresponding analysis if the marker is only found in one analysis.
4. **P-valu**e: Meta-analyzed p value, or p value of the corresponding analysis if the marker is only found in one analysis.
5. **Direction**: '+' positive beta; '-' negative beta; '?' missing in the analysis
   