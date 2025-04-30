Meta-analysis on two input files, sample size weighted. Need p-value and beta.

Check ```./src/equivalent_metal_script/metal_script.sh``` for equivalent METAL settings.

# Usage
Check example calls at the beginning of the ./src/meta_analysis.py file.

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
--shared_cols protein \
--sample_sizes 300 200 \
--overwrite
```

# Flags
1. ```--input_files```: Input file names of regression results, separate by space. Files must have column headers, columns of pvalue and beta. Provide delimiter(s) or infer from the suffix
2. ```--input_delimiter```: (optional) Delimiters of input file. Or provide one value if delimiters are the same in all input files. Use ',', 'tab', 'space' or other values
3. ```--output_path```, ```--output_prefix```: (optional) Output path and prefix
4. ```--pval_cols```: Column names of pvalue in each input file, separated by space. Or provide one value if column names are the same in all input files
5. ```--beta_cols```: Column names of beta in each input file, separated by space. Or provide one value if column names are the same in all input files
6. ```--sample_sizes```: Sample sizes of each regression (integer value). Or provide one value if column names are the same in all input files
7. ```--shared_cols```: Names of shared ID column to merge the regression results, separated by space. Or provide one value if column names are the same in all input files
8. ```--overwrite```: (optional) Overwrite existing output file if True. Default value is False


# Output
Available columns in ```*.meta_output.txt```:

1. **MarkerName**: marker used to merge the two results. Eg. protein id
2. **Weight**: If the marker is meta-analyzed, weight=sample size1 + sample size2. Otherwise, output the sample size of the corresponding analysis
3. **Zscore**: Meta-analyzed z score, or z score of the corresponding analysis if the marker is only found in one analysis.
4. **P-valu**e: Meta-analyzed p value, or p value of the corresponding analysis if the marker is only found in one analysis.
5. **Direction**: '+' positive beta; '-' negative beta; '?' missing in the analysis
   