Meta-analysis on two input files, sample size weighted

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
--se_cols se \
--beta_cols beta \
--shared_cols protein \
--sample_sizes 300 200 \
--overwrite
```

# Flags
1. ```--input_files```: Input file names of regression results, separate by space. Files must have column headers, columns of pvalue and beta. Provide delimiter(s) or infer from the suffix
2. ```--input_delimiter```: Delimiters of input file. Or provide one value if delimiters are the same in all input files. Use ',', 'tab', 'space' or other values
3. ```--output_path```, ```--output_prefix```: Output path and prefix
4. ```--pval_cols```: Column names of pvalue in each input file, separated by space. Or provide one value if column names are the same in all input files
5. ```--se_cols```, ```--beta_cols```: Column names of standard error and beta in each input file, separated by space. Or provide one value if column names are the same in all input files
6. ```--sample_sizes```: Sample sizes of each regression (integer value). Or provide one value if column names are the same in all input files
7. ```--shared_cols```: Names of shared ID column to merge the regression results, separate by space. Or provide one value if column names are the same in all input files
8. ```--overwrite```: Overwrite existing output file if True. Default value is False


   