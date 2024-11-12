## 1. Scripts available
1.1 ```linear_regression_with_multiprocessing.py``` and ```regression.py``` are equivalent when doing OLS linear regression. ```linear_regression_with_multiprocessing.py``` does have multiprocessing option (```--threads```), but it may use a lot of resources and take longer on our server. Use a single thread is fine in most cases
1.2 ```regression.py``` also include linear mixed model option, can be specified by ```--model``` and ```--group_col```

## 2. Example to use the code
2.1 linear_regression_with_multiprocessing.py
```
python linear_regression_with_multiprocessing.py \
--input_file example_data/input.csv \
--output_path output \
--output_fn result_multiprocessing.txt \
--covars AGE GENDER PC1 PC2 PC3 \
--plot \
--condition CONDITION1 \
--ignore_cols ID CONDITION2 \
--id_col ID \
--threads 8 \
--get_residual \
--overwrite \
--permute 100
```

2.2 regression.py
```
python regression.py \
--input_file example_data/input.csv \
--output_path output \
--output_fn result_v2.txt \
--covars AGE GENDER PC1 PC2 PC3 \
--condition CONDITION1 \
--model ols \
--ignore_cols ID CONDITION2 \
--id_col ID \
--get_residual \
--plot \
--overwrite
```