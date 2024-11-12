## 1. Scripts available
    1.1 ```linear_regression_with_multiprocessing.py``` and ```regression.py``` are equivalent when doing OLS linear regression. ```linear_regression_with_multiprocessing.py``` does have multiprocessing option (```--threads```), but it may use a lot of resources and take longer on our server. Use a single thread is fine in most cases
    1.2 ```regression.py``` also include linear mixed model option, can be specified by ```--model``` and ```--group_col```
## 2. Example to use the code
    2.1 linear_regression_with_multiprocessing.py
    2.2 regression.py