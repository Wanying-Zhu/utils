# Functions to process files for CHARGE GWAS
import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np

def winsorization(vals):
    '''
    Winsorize very extreme values:
    If a trait value is more than 6 standard deviations (SD) above or below the mean,
    set it exactly at 6 SDs from the mean
    
    params:
    - vals: array like data to be winsorized
    return:
    - num_winsorized: Number of samples Winsorized
    - vals: array of Winsorized values
    '''
    print('# Winsorize given values')
    vals = np.array(vals)
    avg = np.mean(vals)
    std = np.std(vals)
    upper_threshold, lower_threshold = avg+6*std, avg-6*std
    num_winsorized = (vals>upper_threshold).sum() + (vals<lower_threshold).sum()
    vals = np.where(vals<lower_threshold, lower_threshold, vals)
    vals = np.where(vals>upper_threshold, upper_threshold, vals)
    print('# - Number of values winsorized: %s' % num_winsorized)
    return num_winsorized, vals


def get_residuals(df, exposure_col:str, covar_cols:list):
    '''
    Regress the exposure on required covariates in the sex-combined sample.
    For example, regress out age, sex, and age*sex interaction.
    You may also do this using the R script provided by CHARGE analysis plan.
    
    params:
    - df: dataframe contianing values to be used
    - exposure: a string of column name of the exposure. Eg. depression
    - covar_cols: a list of column names of the covariates to be regressed out. Eg. sex, age, sex*age.
                  Will create interaction column(s) if it does not exist
    return:
    - df:
    '''
    print('# Get residuals of %s' % exposure_col)
    # Sanity checks
    assert exposure_col in df.columns
    for col in covar_cols:
        if '*' not in col:
            assert col in df.columns
        else:
            if col not in df.columns:
                col1, col2 = col.split('*')
                assert col1 in df.columns
                assert col2 in df.columns
                print('# - Create the interaction column: %s*%s' % (col1, col2))
                df[col] = df[col1] * df[col2]
    regr = LinearRegression().fit(df[covar_cols].values, df[exposure_col].values)
    df[exposure_col+'_residuals'] = df[exposure_col].values - regr.predict(df[covar_cols].values)
    percentile_20, percentile_80 = np.quantile(a=df[exposure_col+'_residuals'], q=[0.2, 0.8])
    print('# - 20th percentile:', percentile_20)
    print('# - 80th percentile:', percentile_80)
    return df, percentile_20, percentile_80