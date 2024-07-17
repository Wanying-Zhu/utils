# OLS or logistic regression for GWAS, TWAS or other omic-wide associations
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import os.path

def _progress_bar(current, total):
    '''
    Print progress bar of a loop
    '''
    percent = (current/total) * 100
    bar = '=' * int(percent) + '>' + '-' * int(99 - int(percent))
    print(f'|{bar}| {percent:.2f}%', end='\r')

def _check_input(inputs):
    '''
    Check if input is valid:
    1. a file name (string) ending with .txt or .csv: Read in the file and return a DataFrame. Or
    2. Print error message if input is not a string
    '''
    # Read file into dataframe if not already
    if isinstance(inputs, str):
        if inputs.endswith('.txt'):
            try:
                df = pd.read_csv(inputs, sep='\t')
                return df
            except:
                print(f'#Error: File {inputs} not found\nExit')
        elif inputs.endswith('.csv'):
            try:
                df = pd.read_csv(inputs)
                return df
            except:
                print(f'#Error: File {inputs} not found\nExit')
        else:
            print('#Error: File name must ends with .txt or .csv\Exit')
    else:
        # Print error if input is not file name or dataframe
        print('#Error: Input is not a string or DataFrame\nEixt')

def _check_val_names(var_names, col_names):
    '''
    Check output name and covariate names are valid column names
    - var_name: a list of variable names
    Return False and the invalid variable name
    '''
    for val in var_names:
        if val not in col_names:
            return (False, val)
    return (True, None)

def regressor(inputs,
              outcome:str,
              covariates:list,
              reg_type:str='OLS'):
    '''
    This function performs OLS or logistic regression for GWAS, TWAS or other omic-wide associations
    Param:
    - inputs: input dataframe containing outcome and independent variables
    - outcome: a single column name of dependent variable (y)
    - covariates: a list of names of independent variable(s) (X)
    - reg_type: Default is 'OLS'. OLS (continuous) or logistic regression (binary)
    Return:
    - Fitted regression model
    - 3 lists of results: p-value, beta, SE
    '''
    # ----------------- Sanity checks -----------------
    # Check if variable names exist in column names
    var_name_flag, _ = _check_val_names([outcome] + covariates, inputs.columns)
    if not var_name_flag:
        print(f'#Error: column {_} not found in input dataframe. END')
        exit()
    # ----------------- Regression -----------------
    # Determine regression type
    # y = inputs[outcome]
    # X = inputs[covariates]
    formula = f'{outcome}~' + '+'.join(covariates)
    if reg_type.upper() == 'OLS':
        reg = smf.ols(formula, inputs, missing='drop') # Drop missing values when necessary
        result = reg.fit()
    elif reg_type.upper() == 'LOG':
        reg = smf.logit(formula, inputs, missing='drop')
        result = reg.fit(maxiter=100) # Matching inter number in scikit-learn
    elif reg_type.upper() == 'MIXED':
        # TODO: need to handle formular differently
        reg = smf.mixedlm(formula, inputs, missing='drop')
        # By default in statsmodels source code: method = ['bfgs', 'lbfgs', 'cg']. Add more options in case detault ones failed
        result = reg.fit(method=['powell', 'lbfgs', 'bfgs', 'cg'])
    else:
        print(f'#Error: Unrecognized regression type: {reg_type}')
        return # return None if regression type is not valid
    se = result.bse # Standard error
    betas = result.params # Coefficients
    pvals = result.pvalues # p values
    return result, pvals, betas, se

def regression_loop(inputs,
                    output_fn: str,
                    outcome: str,
                    adj_covariates: list,
                    reg_type: str = 'OLS',
                    p_adj: str='BF'):
    '''
    This function performs OLS or logistic regression for every SNP (can also be protein or other covariates).
    Adjust for given covariates (adj_covariates) in the input data (such as age, sex, PCs, PEER factors, etc.)
    Param:
    - inputs: input file name (use .txt or .csv suffix) or dataframe.
            Each row should contain data of a single subject, each should contain a single variable (accross all subject).

    - outcome: a single column name of dependent variable (y)
    - adj_covariates: a list of independent variable(s) to adjust for. Not include SNP (or protein, RNA expression etc.).
                    Regression is performed iteratively on all other columns, adjust for adj_covariates columns
    - output_fn: output file
    - reg_type: Default is 'OLS'. OLS (continuous) or logistic regression (binary)
    - p_adj: Default is None. Multiple testing adjustment of p values. Bonferroni (BF), FDR (FDR) or None (no adjustment)
    Return:
    - Write result to output file using output_fn
    '''
    # ----------------- Sanity checks -----------------
    print('#Sanity checks')
    # Check if output file is valid
    if output_fn == '':
        print('#Error: Output file name is empty. END')
        exit()
    elif os.path.isfile(output_fn):
        output_fn = output_fn+'.v2'
        print(f'#Output file already exists. Results saved in new output file {output_fn}')

    # Check if multiple testing correction is valid
    if p_adj.upper() == 'BF': # Bonfferoni correction
        p_adj = 'BF'
    elif p_adj.upper() == 'FDR': # FDR correction
        p_adj = 'FDR'
    else: # No correction
        p_adj = 'None'

    # Load file
    print('#Load input file')
    if not isinstance(inputs, pd.DataFrame):  # Read in file if input is not a DataFrame yet
        inputs = _check_input(inputs)

    # Check if variable names exist in column names
    var_name_flag, _ = _check_val_names([outcome] + adj_covariates, inputs.columns)
    if not var_name_flag:
        print(f'#Error: column {_} not found in input dataframe. END')
        exit()

    print(f'\n#Passed sanity checks, starting regression')
    print(f'# - Regression type: {reg_type}')
    print(f'# - Outcome: {outcome}')
    print(f'# - Adjust for covariates: {", ".join(adj_covariates)}')
    print(f'# - Multiple testing correction: {p_adj}')
    print(f'# - Results saved in output file {output_fn}\n')
    # ----------------- Done sanity checks -----------------

    # Loop through all columns except columns in adj_covariates
    # Adjust for covaraites in list adj_covariates
    output_fh = open(output_fn, 'w')
    # Write below results to output file:
    # Variable name, coefficient, raw p value, standard error, number of observations
    output_line = 'variable\tbeta\tpval\tse\tnumber_observations' # header line
    output_fh.write(output_line)
    count = 0
    # Total number of variables to loop through
    total = len(inputs.columns) - len(adj_covariates) - 1
    for col in inputs.columns:
        if (col not in adj_covariates) and (col!=outcome):
            cols_for_regression = [outcome, col] + adj_covariates
            result, pvals, betas, se = regressor(inputs=inputs[cols_for_regression],
                                                 outcome=outcome,
                                                 covariates=[col] + adj_covariates,
                                                 reg_type=reg_type)
            output_line = f'\n{col}\t{betas[col]}\t{pvals[col]}\t{se[col]}\t{result.nobs}'
            output_fh.write(output_line)
            count += 1
            if count%100 == 0:
                _progress_bar(count, total)
    _progress_bar(1, 1) # All done
    output_fh.close()

    # Adjust for multiple testing
    if p_adj != 'None':
        print(f'\n\n#Start multiple testing correction: {p_adj}')
        df_result = pd.read_csv(output_fn, sep='\t') # Load result file
        if p_adj == 'FDR':
            df_result[f'pval_{p_adj}_adj'] = fdrcorrection(df_result['pval'])[1]
        elif p_adj == 'BF':
            df_result[f'pval_{p_adj}_adj'] = df_result['pval'] * df_result['number_observations']
    df_result.to_csv(output_fn, sep='\t', index=False, na_rep='NA')
    print('\n#DONE')