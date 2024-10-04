# Updated 2024/10/03


# Notes: changed calculation of expected p values according to Anna's code

from scipy.stats import chi2
from scipy.stats import beta
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def __read_file(filename):
    '''
    Read a tabular file containing SNPs and P values to be plotted (Check delimiters automatically)
    Parameter:
    - filename: delimiters are tab, comma or whitespace
    Return:
    - a pandas dataframe of the loaded file
    '''
    with open(filename, 'r') as f:
        line = fh.readline()
    if line=='': return None # Empty file, return None
    # Try to split by ','
    tmp = line.strip().split(',')

    # Create an empty dataframe, read data into df later
    df = pd.DataFrame()
    # Check delimiter type, "," or " "
    # One of the read_csv() call should work
    if len(tmp) == 1:
        try: df = pd.read_csv(filename, sep='\t')
        except pd.errors.ParserError: pass 
        
        try: df = pd.read_csv(filename, sep='\s+') # Split by whitespace
        except pd.errors.ParserError: pass
        # File format maybe wrong if both read_csv failed
    else:
        df = pd.read_csv(filename)
    # File only has a line of tile, return None  
    if len(df)==0: return None 
    return df


def __log_obsv_and_expt(df, column_title):
    '''
    Return -log10 transformed observed and expected p vlaues of a dataset (dataframe)
    '''
    scale = 0.5/np.log(10)
    Pobsv = df.loc[:, column_title]  # Observed p values
    Pobsv = np.sort(Pobsv.dropna().values)
    
    # Calculate expected p
    # A method from qqman package. Refer to ppooints function in R. Saved for future reference
    # The results are slightly differetn from Anna's
    # if len(Pobsv)>10:
    #     Pexpt = (np.arange(1, len(Pobsv)+1)-0.5) / (len(Pobsv)+1-1)
    # else:
    #     Pexpt = (np.arange(1, len(Pobsv)+1)-3.0/8) / (len(Pobsv)+1-2*3.0/8)
        
    #     From Anna's code, she used this one to calculate chi-squared test statistics, but not for logP qq plot
    #     Pexpt = chi2.ppf(np.arange(1, len(Pobsv) + 1, 1) / (len(Pobsv) + 1), df=1)

    # Use Pobsv[Pobsv!=0] since log will not work on zeros (avoid warning message)
    ln_Pexpt = 2*np.cumsum(1/np.arange(len(Pobsv[Pobsv!=0]), 0, -1)) # -2*ln(Pexpt), from Anna's code, don't understand why
    logPobsv = -np.log10(Pobsv[Pobsv!=0])
    logPexpt = ln_Pexpt * scale # Convert from -2ln(x) to -log(x). from Anna's code
    return Pobsv, np.sort(logPobsv), logPexpt


def __get_inflation__(observed):
    '''
    Calculate and return lambda (inflation)
    This is different from Anna's code, which used mean(obsvd)/mean(expctd)?
    '''
    obsv_median = np.median(observed)
    Chi = chi2.ppf(1.0 - obsv_median, 1) # 1 refers to freedom of 1
    lmbd = Chi / chi2.ppf(0.5, 1)
    return lmbd


def __plot_ci(ax, ci, logPexpt):
    '''
    Plot confidence interval
    '''
    n = len(logPexpt)
    df = 2 # Degree of freedom, from Anna's code (Why 2?)
    # From Anna's code. Scale is used to convert 2ln(x) to log10(x), but I don't understand it
    scale = 0.5/np.log(10)

    a = list(range(1, n+1))
    b = list(range(n, 0, -1))

    clower = beta.ppf(q=(1-ci)/2, a = a, b = b)
    cupper = beta.ppf(q=1-(1-ci)/2, a = a, b = b)
    lower_bond = chi2.ppf(q=clower, df=df)
    upper_bond = chi2.ppf(q=cupper, df=df)

    ax.fill_between(x=logPexpt, y1=lower_bond*scale, y2=upper_bond*scale,
                     alpha=0.2, color='k', linewidth=0, zorder=0) # Plot at the bottom

def qqplot(data, output='output.png', p_value_col = 'P', plot_ci=True, title='Q-Q plot',
           xlabel='Expected –log10 P-values', ylabel='Observed –log10 P-values', dpi=300,
           selected_SNPs_col = None, ci=0.95, fig=None, ax=None, savefig=True, **kwargs):
    '''
    Parameters:
    - data: name of input file that contains all SNPs discovered, should have column headers.
            Can also be a dataframe, then will skip reading in the file       
    - output: name of the output file to save figure of the QQ plot
    - p_value_col: header of the column that contains P values in the input file
    - title, xlabel, ylabel, dpi=300: label and resolution parameters of the figure
    - selected_SNPs_col: Column title with labels for the selected SNPs (1=selected, 0=not).
                         Will calculate expected pvals of selected SNPs only  and plot
    - plot_ci: plot shaded area of confidence interval if true
    - ci: confidence interval, default is 0.95
    - fig=None, ax=None: figure and ax to plot if provided
    - savefig=True: if figure needs to be saved
    
    Returns: (fig, ax, lambda_original, lambda_novel)
    - fig and ax for more customizations
    
    '''
    if not isinstance(data, pd.DataFrame): # Read in file if data is not a dataframe yet
        # Read in file, calculate -log10 p values and lambda
        df = __read_file(data)
    else: df = data
    
    if df is None:
        print('Empty file, nothing is plotted')
        return(None, None, None) # If file is empty or only has titles, return directly

    if not fig or not ax:
        fig, ax = plt.subplots(dpi=dpi)
    
    # Plot all SNPs as dots
    Pobsv, logPobsv, logPexpt = __log_obsv_and_expt(df, p_value_col)
    infl_original = __get_inflation__(Pobsv)
    
    ax.plot(logPexpt, logPexpt, color='r', linewidth=0.4)
    ax.plot(logPexpt, logPobsv, linestyle='', marker='o', markersize=2, markeredgewidth=0.5,
            fillstyle='none', color='k', zorder=3, **kwargs)
    
    # Plot shaded area of confidence intervel
    if plot_ci==True:
        if ci >= 1: print("Confidence intervel >= 1, no ci plotted")
        elif ci <= 0: print("Confidence intervel <= 0, no ci plotted")
        else: __plot_ci(ax, ci, logPexpt)
    
    infl_novel = ''
    # Plot selected SNPs if column name is provided
    if not selected_SNPs_col:
        print('No selected SNPs provided\nOnly plot all SNPs')
    else:
        # Plot selected SNPs, the same steps as before
        Pobsv_novel, logPobsv_novel, logPexpt_novel = __log_obsv_and_expt(df[df[selected_SNPs_col]==1], p_value_col)
        infl_novel = __get_inflation__(Pobsv_novel)
        
        ax.plot(logPexpt_novel, logPobsv_novel, linestyle='', marker='o', markersize=2, markeredgewidth=0.5,
                fillstyle='none', color='g', zorder=1)
        
        # annotation = "λ (novel) = " + str("{0:.4f}".format(infl_novel))
        # ax.annotate(annotation, xy=(0.7, 0.1), xycoords='axes fraction')


    # Label x and y axis
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    annotation = "λ = " + str("{0:.4f}".format(infl_original))
    ax.annotate(annotation, xy=(0.7, 0.2), xycoords='axes fraction')

    if savefig: fig.savefig(output)
    
    return fig, ax, infl_original, infl_novel  # Return for more custermizations