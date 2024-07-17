import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Volcanot plot
def volcano_plot(df, pval_col, fc_col, size_col=None,
                 log_p=True, log_fc=True, title='',
                 logFC_threshold=1, logP_threshold=6,
                 pval_threshold_col='', **kwargs):
    '''
    Create volcano plot
    Params:
    - df: input result dataframe. Must have fold change (or beta) and p value columns.
            Other columns will be ignored
    - pval_col, fc_col: column names of pvalue and fold change
    - log_p, log_fc: will perform -log10(pvalue) or log2(fold change) if set to true
    - title: title of the plot
    - size_col: column name of marker sizes. Usually use gene abundance (mean expression, etc.)
    - logFC_threshold: value of log2FC to draw vertical lines
    - logP_threshold: value of -log10(P) to draw horizontal lines
    - linewidths: line width of the markers
    - pval_threshold_col: will ignore logP_threshold if provided.
                       Column name of adjusted p values to determine significant p (<0.05).
                       Notice that by default adjP provided by DEseq is FDR correted with alpha=0.1
    - **kwargs: other keyword arguments for plt.scatter, such as linewidths, alpha, etc
    Return:
    fig and ax of volcano plot
    
    '''
    df_result = df.copy() # Copy result for data manipulation
    fig, ax = plt.subplots(figsize=(5,5), dpi=150)
    
    # Transform p values
    if log_p:
        df_result[pval_col] = -np.log10(df_result[pval_col])
    
    # Transform fold change
    if log_fc:
        df_result[fc_col] = np.log2(df_result[fc_col])
        
    if size_col:
        ax.scatter(df_result[fc_col], df_result[pval_col], marker='o', edgecolors='w', s=df_result[size_col], **kwargs)
    else:
        ax.scatter(df_result[fc_col], df_result[pval_col], marker='o', edgecolors='w', **kwargs)
    
    # Color genes with small p values
    mask_fc = (df_result[fc_col]>=logFC_threshold) | (df_result[fc_col]<=-logFC_threshold)
    if len(df_result[mask_fc])>0:
        ax.axvline(x=-logFC_threshold, ls='--', color='k', lw=0.5)
        ax.axvline(x=logFC_threshold, ls='--', color='k', lw=0.5)

    if pval_threshold_col != '': # If adjusted p value is provided, color significant genes
        mask_pval = df_result[pval_threshold_col]<=0.05 # By default DESeq provides FDR corrected p values with alpha=0.1 not 0.05
        if len(df_result[mask_pval])>0:
            if size_col:
                ax.scatter(df_result[mask_pval][fc_col], df_result[mask_pval][pval_col], marker='o', edgecolors='w',
                           **kwargs, color='r', s=df_result[mask_pval][size_col])
            else:
                ax.scatter(df_result[mask_pval][fc_col], df_result[mask_pval][pval_col], marker='o', edgecolors='w',
                           **kwargs, color='r')
        logP_threshold = df_result[mask_pval][pval_col].min()
            
    mask_pval = df_result[pval_col]<=logP_threshold
    if len(df_result[mask_pval])>0:
        ax.axhline(y=logP_threshold, ls='--', color='k', lw=0.5)
          
    if title!='': ax.set_title(title)
    ax.set_xlabel('log2(Fold change)')
    ax.set_ylabel('-log10(pvalue)')
    # ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    ax.grid(True, color='whitesmoke')
    ax.set_axisbelow(True) # Put grid at the bottom layer
    fig.tight_layout()
    return fig, ax