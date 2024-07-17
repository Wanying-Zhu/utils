import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Manhattan plot
def manhattan_plot(data: pd.core.frame.DataFrame, pval: str='pval', position: str='pos',
                   chromosome: str='chr', gene: str='', title: str='manhattan plot',
                   sig_pval: float=-1, annotate: bool=False, colors: list=['black','grey'],
                   markersize: float=2, dpi: float=200, figsize: tuple=(8,4)):
    '''
    Plot Manhattan plot from summary statistics
    Params
    - data: A dataframe containing summary statistics
    - pval: column name of p values
    - position: column name of position
    - chromosome: column name of chromosome
    - gene: column of gene names for annotation (only used if annotate=True)
    - title: figure title
    - sig_pval: threshold of significant p value (No multiple testing correction will be applied on this threshold)
                Default values -1 implies 0.05/number_of_tests is used.
    - annotate: anotate significant points (BF correction)
    - colors: list of colors to plot each chromosome. Default is ['black','grey']
    - markersize: default is 2. Size of markers
    - dpi, figsize: resolution and figure size
    Return
    - fig, ax
    '''
    if annotate: # If name of gene column is provided
        data_copy = data[[gene, chromosome, position, pval]].sort_values(by=[chromosome, position]).copy()
    else:
        data_copy = data[[chromosome, position, pval]].sort_values(by=[chromosome, position]).copy()
        
    data_copy['indx'] = [x for x in range(len(data_copy))] # Create arbitrary indices for plotting
    data_copy['log_pval'] = -np.log10(data_copy[pval]) # Plot -log10 pvalues
    label_text, label_pos = [], [] # Keep track of label text and label positions
    grouped = data_copy.groupby(by=chromosome)
    
    fig, ax = plt.subplots(dpi=dpi, figsize=figsize)
    count=0
    for chr_num, df in grouped: # Plot each group (grouped by chromosome)
        label_text.append(chr_num)
        label_pos.append(df['indx'].mean())
        ax.plot(df['indx'], df['log_pval'], ls='', marker='.',
                markersize=markersize, color=colors[count%len(colors)])
        count += 1
        
    if title is not None: ax.set_title(title)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-$log_{10}(p)$')
    
    # Plot significant line of p value
    if sig_pval == -1: sig_pval = 0.05/len(data_copy)
    ax.axhline(y=-np.log10(sig_pval), lw=0.5)
    ax.set_xticks([label_pos[i] for i in range(0, len(label_pos), 2)])
    ax.set_xticklabels([label_text[i] for i in range(0, len(label_text), 2)], fontsize='8')
    
    if annotate:
        if gene=='': raise ValueError('name of gene column cannot be empty')
        sig_points = data_copy[data_copy[pval]<=(sig_pval)] # Significant data points
        for i in range(len(sig_points)):
            tmp = sig_points.iloc[[i],:] # Cut the temp data for plotting
            ax.annotate(text=tmp[gene].values[0],
                        xy=(tmp['indx'].values[0], tmp['log_pval'].values[0]), fontsize='5')
    
    return fig, ax
