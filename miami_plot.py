'''
Example call:
df_upper = df_lipid
df_lower = df_diamante
fig, ax = miami_plot(upper_df=df_upper, lower_df=df_lower,
                     upper_pval_col='P', upper_chr_col='CHR',
                     upper_pos_col='POS', upper_study='Lipidomics',
                     upper_colors=['navy', 'cornflowerblue'],
                     lower_pval_col='P', lower_chr_col='CHR',
                     lower_pos_col='POS_B38', lower_study='Diamante',
                     lower_colors=['tomato', 'lightsalmon'],
                     figsize=(10, 5), dpi=200, markersize=2)
fig.savefig('test.png')
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def miami_plot(upper_df:pd.DataFrame,
               lower_df:pd.DataFrame,
               upper_pval_col:str='P',
               upper_chr_col:str='CHR',
               upper_pos_col:str='POS',
               upper_sig_pval:float=-1,
               upper_study:str='',
               upper_colors:list=['navy', 'cornflowerblue'],
               lower_pval_col:str='P',
               lower_chr_col:str='CHR',
               lower_pos_col:str='POS',
               lower_sig_pval:float=-1,
               lower_study:str='',
               lower_colors:list=['indianred', 'lightcoral'],
               figsize:tuple=(12, 6), dpi=200,
               fig_title='', markersize:float=2):
    '''
    Create Miami plot
    Params:
    - upper_df, lower_df: data to plot in upper and lower figure
    - upper_pval_col, lower_pval_col: name of pvalue column,
    - upper_chr_col, lower_: name of chromsome column. Chromosome numbers can be integers, or chrx, CHRx
    - upper_pos_col, lower_: name of position column,
    - upper_sig_pval, lower_:float=-1,
    - upper_study, lower_study:str='', label to y axis
    - upper_colors, lower_colors: colors to use
    - figsize, markersize
    - fig_title: figure title if not empty
    Return:
    - fig, ax
    - upper_df, lower_df: Data used to plot.
      The "indx" column contains the x values. The "log_pval_upper" column contains -log10(pval)
    '''
    fig, ax = plt.subplots(dpi=dpi, nrows=2, sharex=True, figsize=figsize)
    # To avoid warning: 'A value is trying to be set on a copy of a slice from a DataFrame'
    upper_df = upper_df.copy()
    lower_df = lower_df.copy()
    upper_df['log_pval_upper'] = -np.log10(upper_df[upper_pval_col])
    lower_df['log_pval_lower'] = -np.log10(lower_df[lower_pval_col])

    # Check the first values, remove 'chr' or 'CHR' from chromosome column if needed
    if 'chr' in str(upper_df[upper_chr_col].iloc[0]):
        upper_df[upper_chr_col] = upper_df[upper_chr_col].apply(lambda x: int(x.split('chr')[-1]))
    elif 'CHR' in str(upper_df[upper_chr_col].iloc[0]):
        upper_df[upper_chr_col] = upper_df[upper_chr_col].apply(lambda x: int(x.split('CHR')[-1]))
        
    if 'chr' in str(lower_df[lower_chr_col].iloc[0]):
        lower_df[lower_chr_col] = lower_df[lower_chr_col].apply(lambda x: int(x.split('chr')[-1]))
    elif 'CHR' in str(lower_df[lower_chr_col].iloc[0]):
        lower_df[lower_chr_col] = lower_df[lower_chr_col].apply(lambda x: int(x.split('CHR')[-1]))

    # Get chromosome numbers to plot, include all chromosomes from both association results
    lst_chr = list(set(list(upper_df[upper_chr_col].unique()) + list(lower_df[lower_chr_col].unique())))
    lst_chr.sort()
    shift = 0 # Shift this value when plot current chromosome
    label_pos, label_text = [], [] # position on x axis to label chromosome number, and label text
    font_size = 10 # Font size of axis tick labels
    upper_df['indx'] = 0 # Positions for plotting
    lower_df['indx'] = 0 # Positions for plotting
    for count, chr_num in enumerate(lst_chr): # Plot both upper and lower of a single chromosome at a time
        mask_upper = upper_df[upper_chr_col]==chr_num
        mask_lower = lower_df[lower_chr_col]==chr_num

        # Max and min position of each dataframe to scale (or shift) the positions
        max_upper, min_upper = upper_df[mask_upper][upper_pos_col].max(), upper_df[mask_upper][upper_pos_col].min()
        max_lower, min_lower = lower_df[mask_lower][lower_pos_col].max(), lower_df[mask_lower][lower_pos_col].min()
        max_in_both = [max_upper, max_lower][np.argmax([max_upper, max_lower])]
        min_in_both = [min_upper, min_lower][np.argmin([min_upper, min_lower])]
        # Create values of x axis to plot
        # - Minus min_in_both so position starts at 0
        # - Plus shift so chromosomes do not overlap
        upper_x = upper_df[mask_upper][upper_pos_col] - min_in_both + shift
        lower_x = lower_df[mask_lower][lower_pos_col] - min_in_both + shift
        upper_df.loc[mask_upper, 'indx'] = upper_x
        lower_df.loc[mask_lower, 'indx'] = lower_x
        
        # Positions if x axis ticks
        label_pos.append((upper_x.sum()+lower_x.sum())/(len(upper_x)+len(lower_x))) # Label text in the middle
        label_text.append(chr_num)

        shift += max_in_both - min_in_both
        # print(f'\n# chr{count+1}: max_in_both={max_in_both}, min_in_both={min_in_both}')
        # print(f'# - shift={shift}')
        ax[0].scatter(upper_x, upper_df[mask_upper]['log_pval_upper'], ls='', marker='.',
                      s=markersize, color=upper_colors[count%len(upper_colors)])
        ax[1].scatter(lower_x, lower_df[mask_lower]['log_pval_lower'], ls='', marker='.',
                      s=markersize, color=lower_colors[count%len(lower_colors)])
    
    # Label y axis
    if upper_study != '':
        ax[0].set_ylabel(upper_study+'\n'+'-$log_{10}(p)$', fontsize=font_size+2)
    else:
        ax[0].set_ylabel('-$log_{10}(p)$', fontsize=font_size+2)
        
    if lower_study != '':
        ax[1].set_ylabel(lower_study+'\n'+'-$log_{10}(p)$', fontsize=font_size+2)
    else:
        ax[1].set_ylabel('-$log_{10}(p)$', fontsize=font_size+2) 
        
    # Remove frames
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].set_xticks([label_pos[i] for i in range(0, len(label_pos), 2)])
    ax[1].set_xticklabels([label_text[i] for i in range(0, len(label_text), 2)])
    ax[1].tick_params(top=True, labeltop=True, bottom=False,
                      labelbottom=False, labelsize=font_size)
    ax[0].tick_params(labelsize=font_size) # Set axis font size of both upper and lower figures
    ax[1].invert_yaxis() # Invert y-axis of bottom figure

    # Plot significant line of p value
    if upper_sig_pval == -1: upper_sig_pval = 0.05/len(df_lower)
    ax[0].axhline(y=-np.log10(upper_sig_pval), lw=0.5)
    if lower_sig_pval == -1: lower_sig_pval = 0.05/len(df_lower)
    ax[1].axhline(y=-np.log10(lower_sig_pval), lw=0.5)

    # label figure title
    if fig_title != '':
        fig.suptitle(fig_title)
    fig.subplots_adjust(hspace=0.22) # Adjust space so that chromosome labels are in the middle
    return fig, ax, upper_df, lower_df