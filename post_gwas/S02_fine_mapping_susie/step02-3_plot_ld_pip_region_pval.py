# Plot some figures:
# 1. regionnal GWAS pvalues (regional Manhattan Plot)
# 2. Fine-mapping PIP
# 3. Regional LD

'''
Example call: TBA

python /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/code/utils/step02-3_plot_ld_pip_region_pval.py \
--output_path output \
--output_prefix test \
--ld_matrix /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/ld_matrix/snps_filtered_by_pval0.05/PC-18:2_18:2-.chr11.region_1.ld \
--gwas_fine_map_result /data100t1/home/wanying/CCHC/lipidomics/post_gwas_finemapping/output/susieR_fine_mapping/output/PC-18:2_18:2-.chr11.region_1.susie.txt  \
--pval_threshold 6.85e-10 \
--colname_pip PIP \
--colname_pval P
'''
import pandas as pd
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os
import logging
import sys
sys.path.append('/data100t1/home/wanying/lab_code/utils')
from setup_logger_and_parse_args import process_args


def heatmap(df_ld_matrix, ax1, ax2):
    '''
    Plot LD matrix r and r2 on ax1 and ax2
    Params:
    - df_data: LD matrix from plink --r square
    - ax1, ax2
    '''
    sns.heatmap(df_ld_matrix, ax=ax1, cmap='RdYlBu_r', square=True) # Reverse color
    ax1.set_title('LD matrix r')
    ax1.set_xlabel('SNP')
    ax1.set_ylabel('SNP')
    
    sns.heatmap(df_ld_matrix**2, ax=ax2, cmap='YlOrRd', square=True)
    ax2.set_title('LD matrix r2')
    ax2.set_xlabel('SNP')

def reginal_manhattan_plot(df_data:pd.DataFrame, ax,
                           lead_snp_pos_mb, lead_snp_p,
                           colors=['blue', 'lightblue', 'limegreen', 'orange', 'red'],
                           sig_pval:float=-1,
                           pval_col:str='P',
                           no_ld:bool=False):
    '''
    Manhattan plot of a given GWAS region
    Params:
    - df_gwas: GWAS summary stats, with p value
    - colors: colors of the dots by LD r2
    - lead_snp_pos_mb, lead_snp_p: position in Mb and p value of the lead SNP
    - sig_pval: threshold of significant p value
    - chr_col, position_col, pval_col: column names of chromosome, position, p values
    - ax: ax to plot
    - no_ld: if LD is not provided (no_ld=True), will not color dots by LD with the lead SNP
    
    '''
    alpha = 1
    size=30
    edge_color = 'k'
    line_width = 0.5
    # colors = ['blue', 'lightblue', 'limegreen', 'orange', 'red']
    # colors = ['darkgrey', 'silver', 'lightgrey', 'whitesmoke' , 'white']
    # colors = ['grey', 'aliceblue', 'azure', 'antiquewhite', 'lightsalmon']
    # colors = ['lightgrey']*5
    if not no_ld: # Color by LD
        mask = df_data['LD']<0.2
        ax.scatter(df_data[mask]['POS_mb'], -np.log10(df_data[mask][pval_col]),
                    c=colors[0], lw=line_width, edgecolors=edge_color, alpha=alpha, s=size, label='0-0.2')
        mask = (df_data['LD']>=0.2) & (df_data['LD']<0.4) 
        ax.scatter(df_data[mask]['POS_mb'], -np.log10(df_data[mask][pval_col]),
                    c=colors[1], lw=line_width, edgecolors=edge_color, alpha=alpha, s=size, label='0.2-0.4')
        mask = (df_data['LD']>=0.4) & (df_data['LD']<0.6) 
        ax.scatter(df_data[mask]['POS_mb'], -np.log10(df_data[mask][pval_col]),
                    c=colors[2], lw=line_width, edgecolors=edge_color, alpha=alpha, s=size, label='0.4-0.6')
        mask = (df_data['LD']>=0.6) & (df_data['LD']<0.8) 
        ax.scatter(df_data[mask]['POS_mb'], -np.log10(df_data[mask][pval_col]),
                    c=colors[3], lw=line_width, edgecolors=edge_color, alpha=alpha, s=size, label='0.6-0.8')
        mask = df_data['LD']>=0.8
        ax.scatter(df_data[mask]['POS_mb'], -np.log10(df_data[mask][pval_col]),
                    c=colors[4], lw=line_width, edgecolors=edge_color, alpha=alpha, s=size, label='0.8-1')
        ax.legend(title='LD $r^2$', frameon=False)
    else: # Do not color dots by LD
        ax.scatter(df_data['POS_mb'], -np.log10(df_data[pval_col]),
                   c='lightgrey', lw=line_width, edgecolors=edge_color, alpha=alpha, s=size)
    # Plot lead SNP
    ax.scatter(lead_snp_pos_mb, -np.log10(lead_snp_p),
                marker='D', c='blueviolet', lw=line_width, edgecolors='k', s=size+10)

    # Plot significant line of pval
    if sig_pval==-1:
        # Default 1e-8
        ax.axhline(-np.log10(1e-8), c='grey', ls='--', alpha=0.5)
    else:
        ax.axhline(-np.log10(sig_pval), c='grey', ls='--', alpha=0.5)
    # ax.set_title(f'{title}')
    ax.set_ylabel('-log10(p value)')
    # ax.label_outer()

def plot_pip(df_data:pd.DataFrame, ax,
             lead_snp_pos_mb, lead_snp_p,
             sig_pval:float=-1, pip_col:str='PIP'):
    '''
    Plot PIP
    Params:
    - df_gwas: GWAS summary stats, with p value
    - colors: colors of the dots by LD r2
    - lead_snp_pos_mb, lead_snp_p: position in Mb and p value of the lead SNP
    - sig_pval: threshold of significant p value
    - pip_col: column name of pip
    - ax: ax to plot
    
    '''
    alpha = 1
    size=30
    edge_color = 'k'
    line_width = 0.5
    ax.scatter(df_data['POS_mb'], df_data[pip_col],
               c='lightgrey', lw=line_width, edgecolors=edge_color, alpha=alpha, s=size)
    
    # Plot significant line of pval
    if sig_pval==-1:
        pass # Do not plot
        # # Default 1e-8
        # ax.axhline(-np.log10(1e-8), c='grey', ls='--', alpha=0.5)
    else:
        ax.axhline(-np.log10(sig_pval), c='grey', ls='--', alpha=0.5)
    ax.set_ylabel('PIP')


if __name__ == "__main__":
    # Already has --output_path and --output_prefix
    args = process_args(True,
                        {'flag_name':'--gwas_fine_map_result', 'type': 'str',
                         'help': 'File name of the fine map and GWAS result to be processed'},
                        {'flag_name':'--ld_matrix', 'type': 'str',
                         'help': 'File name of the plink LD matrix used for fine mapping'},
                        {'flag_name':'--pval_threshold', 'default': 5e-8, 'type': 'float',
                         'help': 'P value threshold to filter when looking for a lead SNP'},
                        {'flag_name':'--colname_pos', 'default': 'POS', 'type': 'str',
                         'help': 'Column name of the positions in the GWAS result'},
                        {'flag_name':'--colname_pip', 'default': 'PIP', 'type': 'str',
                         'help': 'Column name of the PIP from fine mapping'},
                        {'flag_name':'--colname_pval', 'default': 'P', 'type': 'str',
                         'help': 'Column name of the p value in the GWAS result'})
    output_fn = os.path.join(args.output_path, args.output_prefix+'.png')
    logging.getLogger('matplotlib.font_manager').disabled = True # Disable matplotlib font message
    
    logging.info('\n# Load fine-mapping and GWAS result')
    # Load LD matrix
    # Load GWAS summary stats
    # Load fine-mapping output
    df_gwas_fine_map = pd.read_csv(args.gwas_fine_map_result, sep='\t')
    df_gwas_fine_map['indx'] = [x for x in range(len(df_gwas_fine_map))] # Create arbitrary indices for plotting
    logging.info('# - shape (%s, %s)' % df_gwas_fine_map.shape)
    # Get MB position for plotting
    df_gwas_fine_map['POS_mb'] = df_gwas_fine_map[args.colname_pos]/1000000
    logging.info('# Load LD matrix')
    df_ld_matrix = pd.read_csv(args.ld_matrix, sep='\t', header=None)
    logging.info('# - shape (%s, %s)' % df_ld_matrix.shape)
    
    # Merge (concat) the two DataFrames
    # Only take the column of lead SNP
    df_lead_snp = df_gwas_fine_map[df_gwas_fine_map['lead_snp']==1]
    index_lead_snp = df_lead_snp.index[0]
    lead_snp_p = df_lead_snp[args.colname_pval].values[0]
    lead_snp_pos_mb = df_lead_snp['POS_mb'].values[0]
    lead_pip = df_lead_snp[args.colname_pip].values[0]
    
    df_ld = df_ld_matrix.iloc[:, [index_lead_snp]].copy()
    df_ld.columns=['LD']
    df_ld['LD'] = df_ld['LD']**2 # Need r2
    df_merged = pd.concat([df_gwas_fine_map, df_ld], axis=1)
    chromosome = df_merged['CHR'].iloc[0]

    
    fig = plt.figure(figsize=(12,16),dpi=150)
    
    # Plot heatmap of LD matrix
    ax0 = plt.subplot(421)
    ax1 = plt.subplot(422, sharey=ax0)
    heatmap(df_ld_matrix, ax0, ax1)
    
    # ########## Regional plot ##########
    ax2 = plt.subplot(412)
    reginal_manhattan_plot(df_data=df_merged, ax=ax2,
                           lead_snp_pos_mb=lead_snp_pos_mb, lead_snp_p=lead_snp_p,
                           sig_pval=6.85e-10, pval_col=args.colname_pval, no_ld=False)
    
    # ########## Regional plot with credible sets ##########
    ax3 = plt.subplot(413, sharex=ax2)
    reginal_manhattan_plot(df_data=df_merged, ax=ax3,
                           lead_snp_pos_mb=lead_snp_pos_mb, lead_snp_p=lead_snp_p,
                           sig_pval=6.85e-10, pval_col=args.colname_pval, no_ld=True)
    
    line_width, size = 0.5, 30
    # Plot credible set 
    colors = ['salmon', 'yellow', 'green', 'cyan', 'royalblue']*10
    c = 0
    for col in df_merged.columns:
        if '95_credible_set' in col:
            mask = df_merged[col]==1
            ax3.scatter(df_merged[mask]['POS_mb'], -np.log10(df_merged[mask][args.colname_pval]),
                        c=colors[c], lw=line_width, edgecolors='k', s=size, label=f'{col}: N={mask.sum()}')
            c += 1
    ax3.legend()
    ax3.set_ylabel('-log10(p value)')
    
    
    # ########## Plot PIP ##########
    ax4 = plt.subplot(414, sharex=ax2)
    plot_pip(df_data=df_merged, ax=ax4, lead_snp_pos_mb=lead_snp_pos_mb,
             lead_snp_p=lead_snp_p, pip_col=args.colname_pip)
    
    # Plot lead SNP
    ax4.scatter(lead_snp_pos_mb, lead_pip, marker='D', facecolor='None',
                c='blueviolet', lw=line_width, edgecolors='k', s=size+10, label='Lead SNP')
    
    # Highlight credible set
    # Set cycle to automaticlly change color
    colors = ['salmon', 'yellow', 'green', 'cyan', 'royalblue']*10
    c = 0
    for col in df_merged.columns:
        if '95_credible_set' in col:
            mask = df_merged[col]==1
            ax4.scatter(df_merged[mask]['POS_mb'], df_merged[mask][args.colname_pip], c=colors[c],
                        lw=line_width, edgecolors='k', s=size, label=f'{col}: N={mask.sum()}')
            c += 1
    
    # ax4.set_title('Fine-mapping PIP')
    ax4.set_ylabel('PIP')
    ax4.set_xlabel(f'Position on chr{chromosome} (Mb)')
    ax4.legend()
    
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.suptitle(f'{args.output_prefix}', y=1)
    fig.tight_layout()
    fig.savefig(output_fn)