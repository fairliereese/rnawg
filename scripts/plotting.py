import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import upsetplot
from scipy import stats
from .utils import *

def get_talon_nov_colors(cats=None):
    c_dict = {'Known': '#009E73',
              'ISM': '#0072B2',
              'NIC': '#D55E00',
              'NNC': '#E69F00',
              'Antisense': '#000000',
              'Intergenic': '#CC79A7',
              'Genomic': '#F0E442'}
    order = ['Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic', 'Genomic']
    
    if cats:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]            
    return c_dict, order

def plot_det_gene_len(df, opref='figures/'):
    pass

def plot_cell_line_tissue_det_venn(df, how='gene', 
                                   nov='Known', 
                                   opref='figures/'):
    """
    Plot a venn diagram showing how many genes / transcripts
    are detected between the cell line and tissue datasets
        
    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
            or filtered (for gene)
        how (str): Either 'gene' or 'iso'
        nov (str): Novelty, either 'Known', 'NIC', or 'NNC'
        opref (str): Where to save output figures
    """
    sns.set_context('paper', font_scale=1.8)
    
    if how == 'gene':
        df = df.loc[df.gene_novelty == 'Known']
        dataset_cols = get_dataset_cols()
        df = df[['annot_gene_id']+dataset_cols]
        df = df.groupby('annot_gene_id').sum().reset_index()
        
    tissues = get_sample_datasets('tissue')
    cell_lines = get_sample_datasets('cell_line')

    df['tissue'] = df[tissues].sum(1).astype(bool)
    df['cell_line'] = df[cell_lines].sum(1).astype(bool)
    
    if how == 'iso':
        df = df[['transcript_novelty', 'tissue', 'cell_line']]
        df = df.loc[df.transcript_novelty == nov]
    else:
        print(df.head())
        df = df[['annot_gene_id', 'tissue', 'cell_line']]
    
    known_out_green = '#90D6C3'
    known_int_green = '#009E73'
    nnc_out_gold = '#F5DFAE'
    nnc_int_gold = '#E69F00'
    nic_out_orange = '#DEA67A'
    nic_int_orange = '#D55E00'
    
    if nov == 'Known':
        out_color = known_out_green
        int_color = known_int_green
    elif nov == 'NNC':
        out_color = nnc_out_gold
        int_color = nnc_int_gold
    elif nov == 'NIC':
        out_color = nic_out_orange
        int_color = nic_int_orange
        
    df = df.groupby(['tissue', 'cell_line']).count().reset_index()
    if how == 'iso':
        df.rename({'transcript_novelty': 'counts'}, axis=1, inplace=True)
    elif how == 'gene':
        df.rename({'annot_gene_id': 'counts'}, axis=1, inplace=True)
    print(df)
    intersection = df.loc[(df.cell_line)&(df.tissue), 'counts'].values[0]
    tissue = df.loc[(~df.cell_line)&(df.tissue), 'counts'].values[0]
    cell_line = df.loc[(df.cell_line)&(~df.tissue), 'counts'].values[0]
    
    counts = [cell_line, tissue, intersection]
    log_counts = [np.log2(n) for n in counts]
    log_counts = tuple(counts)
    
    v = venn2(subsets=log_counts, set_labels=('',''))
    v.get_patch_by_id('10').set_color(out_color)
    v.get_patch_by_id('01').set_color(out_color)
    v.get_patch_by_id('11').set_color(int_color)
    v.get_patch_by_id('10').set_edgecolor(int_color)
    v.get_patch_by_id('01').set_edgecolor(int_color)
    v.get_patch_by_id('11').set_edgecolor(int_color)
    v.get_patch_by_id('10').set_linewidth(5)
    v.get_patch_by_id('01').set_linewidth(5)
    v.get_patch_by_id('11').set_linewidth(5)
    v.get_patch_by_id('10').set_alpha(1)
    v.get_patch_by_id('01').set_alpha(1)
    v.get_patch_by_id('11').set_alpha(1)
    v.get_label_by_id('10').set_text(counts[0])
    v.get_label_by_id('01').set_text(counts[1])
    v.get_label_by_id('11').set_text(counts[2])
    
    fname = '{}{}_{}_venn.png'.format(opref, how, nov)
    plt.savefig(fname, dpi=300, bbox_inches='tight') 

def plot_n_reps_per_biosamp(df,
                            sample='cell_line',
                            opref='figures/'):
    """
    Plot a bar plot showing the number of libraries
        that went into each sample 
        
    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        sample (str): Either "tissue", "cell_line"
        opref (str): Output prefix to save figure
    """
    
    dataset_cols = get_sample_datasets(sample)
    df = df[dataset_cols]
    df = df.transpose()
    df.reset_index(inplace=True)
    df.rename({'index': 'dataset'}, axis=1, inplace=True)
    df['celltype'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]

    if sample == 'tissue':
        # add in the tissue metadata
        d = os.path.dirname(__file__)
        fname = '{}/../refs/tissue_metadata.csv'.format(d)
    #     fname = '../../refs/tissue_metadata.csv'.format(d)
        tissue = pd.read_csv(fname)
        df = df.merge(tissue[['biosample', 'tissue']],
                        how='left', left_on='celltype',
                        right_on='biosample')
        df.drop('celltype', axis=1, inplace=True)
        df.rename({'tissue': 'celltype'}, axis=1, inplace=True)
        print('Found {} distinct tissues'.format(len(df.celltype.unique())))
    elif sample == 'cell_line':
        print('Found {} distinct cell lines'.format(len(df.celltype.unique()))) 

    df = df[['dataset', 'celltype']]

    ax = sns.countplot(data=df, x='celltype')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    if sample == 'tissue':
        ylabel = '# libraries'
        xlabel = 'Tissue'
    elif sample == 'cell_line':
        ylabel = '# libraries'
        xlabel = 'Cell line'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    
    fname = '{}{}_libs_per_biosamp.png'.format(opref, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
def plot_exons_per_iso(df,
                       nov=['Known'],
                       min_tpm=1,
                       gene_subset=None,
                       opref='figures/'):
    """
    Plots a boxplot showing # exons per isoform.
    
    Parameters:
        df (pandas DataFrame): TALON abundance file
        nov (list of str): Novelty categories to include
        min_tpm (float): Mininmum TPM to include isoform
        gene_subset (str): Choose from 'polya' or None
        opref (str): Output file prefix
    """
    
    t_df = df.copy(deep=True)
    
    # filter based on TPM, novelty, gene subset
    df, tids = get_tpm_table(df,
                             how='iso',
                             nov=nov,
                             min_tpm=min_tpm,
                             gene_subset=gene_subset)

    # add back in novelty and n exon info 
    t_df = t_df[['annot_transcript_id', 'n_exons', 'transcript_novelty']]
    df = df.merge(t_df, how='left', left_index=True, right_on='annot_transcript_id')   
    
    # plot the plot
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(4,6))

    c_dict, order = get_talon_nov_colors(cats=nov)
    ax = sns.boxplot(data=df, x='transcript_novelty', y='n_exons',
                     order=order, palette=c_dict,
                     saturation=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Transcript novelty'
    ylabel = '# exons'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}_exons_per_iso.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_gene_v_iso_sample_det(df,
                               sample='cell_line', 
                               opref='figures/'):
    """
    Plot a hexbin density plot of the number of samples each gene and transcript is detected per librarygenes detected per library, by novelty type

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        sample (str): Either "tissue", "cell_line"
        opref (str): Output prefix to save figure
    """
    df = rm_sirv_ercc(df)
    dataset_cols = get_sample_datasets(sample)
    gridsize = int(len(dataset_cols)/2)
    print('Gridsize: {}'.format(gridsize))

    t_df = df[dataset_cols+['annot_transcript_id']].copy(deep=True)
    t_df.set_index('annot_transcript_id', inplace=True)
    t_df = t_df.astype(bool)
    t_df['n_samples_transcript'] = t_df.sum(1)
    t_df = t_df['n_samples_transcript'].to_frame()
    t_df.reset_index(inplace=True)

    g_df = df.loc[df.gene_novelty == 'Known'].copy(deep=True)
    g_df = g_df[dataset_cols+['annot_gene_id']]
    g_df = g_df.groupby('annot_gene_id').sum()
    g_df = g_df.astype(bool)
    g_df['n_samples_gene'] = g_df.sum(1)
    g_df = g_df['n_samples_gene'].to_frame()
    g_df.reset_index(inplace=True)

    cols = ['annot_gene_id', 'annot_transcript_id', 'transcript_novelty']
    t_df = t_df.merge(df[cols], how='left', on='annot_transcript_id')
    t_df = t_df.merge(g_df, how='left', on='annot_gene_id')

    # c_dict, order = get_talon_nov_colors(['Known', 'NIC', 'NNC'])
    c_dict, order = get_talon_nov_colors()
    sns.set_context('paper', font_scale=1.6)

    for nov in ['Known', 'NIC', 'NNC']:
        temp = t_df.loc[t_df.transcript_novelty == nov].copy(deep=True)
        ax = sns.jointplot(data=temp, x='n_samples_transcript', 
                         y='n_samples_gene',
                         kind='hex',
                         color=c_dict[nov],
                         bins='log',
                         gridsize=gridsize)
        ax = ax.ax_joint

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if sample == 'tissue':
            ylabel = '# tissue libraries gene was detected in'
            xlabel = '# tissue libraries transcript was detected in'
        elif sample == 'cell_line':
            ylabel = '# cell line libraries gene was detected in'
            xlabel = '# cell line libraries transcript was detected in'

        _ = ax.set(xlabel=xlabel, ylabel=ylabel)

        fname = '{}{}_gene_v_{}_iso_n_samp_det.png'.format(opref, sample, nov)
        plt.savefig(fname, dpi=300, bbox_inches='tight')


def plot_gene_v_iso_det(df, filt_df,
                        sample='cell_line',
                        opref='figures/'):
    """
    Plot a scatterplot of the number of genes detected per library
        versus the number of transcripts detected, by novelty type

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        filt_df (pandas DataFrame): TALON abundance, filtered
        sample (str): Either "tissue", "cell_line"
        opref (str): Output prefix to save figure
    """
        
    df = rm_sirv_ercc(df)
    filt_df = rm_sirv_ercc(filt_df)
    dataset_cols = get_sample_datasets(sample)

    # only known genes
    gene_df = df.loc[df.gene_novelty == 'Known'].copy(deep=True)
    gene_df = gene_df[dataset_cols+['annot_gene_id']]
    gene_df = gene_df.groupby('annot_gene_id').sum()

    gene_df = gene_df.astype(bool)
    ind_cols = gene_df.columns
    gene_df = gene_df.transpose()
    gene_df['n_genes'] = gene_df.sum(1)
    gene_df = gene_df['n_genes'].to_frame()

    # known transcripts
    t_df = filt_df.loc[filt_df.transcript_novelty == 'Known'].copy(deep=True)
    t_df = t_df[dataset_cols+['annot_transcript_id']]
    t_df = t_df.astype(bool)
    ind_cols = t_df.columns
    t_df = t_df.transpose()
    t_df['Known'] = t_df.sum(1)
    t_df = t_df['Known'].to_frame()

    gene_df = gene_df.merge(t_df, left_index=True, right_index=True)

    # nic transcripts
    t_df = filt_df.loc[filt_df.transcript_novelty == 'NIC'].copy(deep=True)
    t_df = t_df[dataset_cols+['annot_transcript_id']]
    t_df = t_df.astype(bool)
    ind_cols = t_df.columns
    t_df = t_df.transpose()
    t_df['NIC'] = t_df.sum(1)
    t_df = t_df['NIC'].to_frame()

    gene_df = gene_df.merge(t_df, left_index=True, right_index=True)

    # nnc transcripts
    t_df = filt_df.loc[filt_df.transcript_novelty == 'NNC'].copy(deep=True)
    t_df = t_df[dataset_cols+['annot_transcript_id']]
    t_df = t_df.astype(bool)
    ind_cols = t_df.columns
    t_df = t_df.transpose()
    t_df['NNC'] = t_df.sum(1)
    t_df = t_df['NNC'].to_frame()

    gene_df = gene_df.merge(t_df, left_index=True, right_index=True)
    gene_df.reset_index(inplace=True)
    
    df = gene_df.melt(id_vars=['index', 'n_genes'],
                  value_vars=['Known','NIC','NNC'])
    df.rename({'value': 'transcript_counts'}, axis=1, inplace=True)
    df.rename({'variable': 'novelty'}, axis=1, inplace=True)
    
    c_dict, order = get_talon_nov_colors(['Known', 'NIC', 'NNC'])
    sns.set_context('paper', font_scale=1.6)

    ax = sns.jointplot(data=df, x='transcript_counts', y='n_genes',
                     hue='novelty', palette=c_dict,
    #                  xlim=(0,xlim), ylim=(0,ylim), 
                     joint_kws={'data':df, 's':40, 'alpha':1})
    ax = ax.ax_joint

    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if sample == 'tissue':
        ylabel = 'Known genes per tissue library'
        xlabel = 'Transcripts per tissue library'
    elif sample == 'cell_line':
        ylabel = 'Known genes per cell line library'
        xlabel = 'Transcripts per cell line library'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    
    fname = '{}{}_gene_v_iso_det.png'.format(opref, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight') 

def plot_avg_isos_per_gene(df,
                           min_tpm=1,
                           gene_subset='polya',
                           sample='all', 
                           groupby='sample',
                           nov=['Known', 'NIC', 'NNC'],
                           opref='figures/'):
    """
    Plot the average number of isoforms per gene that are seen in each
    sample or library.
    
    Parameters:
        df (pandas DataFrame): TALON abundance
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library", 
            used to groupby datasets displayed
        nov (str): Novelty category of 
            isoforms to consider
        opref (str): Output prefix to save figure
        
    Returns: 
        df (pandas DataFrame): DataFrame detailing how many samples
            each gene or isoform was seen in
    """
    
    sns.set_context('paper', font_scale=1.6)
    
    # get # isos / gene / sample or library
    df = get_isos_per_gene(df,
                           min_tpm=min_tpm,
                           gene_subset=gene_subset,
                           groupby=groupby, 
                           nov=nov)
    
    # fill 0s (which represent unexpresed genes) with NaNs before
    # calculating the avg. # isos / sample
    df = df.replace(0, np.nan)

    # calculate the average to order the barplots
    avgs = df.mean().to_frame()    
    avgs = avgs.sort_values(by=[0], ascending=False)
    order = avgs.index.tolist()
    
    # melt the df to get an entry for each sample and gene
    df = df.melt()    
    df.rename({'variable': 'tissue',
               'value': 'n_isos'}, axis=1, inplace=True) 
    df = df.loc[~df.n_isos.isnull()]   
    
    # actually make the plot
    ax = sns.barplot(data=df, x='tissue', y='n_isos', order=order)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Sample'
    ylabel = '# isoforms / gene'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(axis="x", rotation=90)


    fname = '{}_isos_per_gene_per_{}.png'.format(opref, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
def plot_biosamp_det(df,
                     how='gene',
                     min_tpm=1, 
                     gene_subset='polya',
                     sample='cell_line',
                     groupby='cell_line',
                     nov='Known',
                     opref='figures/'):
    
    """
    Plot a hist of the number of tissues, cell lines, or datasets each gene or 
    isoform occurs in.
    
    Parameters:
        df (pandas DataFrame): TALON abundance
        how (str): Either "gene" or "iso"
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library", 
            used to groupby datasets displayed
        nov (str): Only used with how='iso', novelty category of 
            isoforms to consider
        opref (str): Output prefix to save figure
        
    Returns: 
        df (pandas DataFrame): DataFrame detailing how many samples
            each gene or isoform was seen in
    """
    
    df = get_det_table(df,
                     how=how,
                     min_tpm=min_tpm,
                     gene_subset=gene_subset,
                     sample=sample,
                     groupby=groupby,
                     nov=[nov])
    
#     # calc TPM per library on desired samples
#     df, tids = get_tpm_table(df,
#                    sample=sample,
#                    how=how,
#                    nov=[nov],
#                    min_tpm=min_tpm,
#                    gene_subset=gene_subset)
    
#     df = df.transpose()
#     df.index.name = 'dataset'
#     df.reset_index(inplace=True)

#     # set up df to groupby sample or library
#     if groupby == 'sample':

#         # add biosample name (ie without rep information)
#         df['biosample'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]
#         df.drop(['dataset'], axis=1, inplace=True)

#         # record the highest TPM value per biosample
#         tissue_df = get_tissue_metadata()
#         tissue_df = tissue_df[['tissue', 'biosample']]

#         df = df.merge(tissue_df, how='left', on='biosample')
#         df.loc[df.tissue.isnull(), 'tissue'] = df.loc[df.tissue.isnull(), 'biosample']
#         df.drop('biosample', axis=1, inplace=True)
#         df.rename({'tissue': 'biosample'}, axis=1, inplace=True)

#         print('Found {} total samples'.format(len(df.biosample.unique().tolist())))
#         df = df.groupby('biosample').max()

#     elif groupby == 'library':
#         df.rename({'dataset': 'library'}, axis=1, inplace=True)
#         print('Found {} total libraries'.format(len(df.library.unique().tolist())))
#         df = df.groupby('library').max()
    
    # finally, calculate the number of biosamples / libraries these 
    # genes or transcripts are expressed >= min TPM
    df = df.transpose()
    df['n_samples'] = df.astype(int).sum(axis=1)
    
    # and make a beautiful plot
    sns.set_context('paper', font_scale=2)
    
    c_dict, order = get_talon_nov_colors()
    color = c_dict[nov]
    ax = sns.displot(data=df, x='n_samples', kind='hist',
                 color=color, binwidth=1)

    # titles
    if how == 'gene':
        ylabel = 'Known genes'
    elif how == 'iso':
        if nov == 'Known':
            ylabel = 'Known transcripts'
        elif nov == 'NIC':
            ylabel = 'NIC transcripts'
        elif nov == 'NNC':
            ylabel = 'NNC transcripts'

    if groupby == 'sample':
        xlabel = 'Number of cell lines and tissues'
    elif groupby == 'tissue':
        xlabel = 'Number of tissues'
    elif groupby == 'cell_line':
        xlabel = 'Number of celltypes'
    elif groupby == 'library':
        if sample == 'cell_line':
            xlabel = 'Number of cell line libraries'
        elif sample == 'tissue':
            xlabel = 'Number of tissue libraries'
        else: 
            xlabel = 'Number of libraries'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    if groupby == sample:
        fname = '{}{}_{}_{}_detection.png'.format(opref, sample, nov, how)
    else:
        fname = '{}{}_{}_{}_library_detection.png'.format(opref, sample, nov, how)

    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    return df

def plot_corr(df, sample='cell_line',
              how='gene', nov='Known', 
              opref='figures/', cluster=False):
    
    corrs = compute_corr(df, how=how, sample=sample)
    if cluster == False:
        g = sns.clustermap(corrs, cmap='viridis',
                        xticklabels=True, yticklabels=True,
                        row_cluster=False, col_cluster=False)
    else:
        g = sns.clustermap(corrs, cmap='viridis',
                xticklabels=True, yticklabels=True)
    
    if how == 'iso':
        if nov == 'Known':
            title_txt = 'Known transcript'
        else:
            title_txt = '{} transcript'.format(nov)
    elif how == 'gene':
        title_txt = 'Known gene'
    title = '{} expression correlations'.format(title_txt)
             
    g.fig.suptitle(title) 
    
    fname = '{}{}_{}_correlation.pdf'.format(opref, nov, how)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    fname = '{}{}_{}_correlation.png'.format(opref, nov, how)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    fname = '{}{}_{}_correlation.tsv'.format(opref, nov, how)
    corrs.to_csv(fname, sep='\t')
    

def plot_ranked_biosamp(df, sample='cell_line', how='iso', nov='known',
                    ylim=None, opref='figures/'):
    sns.set_context('paper', font_scale=2)

    max_df, min_df = compute_detection(df, sample, how, nov)

    c_dict, order = get_talon_nov_colors()
    color = c_dict[nov]
    ax = sns.lineplot(data=max_df, x='rank',
                      y='n_cumulative', linewidth=3, color=color)
    ax = sns.lineplot(data=min_df, x='rank',
                      y='n_cumulative', linewidth=3,
                      ax=ax, color=color)

#     if how == 'gene':
#         # n gencode genes
#         ax.plot([0, 23], [58780,58780], color='gray',
#                 linestyle='--', linewidth=3)
#         font = {'color': 'gray', 'size': 16}
#         plt.text(9.5, 61500, "~60k GENCODE genes", fontdict=font)
#     if how == 'iso' and nov == 'Known':
#         # n gencode genes
#         ax.plot([0, 23], [206761,206761], color='gray',
#                 linestyle='--', linewidth=3)
#         font = {'color': 'gray', 'size': 16}
#         plt.text(5.75, 213000, "~200k GENCODE transcripts", fontdict=font)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if how == 'iso':
        if nov == 'Known':
            ylabel_txt = 'known isoforms'
        else:
            ylabel_txt = '{} isoforms'.format(nov)
    elif how == 'gene':
        ylabel_txt = 'known genes'

    ylabel = 'Cumulative {} detected'.format(ylabel_txt)
    xlabel = 'Ranked {}s'.format(sample)
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    if ylim:
        _ = ax.set(ylim=(-100, ylim))

    fname = '{}cumulative_{}_{}_per_{}.pdf'.format(opref, nov, how, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    fname = '{}cumulative_{}_{}_per_{}.png'.format(opref, nov, how, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    return max_df, min_df

def plot_short_long_det(df, opref, \
                    xlim, ylim, how='gene'):

    sns.set_context('paper', font_scale=2)

    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'ill_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'ill_umi_count'

#     ax = sns.jointplot(data=df, x=c1, y=c2,
#                      xlim=(0,xlim), ylim=(0,ylim),
#                      joint_kws={'data':df, 's':40, 'alpha':1})
    ax = sns.jointplot(data=df, x=c1, y=c2,
                     xlim=(0,xlim), ylim=(0,ylim),
                     joint_kws={'data':df, 's':40, 'alpha':1})
    ax = ax.ax_joint

#     # plot regression lines and equation of regression lines
#     # https://stackoverflow.com/questions/48145924/different-colors-for-points-and-line-in-seaborn-regplot/68135585#68135585
#     # https://stackoverflow.com/questions/45902739/seaborn-annotate-the-linear-regression-equation
#     # https://stackoverflow.com/questions/62705904/add-entry-to-matplotlib-legend-without-plotting-an-object
#     lines = []
#     labels = []
#     for s in df['sample'].unique().tolist():
#         temp = df.loc[df['sample'] == s]
#         color = c_dict[s]
#         line_color = adjust_lightness(color, 0.5)

#         # get coeffs of linear fit
#         slope, intercept, r_value, p_value, std_err = stats.linregress(temp[c1],temp[c2])
#         lines += [mpl.lines.Line2D([0], [0], color=line_color)]
#         labels += ['m={0:.1f}'.format(slope)]

#         print('Slope of {} correlation: {}'.format(s, slope))

#         sns.regplot(data=temp, x=c1, y=c2,
#                     scatter=False, ax=ax, color=color)
#         sns.regplot(data=temp, x=c1, y=c2,
#             scatter=False, ax=ax, color=color, ci=0,
#             line_kws={'color':line_color,
#                       'linestyle':'-',
#                       'label':"m={0:.1f}".format(slope)})

    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')


def plot_reads_per_bc(df, title, oprefix):
    """
    Parameters:
        df (pandas DataFrame): DataFrame of read_annot file
        title (str): Title of plot
        oprefix (str): Output file prefix
    """

    temp = get_reads_per_bc(df)

    sns.set_context('paper', font_scale=1.5)
    counts = temp['counts'].tolist()
    plt.plot(range(len(counts)),
            counts,
            color='lightgray',
            linewidth=2)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_xlabel('Ranked cells by # reads (logscale)')
    ax.set_ylabel('# reads (logscale)')
    ax.set_title(title)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if title == 'Post-TALON':
        title = 'post_talon'

    plt.tight_layout()

    fname = '{}_{}_umis_v_barcodes.png'.format(oprefix, title)
    plt.savefig(fname)
    
def plot_det_len_kde(df,
                     how='gene',
                     subset='polya',
                     min_tpm=1,
                     xlim=None,
                     split_biotypes=False,
                     opref='figures/'):
    """
    Plots dist. of gene or transcript length based on whether 
    gene or transcript was detected at the given TPM value.
    
    Parameters:
        df (pandas DataFrame): TALON abundance file
        how (str): Choose from 'gene' or 'iso'
        subset (str): Choose from None or 'polya'
        min_tpm (float): Min. TPM val for at least one library
        xlim (float): Maximum length to display
        split_biotypes (bool): Split detected and undetected 
            transcripts by biotype
        opref (str): Output file prefix
        
    Returns:
        df (pandas DataFrame): DataFrame used to plot from
    """
    df, ids = get_tpm_table(df,
                   how=how,
                   min_tpm=min_tpm,
                   nov=['Known'],
                   gene_subset=subset)
    gene_df, _, _ = get_gtf_info(how=how, subset=subset)   
    df.reset_index(inplace=True)
    
    if how == 'gene':
        col = 'annot_gene_id'
        ref_col = 'gid'
        x_col = 'length'
    elif how == 'iso':
        col = 'annot_transcript_id'
        ref_col = 'tid'
        x_col = 't_len'
    df = df[col].to_frame()

    df = df.merge(gene_df, how='outer',
                  left_on=col, right_on=ref_col)
    
    df['detected'] = True
    df.loc[df[col].isnull(), 'detected'] = False
    
    sns.set_context('paper', font_scale=2)

    ax = sns.displot(data=df, x=x_col, kind='kde',
                     linewidth=3, hue='detected', common_norm=True)

    if how == 'gene':
        xlabel = 'Gene length'
        ylabel = 'Detected GENCODE genes'
    elif how == 'iso':
        xlabel = 'Transcript length'
        ylabel = 'Detected GENCODE transcripts'
        
    if xlim:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(0,xlim))
    else:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel)
        
    plt.savefig('{}_{}_det_{}_tpm_length.pdf'.format(opref, how, min_tpm), \
                dpi=300, bbox_inches='tight')
    return df

def plot_read_len_kde(df, opref):
    sns.set_context('paper', font_scale=2)

#     ax = sns.displot(data=df, x='read_length', hue=hue,
#                  palette=c_dict, kind='kde', hue_order=order, linewidth=3,
#                  common_norm=common_norm)
    ax = sns.displot(data=df, x='read_length', kind='kde', linewidth=3)
#     ax.set(xlabel='Read length', ylabel='KDE of reads',
#           title='Length distribution of Reads', xlim=(0,7500),
#           xticks=[0, 2500, 5000, 7500])
    plt.savefig('{}_read_length_kde.pdf'.format(opref), dpi=300, bbox_inches='tight')

def plot_cluster_proportions(adata,
                             bar_key,
                             color_key,
                             opref='figures/'):

    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([color_key, bar_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=bar_key, index=color_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)

    fig, ax = plt.subplots(dpi=300)
    fig.patch.set_facecolor("white")
    fig.set_size_inches(15, 5)

    cmap = None
    cluster_palette = '{}_colors'.format(color_key)
    if cluster_palette in adata.uns.keys():
        cluster_palette = adata.uns[cluster_palette]
        cmap = sns.palettes.blend_palette(
            cluster_palette,
            n_colors=len(cluster_palette),
            as_cmap=True)

    props.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        legend=None,
        colormap=cmap
    )

    ax.legend(bbox_to_anchor=(1.1, 2), title=color_key)
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=90)
    ax.set_xlabel(props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    ax.grid(False)

    fname = opref+'{}_by_{}_prop.png'.format(bar_key, color_key)
    plt.savefig(fname, bbox_inches='tight')

def plot_depth_by_tech(adata, how, opref,
                       hue=None, xlim=None, ylim=None):
    df = adata.obs
    sns.set_context('paper', font_scale=2)

    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'sr_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'sr_umi_count'

    if xlim:
        xlim = (0, xlim)
    if ylim:
        ylim = (0, ylim)
    ax = sns.jointplot(data=df, x=c1, y=c2,
                       xlim=xlim, ylim=ylim,
                       joint_kws={'data':df, 's':40, 'alpha':1},
                       hue=hue)
    ax = ax.ax_joint

    # ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')

def add_perc(ax, data, feature):
    total = data[feature].sum()
    ylim = ax.get_ylim()[1]
    n_cats = len(ax.patches)
    for p in ax.patches:
        percentage = '{:.1f}%'.format(100 * p.get_height()/total)
#         x = p.get_x() + p.get_width() / 2 - 0.45
        x = p.get_x() + p.get_width() / 2 - (0.065)*n_cats
        y = p.get_y() + p.get_height() + ylim*0.00625
        ax.annotate(percentage, (x, y), size = 12)

def plot_read_novelty(df, opref, c_dict, order,
                      ylim=None, title=None,
                      datasets='all'):
    sns.set_context("paper", font_scale=1.6)

    temp = df.copy(deep=True)

    # filter on datasets
    if datasets != 'all':
        temp = temp.loc[temp.dataset.isin(datasets)]

    # count number of reads per cat
    temp = temp[['transcript_novelty', 'read_name']].groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'read_name':'counts'}, axis=1, inplace=True)
    print(temp)

    # actual plotting
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar',
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Reads')
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

    # add title
    if not title:
        g.fig.suptitle('Reads per novelty category')
    else:
        g.fig.suptitle('{} reads per novelty category'.format(title))

    # save figure
    fname = '{}_read_novelty'.format(opref)
    g.savefig(fname+'.pdf', dpi=300)


def plot_transcript_novelty(df, oprefix,
                            ylim=None,
                            title=None,
                            whitelist=None,
                            sample=None,
                            novs=None,
                            save_type='pdf'):
    """
    Plot number of transcripts per novelty category.
    
    Parameters: 
        df (pandas DataFrame): TALON read annot file or 
        oprefix (str): Place to save
        ylim (int): y limit of resultant plot
        title (str): Title of resultant plot
        whitelist (list of str): List of transcript IDs to retain
        sample (str): Choose from 'cell_line' or 'tissue'
        save_type (str): Choose from 'pdf' or 'png'
    """
    sns.set_context('paper', font_scale=1.6)

    temp = df.copy(deep=True)
    
    c_dict, order = get_talon_nov_colors(cats=novs)

    # remove transcripts that are not on whitelist
    if whitelist:
        beep = temp.loc[temp.transcript_ID.isin(whitelist)].index.tolist()
        if len(beep) == 0:
            beep = temp.loc[temp.annot_transcript_id.isin(whitelist)].index.tolist()
        temp = temp.loc[beep]
        
    # filter on datasets
    if sample:
        datasets = get_sample_datasets(sample)
    else:
        datasets = get_dataset_cols()
    cols = ['transcript_ID', 'transcript_novelty']
    temp = temp[cols+datasets]
    
    temp['total_counts'] = temp[datasets].sum(1)
    temp = temp.loc[temp.total_counts > 0]

    # count number of isoforms per cat
#     temp = temp[['transcript_ID', 'transcript_novelty', 'read_name']].groupby(['transcript_ID', 'transcript_novelty']).count()
#     temp.reset_index(inplace=True)
#     temp.drop('read_name', axis=1, inplace=True)
    temp = temp[['transcript_ID', 'transcript_novelty']]
    temp = temp.groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'transcript_ID': 'counts'}, axis=1, inplace=True)
    print(temp)
    novs = ['NIC', 'Known', 'NNC']
    complete = temp.loc[temp.transcript_novelty.isin(novs), 'counts'].sum(axis=0)
    print('Number of complete isoforms: {}'.format(complete))
    

    # actual plotting
    sns.set_context('paper', font_scale=1.8)
    plt.figure(figsize=(4,6))
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar',
                saturation=1,
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Isoforms')
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

#     # add title
#     if not title:
#         g.fig.suptitle('Transcript models per novelty category')
#     else:
#         g.fig.suptitle('{} transcript models per novelty category'.format(title))

    # save figure
    fname = '{}{}_isoform_novelty'.format(oprefix,sample)
    if save_type == 'png':
        g.savefig(fname+'.png', dpi=300, bbox_inches='tight')
    elif save_type == 'pdf':
        g.savefig(fname+'.pdf', dpi=300, bbox_inches='tight')

    plt.show()
    plt.clf()


    