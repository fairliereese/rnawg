import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import upsetplot
from scipy import stats
from matplotlib.ticker import ScalarFormatter
import swan_vis as swan
import ternary
from .utils import *

def get_edge_colors():
    """
    Get colors for introns and exons
    """
    c_dict = {'intron': '#CC79A7', 'exon': '#009E73'}
    return c_dict

def get_biosample_colors():
    """
    Get colors for each biosample
    """
    d = os.path.dirname(__file__)
    fname = '{}/../refs/biosample_colors.tsv'.format(d)
    df = pd.read_csv(fname, sep='\t')
    
    c_dict = {}
    for ind, entry in df.iterrows():
        c_dict[entry.biosample] = entry.color
    order = df.biosample.tolist()
    
    return c_dict, order
    
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

def plot_dorito(counts,
                gene=None,
                sources=['GENCODE', 'obs'],
                color='source',
                c_dict=None,
                order=order,
                scale=True,
                top='splicing_ratio',
                opref='figures/'):
    """
    Plot a dorito from the given sources with the given color
    
    Parameters:
        counts (pandas DataFrame): DF of the counts per gene 
            of ic, tss, tes from get_ic_tss_tes or
            compute_triplets (or both!!!)
        gname (str): Gene to plot. Default=None
        sources (list of str): List of columns to include in graph
        color (str): Column from counts to color the points by
        c_dict (dict of str): Dict indexed by unique values
            in the color column of counts
        scale (bool): Whether or not to scale each axis by sum of 
            three numbers
        top (str): Column name to plot as apex of dorito.
            Choose from 'ic' or 'splicing_ratio'
        opref (str): Output file path / prefix to save fig
    """
    
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
    
def plot_region_widths(regions,
                       kind='annot',
                       opref='figures/human'):
    """
    Parameters:
        regions (dict of pandas DataFrame): Output from 
            get_ic_tss_tes
        kind (str): Choose 'annot', 'obs', or 'all'
        opref (str): Output file prefix
    """
    
    # plot histogram of tss / tes region sizes
    sns.set_context('paper', font_scale=1.8)
    for c, temp in regions.items():
        ax = sns.displot(data=temp, x='len', kind='hist', linewidth=0, binwidth=50)
        ax.set(xlabel='{}s'.format(c.upper()))
        fname = '{}_{}_{}_region_widths.png'.format(opref, kind, c)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_genes_n_ic_ends(counts, 
                         kind='annot',
                         opref='figures/human'):
    """
    Parameters:
        counts (pandas DataFrame): DF output from get_ic_tss_tes
        kind (str): Choose from 'annot', 'all', 'obs'
        opref (str): Where to save
    """
    # plot # tss / tes vs ic, color by # genes 
    counts['tss_tes'] = counts.tss+counts.tes
    counts.head()

    # calculate how many genes have n tss+tes and n ics
    temp = counts[['tss_tes', 'intron_chain', 'gid']].groupby(['tss_tes', 'intron_chain']).count().reset_index()
    temp.rename({'gid': 'n_genes'}, inplace=True, axis=1)
    temp['log10_n_genes'] = np.log10(temp.n_genes)

    # plot the figure
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(6,8))

    # ax = sns.scatterplot(data=temp, x='tss_tes', y='intron_chain', hue='n_genes', size='n_genes', palette='viridis')
    ax = sns.scatterplot(data=temp, x='tss_tes', y='intron_chain', hue='log10_n_genes', size='log10_n_genes', palette='viridis')

    norm = plt.Normalize(temp['log10_n_genes'].min(), temp['log10_n_genes'].max())
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar
    ax.get_legend().remove()
    cb = ax.figure.colorbar(sm)
    cb.set_label('log10(# genes)')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    xlabel = '# TSSs + # TESs'
    ylabel = '# intron chains'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel, xscale='log', yscale='log')

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    fname = '{}_{}_n_genes_ic_tss_tes.png'.format(opref, kind)
    plt.savefig(fname, dpi=300, bbox_inches='tight')               
        
def plot_n_ic_tss_tes(counts,
                      label_genes=None,
                      kind='annot',
                      opref='figures/human'):
    """
    Parameters:
        counts (pandas DataFrame): DF output from get_ic_tss_tes
        label_genes (list of str): List of gene names
        kind (str): Choose from 'annot', 'all', 'obs'
        opref (str): Where to save thing
    """
    
    xs = ['tss', 'tss', 'tes']
    ys = ['intron_chain', 'tes', 'intron_chain']
    hues = ['tes', 'intron_chain', 'tss']
    
    for x, y, hue in zip(xs, ys, hues):
        
        # plot the figure
        sns.set_context('paper', font_scale=1.6)
        plt.figure(figsize=(6,8))

        ax = sns.scatterplot(data=counts, x=x, y=y, hue=hue, s=20, palette='viridis')

        norm = plt.Normalize(counts[hue].min(), counts[hue].max())
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])

        # Remove the legend and add a colorbar
        ax.get_legend().remove()
        cb = ax.figure.colorbar(sm)
        if hue == 'tss' or hue == 'tes':
            cb.set_label('# {}s'.format(hue.upper()))
        else:
            cb.set_label('# {}s'.format(hue))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if x == 'tss' or x == 'tes':
            xlabel = '# {}s'.format(x.upper())
        else:
            xlabel = '# {}s'.format(x)
        if y == 'tss' or y == 'tes':
            ylabel = '# {}s'.format(y.upper())
        else:
            ylabel = '# {}s'.format(y)

        # annotate genes that are kinda interesting
        if label_genes:
            xlim = ax.get_xlim()[1]
            ylim = ax.get_ylim()[1]
            for g in label_genes:
                if g in counts.gname.tolist():
                    x_txt = counts.loc[counts.gname == g, x].values[0]+(1/80)*xlim
                    y_txt = counts.loc[counts.gname == g, y].values[0]-(1/80)*ylim
                    plt.annotate(g, (x_txt,y_txt), fontsize='small', fontstyle='italic')
        _ = ax.set(xlabel=xlabel, ylabel=ylabel)
        
        fname = '{}_{}_{}_{}_{}_scatter.png'.format(opref, x,y,hue, kind)
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
    
def plot_major_iso_pis(sg, groupby, 
                       opref='figures/human'):
    """
    Plots a histogram of pi values for the 1st-4th highest-expressed isoform 
    per gene per sample.
    
    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with abundance information
        groupby (str): Column in sg.adata.obs on which to calculat pi
        opref (str): Output file prefix
    """
    obs_col = groupby
    df, counts = swan.calc_pi(sg.adata, sg.t_df, obs_col=obs_col)
    df = df.transpose()

    # merge with gene info
    df = df.merge(sg.t_df['gid'],
                  how='left', left_index=True, right_index=True)

    # get top isoform per gene
    df.reset_index(inplace=True)
    df = df.melt(id_vars=['tid', 'gid'])

    df.rename({'variable': obs_col,
               'value': 'pi'}, axis=1, inplace=True)
    df = df.sort_values(by='pi', ascending=False)

    # remove unexpressed transcripts
    df = df.loc[df.pi != 0]

    # add # isoforms / gene / dataset to use for filtering later
    temp = df[['gid', obs_col, 'tid']].groupby(['gid', obs_col]).count().reset_index()
    temp.rename({'tid': 'iso_counts'}, axis=1, inplace=True)
    df = df.merge(temp, how='left', on=['gid', obs_col])


    for i in range(1,5):
        # https://stackoverflow.com/questions/36310564/pandas-second-max-value-per-group-in-dataframe
        temp = df.groupby(['gid', obs_col]).head(i).groupby(['gid', obs_col]).tail(1).copy(deep=True)

        # make sure we're only looking at isoforms that have at least
        # i expressed isoforms
        temp = temp.loc[temp.iso_counts >= i]

        # plot distribution of max pi value
        # per gene, per tissue / age combo
        sns.set_context('paper', font_scale=2)
        ax = sns.displot(data=temp, x='pi', linewidth=0)

        if i == 1: 
            xlabel = 'Highest pi per gene per {}'.format(obs_col)
        elif i == 2:
            xlabel = '2nd highest pi per gene per {}'.format(obs_col)
        elif i == 3:
            xlabel = '3rd highest pi per gene per {}'.format(obs_col)
        elif i > 3:
            xlabel = '{}th highest pi per gene per {}'.format(i, obs_col)
        ylabel = 'Number of isoforms'
        ax.set(ylabel=ylabel, xlabel=xlabel, xlim=(0,100))

    

def plot_ranked_exon_counts(sg,
                          df, 
                          gene, 
                          min_tpm=1,
                          gene_subset='polya',
                          sample='all',
                          groupby='library',
                          nov=['Known', 'NIC', 'NNC'], 
                          opref='figures/human'):
    
    """
    Plot the ranked counts per exon for a given gene subset according
    to the input expression threshold and novelty categories
    
    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with data from corresponding
    """
    
    # determine which isoforms are actually detected
    df = get_det_table(df, 
                       how='iso',
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample, 
                       groupby=groupby,
                       nov=nov)
    tids = df.columns.tolist()

    # get isoforms from target gene
    df = sg.t_df.loc[sg.t_df.gname == gene]
    df = swan.pivot_path_list(df, 'path')
    df = df.merge(sg.edge_df, how='left', left_on='edge_id', right_index=True, suffixes=(None, '_dupe'))
    df = df.merge(sg.t_df[['tname']], how='left', left_index=True, right_index=True)
    df.drop('edge_id_dupe', axis=1, inplace=True)
    df = df[['edge_id', 'edge_type']]
    df.reset_index(inplace=True)
    tids = list(set(tids)&set(df.tid.tolist()))

    # limit to only detected isoforms and their edges
    df = df.loc[df.tid.isin(tids)]
    eids = df.edge_id.astype('int').tolist()
    df = sg.get_edge_abundance(kind='counts')
    df = df.loc[df.edge_id.isin(eids)]

    # sum up over datasets
    df.set_index('edge_type', inplace=True)
    df = df[sg.datasets]
    df['total_counts'] = df.sum(1)
    df.drop(sg.datasets, axis=1, inplace=True)

    # rank according to exp
    df = df.sort_values(by='total_counts', ascending=True)
    df.reset_index(inplace=True)
    
    sns.set_context('paper', font_scale=1.8)
    c_dict = get_edge_colors()

    for e in df.edge_type.unique():
        temp = df.loc[df.edge_type == e]
        temp['rank'] = [i for i in range(len(temp.index))]
        c = c_dict[e]

        ax = sns.catplot(data=temp, x='rank', y='total_counts', kind='bar',
                         linewidth=0, saturation=1, color=c)

        xlabel = 'Ranked $\it{}$ {}s'.format(gene, e)
        ylabel = 'Total counts'

        _ = ax.set(xlabel=xlabel, ylabel=ylabel, xticks=[])
        fname = '{}_{}_counts_per_{}.png'.format(opref, gene, e)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

    return df 

def plot_exon_hist(sg,
                   df,
                   gene,
                   min_tpm=1,
                   gene_subset='polya',
                   sample='all',
                   groupby='library',
                   nov=['Known', 'NIC', 'NNC'], 
                   opref='figures/human'):
    """
    Plot a histogram of introns and exons and how many isoforms
    they're used in in a particular gene.
    
    Parameters:
        sg (swan_vis SwanGraph): SwanGraph of data
        df (pandas DataFrame): TALON abundance
    """
    # determine which isoforms are actually detected
    df = get_det_table(df, 
                       how='iso',
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample, 
                       groupby=groupby,
                       nov=nov)
    tids = df.columns.tolist()
    
    # get isoforms from target gene
    df = sg.t_df.loc[sg.t_df.gname == gene]
    df = swan.pivot_path_list(df, 'path')
    df = df.merge(sg.edge_df, how='left', left_on='edge_id', right_index=True, suffixes=(None, '_dupe'))
    df = df.merge(sg.t_df[['tname']], how='left', left_index=True, right_index=True)
    df.drop('edge_id_dupe', axis=1, inplace=True)
    df = df[['edge_id', 'edge_type']]
    df.reset_index(inplace=True)
    tids = list(set(tids)&set(df.tid.tolist()))
    
    # limit to only detected isoforms
    df = df.loc[df.tid.isin(tids)]
    
    print('Found {} isoforms for {}'.format(len(tids), gene))
    
    c_dict = get_edge_colors()
    
    # groupby and count the number of isoforms that use each edge
    df = df.groupby(['edge_id', 'edge_type']).count().reset_index()
    df.rename({'tid':'n_isos'}, axis=1, inplace=True)
    
    sns.set_context('paper', font_scale=1.8)
    
    for e_type in df.edge_type.unique():
        temp = df.loc[df.edge_type == e_type]
        ax = sns.displot(data=temp, x='n_isos',
                         kind='hist',
                         color=c_dict[e_type], binwidth=5,
                         linewidth=0)

        xlabel = '# $\it{}$ isoforms'.format(gene)
        ylabel = '# {}s'.format(e_type)

        _ = ax.set(xlabel=xlabel, ylabel=ylabel)
        fname = '{}_{}_isos_per_{}.png'.format(opref, gene, e_type)
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
        
def plot_n_libs_v_avg_isos(df,
                           color='blue',
                           min_tpm=1,
                           gene_subset='polya',
                           sample='all', 
                           nov=['Known', 'NIC', 'NNC'],
                           opref='figures/'):
    """
    Plot a scatterplot with a regression line for the average
        number of isos / gene vs. # libraries / sample
        
    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        filt_df (pandas DataFrame): TALON abundance, filtered
        color (str): Color to plot plot in
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library", 
            used to groupby datasets displayed
        nov (str): Novelty category of 
            isoforms to consider
        opref (str): Output prefix to save figure
    """
    
    # get number of libraries per sample
    n_libs = get_n_libs_per_sample()    
    
    # get avg isos
    df = get_isos_per_gene(df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby='sample', 
                       nov=nov)
    df = df.mean().to_frame().rename({0: 'avg_isos'}, axis=1)
    df.reset_index(inplace=True)
    df.rename({'index': 'biosample'}, axis=1, inplace=True)
    
    # merge with number of libraries
    df = df.merge(n_libs, how='left', on='biosample')
    
    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='n_libraries', y='avg_isos', color='b')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Number of libraries / sample'
    ylabel = 'Average # isoforms / gene'

    # get coeffs of linear fit
    c1 = 'n_libraries'
    c2 = 'avg_isos'
    slope, intercept, r_value, p_value, std_err = stats.linregress(df[c1],df[c2])
    lines = mpl.lines.Line2D([0], [0])
    label = 'm={0:.1f}'.format(slope)

    r2, pval_spear = stats.spearmanr(df[c1],df[c2])

    print('Slope of correlation: {}'.format(slope))
    print('R of correlation: {}'.format(r_value))
    print('R2 of correlation: {}'.format(r2))

    sns.regplot(data=df, x=c1, y=c2,
                scatter=False, ax=ax, 
                color='b')
    sns.regplot(data=df, x=c1, y=c2,
        scatter=False, ax=ax, ci=0, color='b',
        line_kws={'linestyle':'-',
                  'label':"m={0:.1f}".format(slope)})

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    
    fname = '{}_libs_v_avg_isos.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
def plot_n_isos_gene_sample_dist():
    """
    Plot the distribution of # isoforms detected / gene / sample
    """
        
def plot_n_reads_v_avg_isos(df,
                            filt_df,
                            color='blue',
                            min_tpm=1,
                            gene_subset='polya',
                            sample='all', 
                            groupby='sample',
                            nov=['Known', 'NIC', 'NNC'],
                            opref='figures/'):
    """
    Plot a scatterplot with a regression line for the average
        number of isos / gene vs. # reads / sample
        
    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        filt_df (pandas DataFrame): TALON abundance, filtered
        color (str): Color to plot plot in
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library", 
            used to groupby datasets displayed
        nov (str): Novelty category of 
            isoforms to consider
        opref (str): Output prefix to save figure
    """
    
    # get number of reads from unfiltered data
    reads = get_reads_per_sample(df, groupby='sample')
    
    # get avg isos
    df = get_isos_per_gene(filt_df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby=groupby, 
                       nov=nov)
    df = df.mean().to_frame().rename({0: 'avg_isos'}, axis=1)
    df.reset_index(inplace=True)
    df.rename({'index': 'biosample'}, axis=1, inplace=True)
    
    # merge with read depth
    df = df.merge(reads, how='left', on='biosample')
    df.head()
    
    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='n_reads', y='avg_isos', color='b')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Reads per cell line / tissue'
    ylabel = 'Average # isoforms / gene'

    # get coeffs of linear fit
    c1 = 'n_reads'
    c2 = 'avg_isos'
    slope, intercept, r_value, p_value, std_err = stats.linregress(df[c1],df[c2])
    lines = mpl.lines.Line2D([0], [0])
    label = 'm={0:.1f}'.format(slope)
    
    r2, pval_spear = stats.spearmanr(df[c1],df[c2])

    print('Slope of correlation: {}'.format(slope))
    print('R of correlation: {}'.format(r_value))
    print('R2 of correlation: {}'.format(r2))

    sns.regplot(data=df, x=c1, y=c2,
                scatter=False, ax=ax, 
                color='b')
    sns.regplot(data=df, x=c1, y=c2,
        scatter=False, ax=ax, ci=0, color='b',
        line_kws={'linestyle':'-',
                  'label':"m={0:.1f}".format(slope)})

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    
    fname = '{}_reads_v_avg_isos.png'.format(opref)
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
                           sample=sample,
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
    
    # finally, calculate the number of biosamples / libraries these 
    # genes or transcripts are expressed >= min TPM
    df = df.transpose()
    df['n_samples'] = df.astype(int).sum(axis=1)
    
    # and make a beautiful plot
    sns.set_context('paper', font_scale=2)
    
    c_dict, order = get_talon_nov_colors()
    color = c_dict[nov]
    ax = sns.displot(data=df, x='n_samples', kind='hist',
                 color=color, binwidth=1, linewidth=0)

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
    
    
def plot_gene_tpm_v_n_isos(df, filt_df,
                           min_tpm=1,
                           groupby='sample', 
                           gene_subset='polya', 
                           nov=['Known', 'NIC', 'NNC'],
                           opref='figures/human'):
    """
    Plot a scatterplot of gene tpm / library or sample vs.
    number of isoforms per gene / library or sample
    
    Parameters:
        df (pandas DataFrame): TALON abundance file
        filt_df (pandas DataFrame): filtered talon abundance file
    """
    
    tpm_df, _ = get_tpm_table(df, 
                       how='gene', 
                       min_tpm=min_tpm,
                       groupby=groupby,
                       gene_subset=gene_subset)
    
    iso_df = get_isos_per_gene(filt_df, 
                               min_tpm=min_tpm,
                               gene_subset=gene_subset, 
                               groupby=groupby, 
                               nov=nov)
    tpm_df = tpm_df.melt(ignore_index=False, value_name='tpm')
    iso_df = iso_df.melt(ignore_index=False, value_name='n_iso', var_name='biosample').fillna(0)
    
    tpm_df.index.name = 'annot_gene_id'
    tpm_df.reset_index(inplace=True)
    iso_df.reset_index(inplace=True)

    df = tpm_df.merge(iso_df, how='outer', on=['annot_gene_id', 'biosample'])
    
    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='n_iso', y='tpm')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = '# isoforms / gene / sample'
    ylabel = 'TPM / gene / sample'

    # _ = ax.set(xlabel=xlabel, ylabel=ylabel, yscale='log')
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)


    fname = '{}_isos_v_tpm.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    # add gname
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname']]
    df = df.merge(gene_df, how='left', left_on='annot_gene_id', right_on='gid')
    
    return df
    
def plot_det_vs_gencode_isos(df,
                         min_tpm=1,
                         gene_subset='polya',
                         nov=['Known', 'NIC', 'NNC'],
                         label_genes=None,
                         opref='figures/',
                         ylim=None, 
                         xlim=None):
    """
    Plot a scatterplot of the the total number of isoforms detected
    vs the number of isoforms in gencode
    """
    
    # detected isoforms
    det_df = get_isos_per_gene(df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       groupby='all', 
                       nov=nov)
    det_df.rename({'all': 'n_isos_det'}, axis=1, inplace=True)
    
    # annotated isoforms
    gc_df = get_n_gencode_isos(subset='polya')
    gc_df = gc_df[['gid', 'n_isos_gencode']]
    
    df = det_df.merge(gc_df, how='left', left_index=True, right_on='gid')    
    
    # add gene name 
    gene_df, _, _ = get_gtf_info(how='gene', subset='polya')
    gene_df = gene_df[['gid', 'gname']]
    df = df.merge(gene_df, how='left', on='gid')
    
    # plot the figure
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(data=df, x='n_isos_det', y='n_isos_gencode')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    xlabel = 'Total # isoforms / gene'
    ylabel = '# isoforms / gene in GENCODE'
    _ = ax.set(xlabel=xlabel, ylabel=ylabel, xscale='log', yscale='log')
        
    # set x and y lims if provided
    if xlim:
        xlim = (0, xlim)
        ax.set(xlim=xlim)
    if ylim:
        ylim = (0, ylim)
        ax.set(ylim=ylim)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # annotate genes that are kinda interesting
    if label_genes:
        xlim = ax.get_xlim()[1]
        ylim = ax.get_ylim()[1]
        for g in label_genes:
            if g in df.gname.tolist():
                x = df.loc[df.gname == g, 'n_isos_det'].values[0]+(2/75)*xlim
                y = df.loc[df.gname == g, 'n_isos_gencode'].values[0]-(1.5/75)*ylim
                plt.annotate(g, (x,y), fontsize='small', fontstyle='italic')


    fname = '{}_total_v_gencode_isos_per_gene.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    return df
    
def plot_max_vs_all_isos(df,
                         min_tpm=1,
                         gene_subset='polya',
                         groupby='sample', 
                         nov=['Known', 'NIC', 'NNC'],
                         label_genes=None,
                         opref='figures/',
                         ylim=None, 
                         xlim=None):
    """
    Plot a scatterplot of the maximum number of isoforms detected 
    per sample or library vs. the total number of isoforms detected
    """
    df_copy = df.copy(deep=True)
    
    # get maximum number of detected isoforms per library or sample
    max_isos = get_isos_per_gene(df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       groupby=groupby, 
                       nov=nov)
    max_isos = max_isos.max(axis=1).to_frame()
    max_isos.rename({0: 'max_isos'}, axis=1, inplace=True)
    
    # get total number of detected isoforms overall
    total = get_isos_per_gene(df_copy,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       groupby='all', 
                       nov=nov)
    total.rename({'all': 'total_isos'}, axis=1, inplace=True)
                 
    # merge 
    df = max_isos.merge(total, left_index=True, right_index=True)
    
    # add gene name
    gene_df, _, _ = get_gtf_info(how='gene')
    df = df.merge(gene_df, how='left', left_index=True, right_on='gid')
    
    # plot the figure
    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='total_isos', y='max_isos')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # set x and y lims if provided
    if xlim:
        xlim = (0, xlim)
        ax.set(xlim=xlim)
    if ylim:
        ylim = (0, ylim)
        ax.set(ylim=ylim)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # annotate genes that are kinda interesting
    if label_genes:
        xlim = ax.get_xlim()[1]
        ylim = ax.get_ylim()[1]
        for g in label_genes:
            if g in df.gnames.tolist():
                x = df.loc[df.gname == g, 'total_isos'].values[0]+(2/75)*xlim
                y = df.loc[df.gname == g, 'max_isos'].values[0]-(1.5/75)*ylim
                plt.annotate(g, (x,y), fontsize='small', fontstyle='italic')

    xlabel = 'Total # isoforms / gene'
    ylabel = 'Max. # isoforms / gene in one sample'
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}_max_v_all_isos_per_gene.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    return df

def plot_isos_per_gene_hist(df,
                            min_tpm=1,
                            gene_subset='polya', 
                            sample='all',
                            groupby='sample', 
                            nov=['Known'],
                            rm_1=False,
                            pseudocount=None,
                            opref='figures/'):
    """
    Plots dist. of # isos / gene across the different samples
    
    Parameters:
        df (pandas DataFrame): TALON abundance file
        gene_subset (str): Choose from None or 'polya' or 'tf'
        sample (str): Choose from 'all', 'tissue', or 'cell_line'
        min_tpm (float): Min. TPM val for at least one library
        groupby (str): Choose from 'sample' or 'library'
        nov (list of str): Novelty categories to consider
        rm_1 (bool): Whether or not to remove 1-count gene / sample
        pseudocount (int): Add a pseudocount
        opref (str): Output file prefix 
    """
    
    df = get_isos_per_gene(df,
                           min_tpm=min_tpm,
                           gene_subset=gene_subset,
                           groupby=groupby,
                           sample=sample,
                           nov=nov)
    
    # get long form dataframe
    df = df.melt(ignore_index=False) 
    
    # remove 0-count gene / sample combos (as these
    # are genes that were not detected)
    df = df.loc[df.value != 0]
    
    # remove 1-count gene / sample combos too
    if rm_1:
        df = df.loc[df.value != 1]
    
    # add pseudocount if wanted
    if pseudocount:
        df.value = df.value + pseudocount
            
    sns.set_context('paper', font_scale=2)

    ax = sns.displot(data=df, x='value', kind='hist', binwidth=1, linewidth=0)
    xlabel = '# isoforms / gene / sample'
    ylabel = 'Number of genes'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel, yscale='log')
    for a in ax.axes.flat:    
        a.yaxis.set_major_formatter(ScalarFormatter())

    plt.savefig('{}_hist_isos_per_gene_per_sample.png'.format(opref), \
                dpi=300, bbox_inches='tight')
    return df, ax
    
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
    

def plot_transcript_novelty_per(df, 
                                gene='ELN',
                                min_tpm=1,
                                gene_subset='polya',
                                groupby='sample',
                                nov=['Known', 'NIC', 'NNC'], 
                                opref='figures/'):
    
    """
    Plot the # of detected isoforms / novelty category for a particular
    gene across samples or libraries
    """
    
    # get metadata about the transcripts that we need
    t_df = df.copy(deep=True)
    t_df = t_df[['annot_transcript_id', 'annot_gene_name',
                 'annot_gene_id', 'transcript_novelty']]  
    
    df = get_det_table(df,
              how='iso',
              min_tpm=min_tpm, 
              gene_subset=gene_subset,
              groupby=groupby,
              nov=nov)
    df = df.transpose()
    
    # isolate isoforms that belong to gene of interest
    df = df.merge(t_df, how='left', left_index=True, right_on='annot_transcript_id')
    df = df.loc[df.annot_gene_name == gene]
    
    # calc # isos detected per nov
    df.drop(['annot_gene_id', 'annot_gene_name'], axis=1, inplace=True)
    df = df.melt(id_vars=['annot_transcript_id', 'transcript_novelty'])
    df.rename({'variable': groupby, 
               'value': 'detected'}, axis=1, inplace=True)
    df = df.groupby([groupby, 'transcript_novelty', 'detected']).count().reset_index()
    df.rename({'annot_transcript_id':'counts'}, axis=1, inplace=True)
    df = df.loc[df.detected == True]
    
    c_dict, order = get_talon_nov_colors(nov)

    # actual plotting
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(4,6))
    g = sns.catplot(data=df, 
                x=groupby,
                hue='transcript_novelty',
                y='counts', kind='bar',
                palette=c_dict, hue_order=order,
                saturation=1)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]

    g.set_ylabels('$\it{}$ isoforms'.format(gene))
    g.set_xlabels('Transcript novelty')
    
    fname = '{}_{}_novelty_per_{}.png'.format(opref, gene, groupby)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
def plot_transcript_novelty_per_1(df, 
                                  gene='ELN',
                                  dataset='h9_chondro', 
                                  min_tpm=1,
                                  gene_subset='polya', 
                                  groupby='sample',
                                  nov=['Known', 'NIC', 'NNC'],
                                  opref='figures/'):
    
    # get metadata about the transcripts that we need
    t_df = df.copy(deep=True)
    t_df = t_df[['annot_transcript_id', 'annot_gene_name',
                 'annot_gene_id', 'transcript_novelty']]  
    
    df = get_det_table(df,
              how='iso',
              min_tpm=min_tpm, 
              gene_subset=gene_subset,
              groupby=groupby,
              nov=nov)
    
    # limit to only those detected in sample / library of interest
    df = df.loc[dataset].to_frame()
    df = df.loc[df[dataset] == True]
        
    # isolate isoforms that belong to gene of interest
    df = df.merge(t_df, how='left', left_index=True, right_on='annot_transcript_id')
    df = df.loc[df.annot_gene_name == gene]
    
    df = df[['annot_transcript_id', 'transcript_novelty']].groupby('transcript_novelty').count().reset_index()
    df.rename({'annot_transcript_id': 'counts'}, axis=1, inplace=True)
    c_dict, order = get_talon_nov_colors(nov)
    
    # actual plotting
    sns.set_context('paper', font_scale=1.8)
    plt.figure(figsize=(4,6))
    g = sns.catplot(data=df, x='transcript_novelty',
                y='counts', kind='bar',
                saturation=1,
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('$\it{}$ isoforms in {}'.format(gene, dataset))
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, df, 'counts')
    
    fname = '{}_{}_{}_novelty.png'.format(opref, gene, dataset)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    



    