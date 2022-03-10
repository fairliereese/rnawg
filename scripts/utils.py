import scanpy as sc
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import pyranges.pyranges as pr
import pyranges as pyranges
import scipy.stats as st

def rm_sirv_ercc(df):
    """From TALON ab file"""
    df = df.loc[~df.annot_gene_id.str.contains('SIRV')]
    df.loc[~df.annot_gene_id.str.contains('ERCC-')]
    return df

def get_dataset_cols():
    d = os.path.dirname(__file__)
    # fname = '{}/../lr_bulk/hr_to_biosample_type_back.tsv'.format(d)
    # print('Warning: using old version of hr_to_biosample_type. Is this ok?')
    fname = '{}/../lr_bulk/hr_to_biosample_type.tsv'.format(d)
    df = pd.read_csv(fname, sep='\t')
    datasets = df.hr.tolist()
    return datasets

def get_sample_datasets(sample=None):
    """
    Get the human-readable names of the datasets belonging
    to the input sample type.
    
    Parameters:
        sample (str): 'cell_line' or 'tissue'
        
    Returns:
        datasets (list of str): List of datasets belonging to that specific sample type
    """
    d = os.path.dirname(__file__)
    # fname = '{}/../lr_bulk/hr_to_biosample_type_back.tsv'.format(d)
    # print('Warning: using old hr_to_biosample_type, is this OK?')
    fname = '{}/../lr_bulk/hr_to_biosample_type.tsv'.format(d)
    df = pd.read_csv(fname, sep='\t')
    if sample == 'all':
        datasets = df.hr.tolist()
    elif sample:
        datasets = df.loc[df.biosample_type == sample, 'hr'].tolist()
    else: 
        datasets = df.hr.tolist()
        
    return datasets

def compute_detection(df, sample='cell_line',
                      how='iso', nov='Known'):
    
    df = rm_sirv_ercc(df)
    
    dataset_cols = get_sample_datasets(sample)
    
    if how == 'iso':
        df.set_index('annot_transcript_id', inplace=True)
        df = df.loc[df.transcript_novelty == nov]
        df = df[dataset_cols]

    # sum up counts across the same gene
    if how == 'gene':
        # only known genes
        df = df.loc[df.gene_novelty == 'Known']
        df = df[dataset_cols+['annot_gene_id']]
        df = df.groupby('annot_gene_id').sum()

    df = df.transpose()
    df.reset_index(inplace=True)
    df.rename({'index': 'dataset'}, axis=1, inplace=True)

    # get the celltype
    df['celltype'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]
    
    if sample == 'tissue':

        # add in the tissue metadata
        d = os.path.dirname(__file__)
        fname = '{}/../refs/tissue_metadata.csv'.format(d)
        tissue = pd.read_csv(fname)
        df = df.merge(tissue[['biosample', 'tissue']],
                        how='left', left_on='celltype',
                        right_on='biosample')
        df.drop('celltype', axis=1, inplace=True)
        df.rename({'tissue': 'celltype'}, axis=1, inplace=True)
        print('Found {} distinct tissues'.format(len(df.celltype.unique())))
    else:
        print('Found {} distinct cell lines'.format(len(df.celltype.unique())))   

    df.drop(['dataset'], axis=1, inplace=True)

    # sum over celltype
    df = df.groupby('celltype').sum()
    temp = df.copy(deep=True)

    max_df = get_rank_order(temp, 'max')
    temp = df.copy(deep=True)

    min_df = get_rank_order(temp, 'min')

    return max_df, min_df


def get_rank_order(df, how='max'):
    rank_order = ['n/a']
    rank_exp = [0]
    rank_cumulative = [0]

    n_celltypes = len(df.index.unique().tolist())

    while len(rank_order) < n_celltypes+1:

        # how many are expressed?
        df['n_expressed'] = df.gt(0).sum(axis=1)

        # which celltype expresses most?
        if how == 'max':
            celltype = df.n_expressed.idxmax()
        elif how == 'min':
            celltype = df.n_expressed.idxmin()

        n_exp = df.loc[celltype, 'n_expressed']

        if len(rank_order) == 1:
            rank_cumulative += [n_exp]
        else:
            rank_cumulative += [rank_cumulative[-1]+n_exp]

        rank_order += [celltype]
        rank_exp += [n_exp]

        # subset matrix by those that aren't expressed in the stashed celltype
        df.drop('n_expressed', axis=1, inplace=True)
        temp = df.loc[celltype].gt(0).to_frame()
        temp.rename({celltype: 'expressed'}, axis=1, inplace=True)
        remove_cols = temp.loc[temp.expressed == True].index.tolist()
        df.drop(remove_cols, axis=1, inplace=True)

        # also remove the celltype that was just analyzed
        df.drop(celltype, axis=0, inplace=True)

    temp = pd.DataFrame(data=rank_cumulative, columns=['n_cumulative'])
    temp['rank'] = temp.index+1
    temp['celltype'] = rank_order

    return temp

def get_tissue_metadata():
    """
    Get the biosample <--> higher level biosample mapping
    
    Returns:
        tissue (pandas DataFrame): DataFrame containing original
            biosample_term_name as well as higher level version
    """
    
    d = os.path.dirname(__file__)
    fname = '{}/../refs/tissue_metadata.csv'.format(d)
    tissue = pd.read_csv(fname)
    return tissue

def get_n_gencode_isos(subset=None):
    """
    Get a DataFrame of the number of annotated isos / gene in GENCODE
    """
    
    df, _, _ = get_gtf_info(how='iso',
                            subset=subset)    
    df = df[['gid', 'tid']]
    df = df.groupby('gid').count().reset_index()
    df.rename({'tid': 'n_isos_gencode'}, axis=1, inplace=True)
    df.sort_values(by='n_isos_gencode', ascending=False, inplace=True)
    gene_df, _, _ = get_gtf_info(how='gene', subset=subset)
    df = df.merge(gene_df, how='left', on='gid')
       
    return df

def add_tss_ic_tes(sg):
    """
    Adds the intron chain, tss, tes, and first splice donor
    of each transcript to the t_df object. Also adds coords
    of tss and tes
    
    Parameters: 
        sg (swan_vis SwanGraph): SwanGraph with annotated and observed transcripts
        
    Returns
        df (pandas DataFrame): DF with start / end vertex / coord
            info, ic, and first splice donor vertex info 
    """
    df = sg.t_df.copy(deep=True)
    
    # add intron chains
    paths = df.path.values.tolist()
    paths = [tuple(path[1:-1]) for path in paths]
    df['intron_chain'] = paths

    # add tss
    paths = df.loc_path.values.tolist()
    tsss = [path[0] for path in paths]
    df['tss'] = tsss

    # add tes
    paths = df.loc_path.values.tolist()
    tess = [path[-1] for path in paths]
    df['tes'] = tess

    # add first splice donor
    paths = df.loc_path.values.tolist()
    first_sds = [path[1] for path in paths]
    df['first_sd'] = first_sds
    
    # add coordinates
    cols = ['tss', 'tes']
    for c in cols:
        # first, add tss / tes coords
        df = df.merge(sg.loc_df[['vertex_id', 'chrom', 'coord']],
                  how='left', left_on=c, right_index=True) 
        df.drop(['vertex_id'], axis=1, inplace=True)
        df.rename({'chrom': '{}_chrom'.format(c),
                  'coord': '{}_coord'.format(c)},
                  axis=1, inplace=True)
        
    return df

def count_tss_ic_tes(df, subset=None):
    """
    Count up unique tsss, ics, and tess for 
    a given subset of transcript ids
    
    Parameters:
        df (pandas DataFrame): t_df from get_ic_tss_tes
        subset (list of str): List of transcript ids
        
    Returns:
        counts (pandas DataFrame): df w/ an entry detailing
            how many tss, ics, tes there are for each gene
    """
    df = df.copy(deep=True)
    df = df.loc[df.tid.isin(subset)]
    
    # raw tss, ic, tes count
    cols = ['tss', 'intron_chain', 'tes']
    for i, col in enumerate(cols):
        if col in ['tss', 'tes']:
            id_col = '{}_cluster'.format(col)
        else:
            id_col = col
        temp = df[[id_col, 'gid']].groupby('gid').nunique()
        if i == 0:
            counts = temp
        else:
            counts = counts.merge(temp, left_index=True, right_index=True)
            
    # unique combinations of tss, ic, tes
    df['tss_ic_tes'] = df.tss_cluster.astype(str)+'_'+\
                       df.intron_chain.astype(str)+'_'+\
                       df.tes_cluster.astype(str)
    temp = df[['tss_ic_tes', 'gid']].groupby('gid').nunique()
    counts = counts.merge(temp, how='left', left_index=True, right_index=True)
    
    for col in counts.columns:
        if col == 'tss_cluster':
            counts.rename({col: 'tss'}, axis=1, inplace=True)
        elif col == 'tes_cluster':
            counts.rename({col: 'tes'}, axis=1, inplace=True)
            
    # compute splicing ratio
    counts['splicing_ratio'] = counts.intron_chain/((counts.tes+counts.tss)/2)
    
    return counts

def df_to_pyranges(ends, kind='tss'):
    
    # reformat column names if needed
    cols = ends.columns
    if 'Start' in cols and 'End' in cols and 'Chromosome' in cols:
        pass
    else:
        coord = '{}_coord'.format(kind)
        chrom = '{}_chrom'.format(kind)
        ends = ends[cols].copy(deep=True)
        ends.rename({coord: 'Start', 
                     chrom: 'Chromosome'},
                     axis=1, inplace=True)
        ends['End'] = ends.Start
        
    # turn into a pyranges object
    cols = ['gid', 'gname',
            'Start', 'End',
            'Chromosome', kind]
    if kind == 'tss':
        cols.append('first_sd')
    ends = ends[cols]
    ends.drop_duplicates(inplace=True)
    ends = pr.PyRanges(df=ends)
    
    return ends
        
def cluster_ends(ends,
                 slack,
                 cluster_start=1,
                 kind='tss'):
    """
    Cluster TSSs / TESs.
    
    Parameters:
        ends (pandas DataFrame): Slice of dataframe from add_tss_ic_tes
        slack (int): Allowable distance for merging ends
        cluster_start (int): # to start numbering clusters from
        kind (str): 'tss' or 'tes'
        
    Returns:
        reg (pandas DataFrame): DF describing found regions
        clust (pandas DataFrame): DF describing which region 
            each end observation corresponds to
    """
    
    ends = df_to_pyranges(ends, kind=kind)

    # get regions and region assignments
    cols = ['gid', 'gname']
    if kind == 'tss':
        cols.append('first_sd')
    
    # merge to get regions
    reg = ends.merge(strand=None, by=cols, slack=slack)
    reg = reg.as_df()
    reg['len'] = reg.End - reg.Start
    reg['Cluster'] = [i for i in range(cluster_start, len(reg.index)+cluster_start)]
    
    # cluster to get region assignment per end
    clust = ends.cluster(strand=None, by=cols, slack=slack)
    clust = clust.as_df()
    clust['Cluster_new'] = clust.Cluster+cluster_start-1
    clust.drop('Cluster', axis=1, inplace=True)
    clust.rename({'Cluster_new': 'Cluster'}, axis=1, inplace=True)
    
    return reg, clust

def compute_triplets(t_df,
                     df,
                     min_tpm=1,
                     groupby='sample',
                     bool_col=None):
    """
    Compute the triplets on the sample or library level
    
    Parameters:
        t_df (pandas DataFrame): t_df output from get_ic_tss_tes
        df (pandas DataFrame): Filtered TALON abundance or 90% set 
            dataframe 
        min_tpm (int): Min TPM to be considered detected
        groupby (str): Choose 'library' or 'sample'
        
    Returns:
        counts (pandas DataFrame): DF w/ n tss, ic, tes, and 
            splicing ratio for each gene in each sample / lib
    """
    
    # get table of which isoforms are detected in 
    # which samples / libraries, or which isoforms 
    # are part of the 90% set / sample
    if 'transcript_novelty' in df.columns:
        df = get_det_table(df,
                           how='iso',
                           min_tpm=min_tpm,
                           groupby=groupby,
                           nov=['Known', 'NIC', 'NNC'])
    # otherwise, expect 90% set file format and coerce into 
    # boolean 90% set detection table format
    elif 'pi' in df.columns:
        df = df[['tid', 'biosample']]
        df['in_90_set'] = True
        df = df.pivot(index='biosample', columns='tid', values='in_90_set').fillna(value=False)
        df.columns.name = ''

    # loop through samples / libraries compute triplet
    # for detected transcripts in each sample
    counts = pd.DataFrame()
    for ind, entry in df.iterrows():
        entry = entry.to_frame()
        tids = entry.loc[entry[ind] == True].index.tolist()
        temp = count_tss_ic_tes(t_df, subset=tids)
        temp['source'] = ind
        if bool_col:
            temp[bool_col] = True
        counts = pd.concat([counts, temp])
    
    # get gene info and add
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname', 'biotype',
                  'biotype_category', 'tf']]
    gene_df.drop_duplicates(inplace=True)
    counts = counts.merge(gene_df, how='left', left_index=True, right_on='gid')
    
    return counts
    
def get_ic_tss_tes(sg,
                   df,
                   min_tpm=1,
                   gene_subset='polya',
                   annot_slack=200,
                   novel_slack=100,
                   verbose=False):
    """
    Extract information about annotaed and observed tss, tes, 
    and intron chain usage from a SwanGraph t_df. 
    
    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with both annotation
            and observed transcript data added
        df (pandas DataFrame): Filtered TALON abundance file
        min_tpm (int): Min TPM to consider a transcript detected 
            and therefore to include ends in end calling
        gene_subset (str): Choose from 'polya', 'tf'
        annot_slack (int): Distance b/w which to merge annotated ends
        novel_slack (int): Distance b/w which to merge observed ends
        verbose (bool): Whether or not to print output
        
    Returns:
        all_df (pandas DataFrame): sg.t_df modified to include 
            information about intron chain, tss, and tes
        regions (dict of pandas DataFrames): Indexed by
            'tss' and 'tes'. Bed regions for each end cluster
            as annotated in all_df
        counts (pandas DataFrame): DF of counts for intron 
            chains, TSSs, TESs, and unique combinations of the three
            for annotated, observed, and both
    """
    
    all_df = add_tss_ic_tes(sg)
    
    # limit to those annotated or in list of detected tids that we allow
    # only ever allow known, nic, nnc
    _, inds = get_tpm_table(df,
                             how='iso',
                             min_tpm=min_tpm,
                             gene_subset=gene_subset,
                             nov=['Known', 'NIC', 'NNC'])
    novel_tids = inds
    
    if type(novel_tids) == list:
        all_df = all_df.loc[(all_df.annotation == True)|(all_df.tid.isin(novel_tids))]

    end_types = ['tss', 'tes']
    end_regions = dict()
    for c in end_types:

        if verbose:
            print()

        #### annotated transcripts ####
        
        t_df = all_df.loc[all_df.annotation == True].copy(deep=True)
        if verbose: 
            n = len(t_df.index)
            print('Finding {}s for {} annotated transcripts'.format(c, n))

        reg, clust = cluster_ends(t_df,
                                  slack=annot_slack,
                                  cluster_start=1,
                                  kind=c)
        reg['annotation'] = True
        reg['source'] = 'GENCODE'
        clust['annotation'] = True
        if verbose:
            n = len(reg.index)
            print('Found {} annotated {} clusters'.format(n,c))

        #### novel transcripts ####    

        # assign ends from novel transcripts to annotated ends
        t_df = all_df.loc[(all_df.annotation == False)&(all_df.tid.isin(sg.adata.var.index.tolist()))]

        if verbose:
            n = len(t_df.index)
            print('Finding {}s for {} novel transcripts'.format(c, n))

        # case 1: ends from novel transcripts are w/i annotated regions
        ends = df_to_pyranges(t_df, kind=c)
        reg = pr.PyRanges(df=reg)
        ends = ends.join(reg, how='left',
                         slack=1,
                         strandedness=None, suffix='_annot')
        ends = ends.as_df()

        # limit to those w/ matching gid, first sd and add to cluster df
        if c == 'tss':
            inds = ends.loc[(ends.gid==ends.gid_annot)&(ends.first_sd==ends.first_sd_annot)].index.tolist()
        else:
            inds = ends.loc[ends.gid == ends.gid_annot].index.tolist()

        if verbose:
            n = len(inds)
            print('Found {} novel {}s that are already in the annotation'.format(n,c))        
        clust = pd.concat([clust, ends.loc[inds]])

        # case 2: ends from novel transcripts need to be clustered
        # on their own

        # remove duplicates that arise from imperfect merging
        ends['in_region'] = False
        ends.loc[inds, 'in_region'] = True
        ends.sort_values(by='in_region', inplace=True, ascending=False)
        cols = ['gid', c]
        if c == 'tss':
            cols.append('first_sd')
        ends.drop_duplicates(subset=cols, keep='first', inplace=True)
        inds = ends.loc[ends.in_region == True].index.tolist()

        # subset based on ends that are unsupported by ref regions
        inds = list(set(ends.index.tolist())-set(inds))
        t_df = ends.loc[inds]
        if verbose:
            n = len(t_df.index)
            print('Finding {}s for {} novel ends'.format(c,n))
        n = clust.Cluster.max()+1
        nov_reg, nov_clust = cluster_ends(t_df, 
                                          slack=novel_slack,
                                          cluster_start=n, 
                                          kind=c)
        nov_reg['annotation'] = False
        nov_reg['source'] = 'obs'
        nov_clust['annotation'] = False
        if verbose:
            n = len(nov_reg.index)
            print('Found {} novel {} clusters'.format(n,c))

        # check how many novel clusters fall into already 
        # annotated regions
        nov_reg = pr.PyRanges(df=nov_reg)
        temp = nov_reg.join(reg, how=None, strandedness=None, suffix='_annot')
        temp = temp.as_df()
        if verbose:
            if c == 'tss':
                temp = temp.loc[(temp.first_sd == temp.first_sd_annot)&(temp.gid == temp.gid_annot)]
            else:
                temp = temp.loc[temp.gid == temp.gid_annot]
            cols = ['gid']
            if c == 'tss':
                cols.append('first_sd')
            temp = temp.drop_duplicates(subset=cols)
            n = len(temp.index)
            print('{} new {} regions overlap annotated regions'.format(n,c))

        # finally, add novel regions to clust and reg dfs
        reg = reg.as_df()
        nov_reg = nov_reg.as_df()
        clust = pd.concat([clust, nov_clust])
        reg = pd.concat([reg, nov_reg])

        # some final formatting for these dfs
        cols = ['gid', 'gname', c,
                'Cluster', 'annotation']
        if c == 'tss':
            cols.append('first_sd')
        clust = clust[cols]
        clust.rename({'Cluster': '{}_cluster'.format(c),
                      'annotation': '{}_annotation'.format(c)},
                      axis=1, inplace=True)
        clust.drop_duplicates(inplace=True)
        end_regions[c] = reg

        # add cluster ids back into the original df
        cols = ['gid', 'gname', c]
        if c == 'tss':
            cols.append('first_sd')
        all_df = all_df.merge(clust, how='left', on=cols)
        
    # counts for all, annotated, and observed go into the same df
    # with a different source
    counts = pd.DataFrame()

    # annotated counts
    tids = all_df.loc[all_df.novelty == 'Known'].tid.tolist()
    annot_counts = count_tss_ic_tes(all_df, subset=tids)
    annot_counts['source'] = 'GENCODE'
    counts = pd.concat([counts, annot_counts])
    
    # annotated + novel counts
    tids = list(set(tids)|set(novel_tids))
    all_counts = count_tss_ic_tes(all_df, subset=tids)
    all_counts['source'] = 'all'
    counts = pd.concat([counts, all_counts])
    
    # observed counts
    tids = list(set(novel_tids)&set(sg.adata.var.index.tolist()))
    obs_counts = count_tss_ic_tes(all_df, subset=tids)
    obs_counts['source'] = 'obs'
    counts = pd.concat([counts, obs_counts])
    
    # get gene info and add
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname', 'biotype',
                  'biotype_category', 'tf']]
    gene_df.drop_duplicates(inplace=True)
    counts = counts.merge(gene_df, how='left', left_index=True, right_on='gid')
    
    t_df = all_df.copy(deep=True)
        
    return t_df, end_regions, counts

def get_gtf_info(how='gene',
                 subset=None):
    """
    Gets the info from the annotation about genes / transcripts
    
    Parameters:
        how (str): 'gene' or 'iso'
        subset (str): 'polya', 'tf', 'protein_coding' or None
        
    Returns:
        df (pandas DataFrame): DataFrame with info for gene / transcript
        biotype_counts (pandas DataFrame): DataFrame with the counts 
            per biotype reported in gencode
        biotype_cat_counts (pandas DataFrame): DataFrame with the counts
            per meta biotype reported in gencode
    """
    d = os.path.dirname(__file__)
    if how == 'gene':
        fname = '{}/../refs/gencode_v29_gene_metadata.tsv'.format(d)
    elif how == 'iso':
        fname = '{}/../refs/gencode_v29_transcript_metadata.tsv'.format(d)
            
    df = pd.read_csv(fname, sep='\t')
    
    if how == 'gene':
        id_col = 'gid'
    elif how == 'iso':
        id_col = 'tid'
    
    if subset == 'polya':
        polya_cats = ['protein_coding', 'lncRNA', 'pseudogene']
        df = df.loc[df.biotype_category.isin(polya_cats)]
    elif subset == 'protein_coding':
        df = df.loc[df.biotype_category == 'protein_coding']
    elif subset == 'pseudogene':
        df = df.loc[df.biotype_category == 'pseudogene']
    elif subset == 'tf':
        df = df.loc[df.tf == True]

    biotype_counts = df[[id_col, 'biotype']].groupby('biotype').count()
    biotype_counts.reset_index(inplace=True)
    biotype_counts.rename({id_col: 'gencode_counts'}, axis=1, inplace=True)
    
    biotype_cat_counts = df[[id_col, 'biotype_category']].groupby('biotype_category').count()
    biotype_cat_counts.reset_index(inplace=True)
    biotype_cat_counts.rename({id_col: 'gencode_counts'}, axis=1, inplace=True)
    
    return df, biotype_counts, biotype_cat_counts

def get_det_table(df,
                  how='gene',
                  min_tpm=1, 
                  gene_subset='polya',
                  sample='all',
                  groupby='library',
                  nov='Known'):
    """
    Get a dataframe of True / False whether or not a gene / isoform
    was detected in a specific library or sample
    
    Parameters:
        df (pandas DataFrame): TALON abundance
        how (str): Either "gene" or "iso"
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", 'library', or 'all' 
            used to groupby datasets displayed
        nov (list of str): Only used with how='iso', novelty categories of 
            isoforms to consider
        
    Returns: 
        df (pandas DataFrame): DataFrame with True / False entries 
            for each isoform / gene per library / sample
    """
    
    # calc TPM per library on desired samples
    df, tids = get_tpm_table(df,
                   sample=sample,
                   how=how,
                   nov=nov,
                   min_tpm=min_tpm,
                   gene_subset=gene_subset)
    
    df = df.transpose()
    df.index.name = 'dataset'
    df.reset_index(inplace=True)
    
    # set up df to groupby sample or library
    if groupby == 'sample':

        # add biosample name (ie without rep information)
        df['biosample'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]
        df.drop(['dataset'], axis=1, inplace=True)

        # record the highest TPM value per biosample
        tissue_df = get_tissue_metadata()
        tissue_df = tissue_df[['tissue', 'biosample']]

        df = df.merge(tissue_df, how='left', on='biosample')
        df.loc[df.tissue.isnull(), 'tissue'] = df.loc[df.tissue.isnull(), 'biosample']
        df.drop('biosample', axis=1, inplace=True)
        df.rename({'tissue': 'biosample'}, axis=1, inplace=True)

        print('Found {} total samples'.format(len(df.biosample.unique().tolist())))
        df = df.groupby('biosample').max()

    elif groupby == 'library':
        df.rename({'dataset': 'library'}, axis=1, inplace=True)
        print('Found {} total libraries'.format(len(df.library.unique().tolist())))
        df = df.groupby('library').max()
    
    elif groupby == 'all':
        df['dataset'] = 'all'
        df = df.groupby('dataset').max()
        
    df = (df >= min_tpm)
    return df

def get_reads_per_sample(df,
                         groupby='sample'):
    """
    Calculate the number of reads per sample
    
    Parameters:
        df (pandas DataFrame): Unfiltered TALON abundance file
        groupby (str): 'sample' or 'library'
    """
    
    # remove irrelevant columns
    dataset_cols = get_dataset_cols()
    cols = ['annot_transcript_id']+dataset_cols
    df = df[cols]
    df.set_index('annot_transcript_id', inplace=True)
    df = df.transpose()
    df.index.name = 'dataset'
    df.reset_index(inplace=True)
    df.columns.name = ''
    
    # calculate the number of reads per library
    datasets = df.dataset.tolist()
    df = df.sum(axis=1).to_frame()
    df['dataset'] = datasets
    df.rename({0: 'n_reads'}, axis=1, inplace=True)
    
    if groupby == 'sample':
        # add biosample name (ie without rep information)
        df['biosample'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]
        df.drop(['dataset'], axis=1, inplace=True)

        # record the highest TPM value per biosample
        tissue_df = get_tissue_metadata()
        tissue_df = tissue_df[['tissue', 'biosample']]

        df = df.merge(tissue_df, how='left', on='biosample')
        df.loc[df.tissue.isnull(), 'tissue'] = df.loc[df.tissue.isnull(), 'biosample']
        df.drop('biosample', axis=1, inplace=True)
        df.rename({'tissue': 'biosample'}, axis=1, inplace=True)

        print('Found {} total samples'.format(len(df.biosample.unique().tolist())))
        
        df = df.groupby('biosample').sum().reset_index()        

    return df

def get_n_libs_per_sample():
    """
    Calculate the number of libraries that makes up each sample
    
    Returns
        df (pandas DataFrame): DataFrame where one column is 
            the biosample and second column is # of libraries
    """
    
    datasets = get_dataset_cols()
    df = pd.DataFrame(data=datasets, columns=['dataset'])

    # add biosample name (ie without rep information)
    df['biosample'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]

    # record the highest TPM value per biosample
    tissue_df = get_tissue_metadata()
    tissue_df = tissue_df[['tissue', 'biosample']]

    df = df.merge(tissue_df, how='left', on='biosample')
    df.loc[df.tissue.isnull(), 'tissue'] = df.loc[df.tissue.isnull(), 'biosample']
    df.drop('biosample', axis=1, inplace=True)
    df.rename({'tissue': 'biosample'}, axis=1, inplace=True)

    df = df.groupby('biosample').count().reset_index()
    df.rename({'dataset': 'n_libraries'}, axis=1, inplace=True)
    
    return df

def get_isos_per_gene(df,
                      min_tpm=1,
                      gene_subset='polya',
                      sample='all', 
                      groupby='sample',
                      nov=['Known', 'NIC', 'NNC']):
    """
    Compute the number of isoforms expressed per gene per 
    sample or library
    
    Parameters:
        df (pandas DataFrame): TALON abundance
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library", 
            used to groupby datasets displayed
        nov (str): Novelty category of 
            isoforms to consider
        
    Returns: 
        df (pandas DataFrame): DataFrame detailing how many samples
            isoforms / gene / sample or library are detected
    """
    g_df = df.copy(deep=True)
    df = get_det_table(df, 
              how='iso',
              min_tpm=min_tpm,
              sample=sample,
              gene_subset=gene_subset,
              groupby=groupby,
              nov=nov)
    
    # merge with gene info
    df = df.transpose()
    g_df = g_df[['annot_gene_id', 'annot_transcript_id']]
    df = df.merge(g_df, how='left', left_index=True, right_on='annot_transcript_id')
    
    # count number of expressed isoforms / gene
    df = df.drop(['annot_transcript_id'], axis=1)
    df.set_index('annot_gene_id', inplace=True)
    df = df.astype(int)
    df.reset_index(inplace=True)
    df = df.groupby('annot_gene_id').sum()
    
    # convert 0s into nans
    df.replace(0, np.nan, inplace=True)
    
    return df

def get_gene_iso_det_table(df, filt_df, 
                           min_isos=2,
                           iso_nov=['Known', 'NIC', 'NNC'],
                           gene_nov=['Known'],
                           min_tpm=1,
                           gene_subset='polya', 
                           sample='all', 
                           groupby='sample'):
    
    """
    Compute a DataFrame which tells you whether genes
    contain more than a certain number of detected isoforms.
    
    
    
    """ 
    # get expressed genes
    gene_df= get_det_table(df,
                       how='gene',
                       nov=gene_nov,
                       min_tpm=min_tpm,
                       groupby=groupby,
                       gene_subset=gene_subset)
    gene_df = gene_df.transpose()
    gene_df.columns.name = ''  
    
    # get number of isoforms per gene 
    df = get_isos_per_gene(filt_df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby=groupby, 
                       nov=iso_nov)
    
    # >= n isoforms detected
    df = (df >= min_isos)  
    
    # left merge gene_df with df so we get all the expressed genes
    gene_df = gene_df.merge(df, how='left',
                            left_index=True, right_index=True,
                            suffixes=('_gene_det', '_2_iso_det'))
    
    # subset based on relevant columns and reformat
    gene_det_cols = [c for c in gene_df.columns if '_gene_det' in c]
    iso_det_cols = [c for c in gene_df.columns if '_2_iso_det' in c]

    iso_df = gene_df[iso_det_cols]
    gene_df = gene_df[gene_det_cols]

    iso_df.columns = [c.replace('_2_iso_det', '') for c in iso_df.columns]
    gene_df.columns = [c.replace('_gene_det', '') for c in gene_df.columns]
    
    # sort both dataframes by gene name
    gene_df.sort_index(inplace=True)
    df.sort_index(inplace=True)
    
    # make into ints
    iso_df.fillna(False, inplace=True)
    gene_df.fillna(False, inplace=True)

    iso_df = iso_df.astype(int).astype(str)
    gene_df = gene_df.astype(int).astype(str)
    
    # key:
    # 00: <2 isos detected, no gene detected
    # 01: <2 isos detected, gene detected
    # 10: >=2 isos detected, no gene detected (should be infrequent or never)
    # 11: >=2 isos detected, gene detected
    df = iso_df+gene_df
    df = df.transpose()
    
    return df

def get_tpm_table(df,
                    sample='all',
                    how='gene',
                    groupby='library',
                    nov=None,
                    min_tpm=None,
                    gene_subset=None,
                    save=False):
    """
    Parameters:
        df (pandas DataFrame): TALON abundance table
        sample (str): Choose from 'cell_line', 'tissue', or None
        how (str): Choose from 'gene' or 'iso'
        groupby (str): Choose from 'library' or 'sample'. Sample will avg.
        nov (list of str): List of accepted novelty types (w/ how='iso')
        min_tpm (float): Keep only genes / isos that have at least one
            TPM >= the value across the libraries
        gene_subset (str): Choose from 'polya' or None 
        save (bool): Whether or not to save the output matrix
        
    Returns:
        df (pandas DataFrame): TPMs for gene or isoforms in the requested
            samples above the input detection threshold.
        ids (list of str): List of str indexing the table 
    """
    print('Calculating {} TPM values'.format(how))
    
    if sample == 'cell_line' or sample == 'tissue':
        print('Subsetting for {} datasets'.format(sample))

    dataset_cols = get_sample_datasets(sample)
    df = rm_sirv_ercc(df)
    
    # merge with information about the gene
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'biotype_category', 'tf']]
    df = df.merge(gene_df, how='left', left_on='annot_gene_id', right_on='gid')
    
    # get indices that we'll need to subset on 
    if how == 'gene':
        id_col = 'annot_gene_id'
        nov_col = 'gene_novelty'
        nov = ['Known']
    elif how == 'iso':
        id_col = 'annot_transcript_id'
        nov_col = 'transcript_novelty'

    # filter on novelty 
    if nov: 
        print('Subsetting for novelty categories {}'.format(nov))
        nov_inds = df.loc[df[nov_col].isin(nov), id_col].tolist()
    else:
        nov_inds = df[id_col].tolist()

    # filter on gene subset
    if gene_subset:
        print('Subsetting for {} genes'.format(gene_subset))
        if gene_subset == 'polya':
            polya_cats = ['protein_coding', 'pseudogene', 'lncRNA']
            gene_inds = df.loc[df.biotype_category.isin(polya_cats), id_col].tolist()
        elif gene_subset == 'tf':
            gene_inds = df.loc[df.tf == True, id_col].tolist()
    else:
        gene_inds = df[id_col].tolist()

    # get intersection of both
    subset_inds = list(set(nov_inds)&set(gene_inds))

    # sum up counts across the same gene
    if how == 'gene':
        df = df[dataset_cols+[id_col]]
        df = df.groupby(id_col).sum().reset_index()

    # set index so that all values in df reflect
    # counts per transcript or gene
    df.set_index(id_col, inplace=True)

    # compute TPM
    tpm_cols = []
    for d in dataset_cols:
        tpm_col = '{}_tpm'.format(d)
        total_col = '{}_total'.format(d)
        df[total_col] = df[d].sum()
        df[tpm_col] = (df[d]*1000000)/df[total_col]
        tpm_cols.append(tpm_col)
    df = df[tpm_cols]

    # reformat column names
    df.columns = [c.rsplit('_', maxsplit=1)[0] for c in df.columns] 

    # enforce tpm threshold
    if min_tpm:
        print('Enforcing minimum TPM')
        print('Total # {}s detected: {}'.format(how, len(df.index)))
        df = df.loc[(df >= min_tpm).any(axis=1)]
        print('# {}s >= {} tpm: {}'.format(how, min_tpm, len(df.index)))

    # subset if necessary
    if gene_subset or nov:
        print('Applying gene type and novelty subset')
        df = df.loc[df.index.isin(subset_inds)]
        
    # average over biosample
    if groupby == 'sample':
        print('Averaging over biosample')
        df = df.transpose()
        df.reset_index(inplace=True)
        
        # add biosample name (ie without rep information)
        df['biosample'] = df['index'].str.rsplit('_', n=2, expand=True)[0]
        df.drop(['index'], axis=1, inplace=True)

        # record the avg TPM value per biosample
        tissue_df = get_tissue_metadata()
        tissue_df = tissue_df[['tissue', 'biosample']]

        df = df.merge(tissue_df, how='left', on='biosample')
        df.loc[df.tissue.isnull(), 'tissue'] = df.loc[df.tissue.isnull(), 'biosample']
        df.drop('biosample', axis=1, inplace=True)
        df.rename({'tissue': 'biosample'}, axis=1, inplace=True)

        print('Found {} total samples'.format(len(df.biosample.unique().tolist())))

        df = df.groupby('biosample').mean()
        df = df.transpose()
       
    print('Number of {}s reported: {}'.format(how, len(df.index)))

    if save:
        fname = '{}_{}_tpm.tsv'.format(sample, how)
        df.to_csv(fname, sep='\t')

    ids = df.index.tolist()    
        
    return df, ids

def compute_corr(df, how='gene', nov='Known', sample='cell_line'):

    dataset_cols = get_sample_datasets(sample)
    df = rm_sirv_ercc(df)

    if how == 'iso':
        df.set_index('annot_transcript_id', inplace=True)
        df = df.loc[df.transcript_novelty == nov]
        df = df[dataset_cols]

    # sum up counts across the same gene
    if how == 'gene':
        # only known genes
        df = df.loc[df.gene_novelty == 'Known']
        df = df[dataset_cols+['annot_gene_id']]
        df = df.groupby('annot_gene_id').sum()

    # sanity check
    print(len(df.index))

    # compute TPM
    tpm_cols = []
    for d in dataset_cols:
        tpm_col = '{}_tpm'.format(d)
        total_col = '{}_total'.format(d)
        df[total_col] = df[d].sum()
        df[tpm_col] = (df[d]*1000000)/df[total_col]
        tpm_cols.append(tpm_col)
    df = df[tpm_cols]

    # compute correlation between each set of datasets
    data = [[np.nan for i in range(len(df.columns))] for j in range(len(df.columns))]
    corrs = pd.DataFrame(data=data, index=df.columns, columns=df.columns)

    tested = []
    for d1 in df.columns.tolist():
        for d2 in df.columns.tolist():
            if [d1, d2] in tested:
                continue
            tested.append([d1, d2])
            tested.append([d2, d1])
            corr = st.pearsonr(df[d1].tolist(), df[d2].tolist())
            corrs.at[d1, d2] = corr[0]
            corrs.at[d2, d1] = corr[0]

    corrs.reset_index(inplace=True)

    if sample == 'tissue':

        # add in the tissue metadata
        d = os.path.dirname(__file__)
        fname = '{}/../refs/tissue_metadata.csv'.format(d)
        tissue = pd.read_csv(fname)
        corrs['celltype'] = corrs['index'].str.rsplit('_', n=2, expand=True)[0]
        corrs = corrs.merge(tissue[['biosample', 'tissue']],
                        how='left', left_on='celltype',
                        right_on='biosample')

        corrs.sort_values(by='tissue', inplace=True)
        corrs.drop(['tissue', 'celltype'], axis=1, inplace=True)
        corrs.set_index('index', inplace=True)
        corrs = corrs[corrs.index.tolist()]
    else:
        corrs.sort_values(by='index', inplace=True)
        corrs.set_index('index', inplace=True)
        corrs = corrs[corrs.index.tolist()]
    return corrs


def filter_cells(adata, min_umi,
                 max_umi,
                 max_mt,
                 min_genes,
                 depth, verbose=False):
    """
    """
    if verbose:
        n = len(adata.obs.loc[adata.obs.depth == depth].index)
        print('# cells for depth {}: {}'.format(depth, n))

    # either has to satisfy the cutoff or be at a different depth

    # > ### UMI filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.umi_count > min_umi)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
        n = len(adata.obs.loc[adata.obs.depth == depth].index)
        print('# cells for depth {} after removing cells with < {} UMI: {}'.format(depth, min_umi, n))

    # < ### UMI filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.umi_count < max_umi)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
       n = len(adata.obs.loc[adata.obs.depth == depth].index)
       print('# cells for depth {} after removing cells with > {} UMI: {}'.format(depth, max_umi, n))

    # > ### genes filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.n_genes_by_counts > min_genes)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
       n = len(adata.obs.loc[adata.obs.depth == depth].index)
       print('# cells for depth {} after removing cells with < {} genes: {}'.format(depth, min_genes, n))

    # < % MT filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.pct_counts_mt < max_mt)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
       n = len(adata.obs.loc[adata.obs.depth == depth].index)
       print('# cells for depth {} after removing cells with < {} % MT: {}'.format(depth, max_mt, n))

    return adata

def read_raw_data(meta, genes, mtx, depth):
    """
    """
    obs = pd.read_csv(meta)
    var = pd.read_csv(genes)
    adata = sc.read_mtx(mtx)
    X = adata.X
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    adata.obs['depth'] = depth
    adata.obs.set_index('cell_barcode', inplace=True)
    adata.var.set_index('gene_id', inplace=True)

    return adata


def load_raw_data(dataset):
    """
    """
    # different files for each of the tissues
    # TODO - can probably automate this a lil more
    if dataset == 'cortex':

        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/cortex/sr_splitseq/splitpipe/'

        # shallow
        meta = d+'cortex_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'cortex_12k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'cortex_12k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '12k'
        adata_shallow = read_raw_data(meta, genes, mtx, depth)

        # deep
        meta = d+'cortex_2k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'cortex_2k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'cortex_2k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        adata_deep = read_raw_data(meta, genes, mtx, depth)

        adata = adata_shallow.concatenate(adata_deep)

        # drop some trash
        adata.obs.drop(['batch', 'Unnamed: 0'], axis=1, inplace=True)
        adata.var.drop(['Unnamed: 0-0', 'Unnamed: 0-1'], axis=1, inplace=True)

    elif dataset == 'hippocampus':

        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/hippocampus/sr_splitseq/splitpipe/'

        # shallow
        meta = d+'hippocampus_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_12k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_12k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '12k'
        adata_shallow = read_raw_data(meta, genes, mtx, depth)

        # deep
        meta = d+'hippocampus_2k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_2k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_2k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        adata_deep = read_raw_data(meta, genes, mtx, depth)

        adata = adata_shallow.concatenate(adata_deep)

        # drop some trash
        adata.obs.drop(['batch', 'Unnamed: 0'], axis=1, inplace=True)
        adata.var.drop(['Unnamed: 0-0', 'Unnamed: 0-1'], axis=1, inplace=True)

    # add metadata
    # add additional metadata
    adata.obs['brain_region'] = dataset
    adata.obs['age'] = adata.obs['sample'].str.split('_', expand=True)[1]
    adata.obs['sex'] = adata.obs['sample'].str.split('_', expand=True)[2]
    adata.obs['rep'] = adata.obs['sample'].str.split('_', expand=True)[3]

    # calc some metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = adata.var.gene_name.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    return adata

def get_reads_per_bc(df):
    """
    Parameters:
        df (pandas DataFrame): read_annot df
    """
    df = df[['read_name', 'dataset']]
    temp = df.groupby('dataset').count().reset_index()
    temp.rename({'read_name': 'counts',
                 'dataset': 'barcode'}, axis=1, inplace=True)
    temp.sort_values(by='counts', ascending=False, inplace=True)
    return temp

def get_transcript_exp(df):
    """
    Parameters:
        df (pandas DataFrame): talon ab
    """
    # get gene level
    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                           'annot_transcript_id', 'annot_gene_name',
                           'annot_transcript_name', 'n_exons', 'length',
                           'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_columns ]
    id_col = ['annot_transcript_id', 'annot_transcript_name', \
              'annot_gene_id', 'annot_gene_name', 'transcript_novelty']

    # aggregate by gene
    t_df = df[id_col+dataset_cols]

    return t_df


def get_gene_exp(df, filter_novel=True):
    """
    Parameters:
        df (pandas DataFrame): talon ab
    """
    # get gene level
    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                           'annot_transcript_id', 'annot_gene_name',
                           'annot_transcript_name', 'n_exons', 'length',
                           'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_columns ]
    id_col = ['annot_gene_id', 'annot_gene_name']
    novelty_col = 'gene_novelty'

    if filter_novel:
        # only use known genes
        df = df.loc[df.gene_novelty == 'Known']

    # aggregate by gene
    gene_df = df[id_col+dataset_cols].groupby(id_col).sum()
    gene_df.reset_index(inplace=True)

    return gene_df

def get_bc1_matches():
    # from spclass.py - barcodes and their well/primer type identity
    d = os.path.dirname(__file__)
    bc_file = '{}/../refs/bc_8nt_v2.csv'.format(d)
    bc_df = pd.read_csv(bc_file, index_col=0, names=['bc'])
    bc_df['well'] = [i for i in range(0, 48)]+[i for i in range(0, 48)]
    bc_df['primer_type'] = ['dt' for i in range(0, 48)]+['randhex' for i in range(0, 48)]

    # pivot on well to get df that matches bcs with one another from the same well
    bc_df = bc_df.pivot(index='well', columns='primer_type', values='bc')
    bc_df = bc_df.rename_axis(None, axis=1).reset_index()
    bc_df.rename({'dt': 'bc1_dt', 'randhex': 'bc1_randhex'}, axis=1, inplace=True)

    return bc_df

def get_sample_metadata(samples):

    d = os.path.dirname(__file__)

    fname = '{}/../refs/age_metadata.tsv'.format(d)
    age = pd.read_csv(fname, sep='\t')

    fname = '{}/../refs/sex_metadata.tsv'.format(d)
    sex = pd.read_csv(fname, sep='\t')

    fname = '{}/../refs/tissue_metadata.tsv'.format(d)
    tissue = pd.read_csv(fname, sep='\t')

    samples[['tissue', 'age', 'sex', 'rep']] = samples['sample'].str.split('_', expand=True)

    samples = samples.merge(age, how='left', left_on='age', right_on='short')
    samples = samples.merge(tissue, how='left', left_on='tissue', right_on='short')

    samples = samples[['sample', 'tissue_desc', 'age_desc', 'sex', 'rep']]
    samples.rename({'tissue_desc': 'tissue', 'age_desc': 'age'},
                  axis=1, inplace=True)
    return samples


def get_illumina_metadata(dataset):

    # we want info about primer type as well
    bc_df = get_bc1_matches()

    # read in illumina bcs
    d = os.path.dirname(__file__)
    fname = '{}/../{}/sr_splitseq/scanpy/illumina_raw_metadata.csv'.format(d, dataset)
    ill_df = pd.read_csv(fname)

    # get 24nt barcode
    ill_df = ill_df.merge(bc_df, how='left', \
                 left_on='rnd1_well', right_on='well')

    ill_df['bc3'] = ill_df.cell_barcode.str.slice(start=0, stop=8)
    ill_df['bc2'] = ill_df.cell_barcode.str.slice(start=8, stop=16)
    ill_df['bc1'] = ill_df.bc1_dt
    ill_df['bc'] = ill_df.bc3+ill_df.bc2+ill_df.bc1

    # get metadata
    if ill_df['sample'].values[0].count('_') < 3:
        ill_df['sample_pref'] = ill_df['sample'].str.slice(0, -1)
        ill_df['sample_suff'] = ill_df['sample'].str.slice(-1, -2, -1)
        ill_df['sample'] = ill_df.sample_pref+'_'+ill_df.sample_suff
        ill_df.drop(['sample_pref', 'sample_suff'], axis=1, inplace=True)

    # port metadata into separate columns
    samples = pd.DataFrame(ill_df['sample'].unique())
    samples.columns = ['sample']
    samples[['tissue', 'age', 'sex', 'rep']] = samples['sample'].str.split('_', expand=True)
    samples = get_sample_metadata(samples)
    ill_df = ill_df.merge(samples, how='left', on='sample')
    ill_df.rename({'umi_count': 'sr_umi_count',
               'gene_count': 'sr_gene_count'},
               axis=1, inplace=True)

    # add annotation info
    if dataset == 'adrenal':
        fname = '{}/../{}/sr_splitseq/scanpy/illumina_processed_metadata.csv'.format(d, dataset)
        sr_df = pd.read_csv(fname)
        sr_df.rename({'barcode': 'bc',
               'cellType': 'sr_celltype',
               'seurat_clusters': 'sr_clusters'}, axis=1, inplace=True)
        sr_df.sr_clusters = sr_df.sr_clusters.astype(str)

    elif dataset == 'cortex':
        fname = '{}/../{}/sr_splitseq/scanpy/illumina_processed_metadata.tsv'.format(d, dataset)
        sr_df = pd.read_csv(fname, sep='\t')
        sr_df.rename({'cell_barcode': 'bc'}, axis=1, inplace=True)
        sr_df[['cell_barcode', 'rnd1_well']] = sr_df.bc.str.split('_', expand=True)
        sr_df.rnd1_well = sr_df.rnd1_well.astype(int)

        # get 24nt barcode
        sr_df = sr_df.merge(bc_df, how='left', \
                     left_on='rnd1_well', right_on='well')

        sr_df['bc3'] = sr_df.cell_barcode.str.slice(start=0, stop=8)
        sr_df['bc2'] = sr_df.cell_barcode.str.slice(start=8, stop=16)
        sr_df['bc1'] = sr_df.bc1_dt
        sr_df['bc'] = sr_df.bc3+sr_df.bc2+sr_df.bc1

        sr_df = sr_df[['bc', 'celltype', 'leiden']]
        sr_df.rename({'celltype': 'sr_celltype',
                      'leiden': 'sr_clusters'},
                      axis=1, inplace=True)
        sr_df.sr_clusters = sr_df.sr_clusters.astype(str)

    if dataset != 'hippocampus':
        ill_df = ill_df.merge(sr_df, how='left', on='bc')

    return ill_df

def make_adata(df, ill_df=None, verbose=False, how='gene'):

    if how == 'gene':
        var = df[['annot_gene_id', 'annot_gene_name']]
        df.drop(['annot_gene_name'], axis=1, inplace=True)
        df.set_index('annot_gene_id', inplace=True)
    elif how == 'transcript':
        var = df[['annot_transcript_id', 'annot_transcript_name', \
                'annot_gene_id', 'annot_gene_name', 'transcript_novelty']]
        df.drop(['annot_transcript_name', 'annot_gene_id', \
                 'annot_gene_name', 'transcript_novelty'], axis=1, inplace=True)
        df.set_index('annot_transcript_id', inplace=True)

    df = df.transpose()
    df.index.name = 'bc'
    X = df.values
    df.reset_index(inplace=True)
    obs = df.bc.to_frame()
    obs = df.bc.to_frame()
    obs['bc3_long'] = obs['bc'].str.slice(0,8)
    obs['bc2_long'] = obs['bc'].str.slice(8,16)
    obs['bc1_long'] = obs['bc'].str.slice(16,-1)

    if verbose:
        print('Found {} unique bc3s'.format(len(obs.bc3_long.unique())))
        print('Found {} unique bc2s'.format(len(obs.bc2_long.unique())))
        print('Found {} unique bc1s'.format(len(obs.bc1_long.unique())))


    if ill_df is not None:
        obs = obs.merge(ill_df, how='left', on='bc')
    #     obs.set_index('bc', inplace=True)
#         obs.seurat_clusters = obs.seurat_clusters.astype('category')
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.obs.set_index('bc', inplace=True)

    if how == 'gene':
        adata.var.set_index('annot_gene_id', inplace=True)
    elif how == 'transcript':
        adata.var.set_index('annot_transcript_id', inplace=True)
#     adata.raw = adata

    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var.annot_gene_name.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    return adata


def write_seurat_tables(adata, opref):
    obs = adata.obs
    x = adata.X

    x_coords = [n[0] for n in adata.obsm['X_umap']]
    y_coords = [n[1] for n in adata.obsm['X_umap']]

    obs['umap_x'] = x_coords
    obs['umap_y'] = y_coords

    row_labels = obs.index.tolist()
    col_labels = adata.var.index.tolist()

    x = pd.DataFrame(data=x, index=row_labels, columns=col_labels)

    obs.to_csv('{}_obs.tsv'.format(opref), sep='\t')
    x.to_csv('{}_x.tsv'.format(opref), sep='\t')

def write_swan_input_tables(adata, opref):

    # dump to formats that swan can handle
    meta = adata.obs.copy(deep=True)
    meta.reset_index(inplace=True)
    meta.rename({'bc':'dataset'}, axis=1, inplace=True)
    fname = '{}_swan_metadata.tsv'.format(opref)
    meta.to_csv(fname, sep='\t', index=False)

    index = adata.obs.index.tolist()
    cols = adata.var.index.tolist()
    data = adata.raw.X

    df = pd.DataFrame(data=data, index=index, columns=cols)
    df.reset_index(inplace=True)
    df.rename({'index': 'dataset'}, axis=1, inplace=True)
    df.set_index('dataset', inplace=True)
    df = df.transpose()
    df.index.name = 'transcript_id'
    df.reset_index(inplace=True)
    fname = '{}_swan_abundance.tsv'.format(opref)
    df.to_csv(fname, sep='\t', index=False)

def write_bc_leiden(adata, opref):
    obs = adata.obs['leiden']
    obs.to_csv('{}_leiden.tsv'.format(opref), sep='\t')
