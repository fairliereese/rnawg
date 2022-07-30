import pdb
import scanpy as sc
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
import pyranges as pr
import pyranges as pyranges
import cerberus
import scipy.stats as st

def get_biotype_map():
    """
    Get a dictionary mapping each gene type to a more general biotype
    """
    map = {'protein_coding': ['protein_coding'],
           'lncRNA': ['lincRNA',
                      'lncRNA',
                      'processed_transcript',
                      'sense_intronic',
                      '3prime_overlapping_ncRNA',
                      'bidirectional_promoter_lncRNA',
                      'sense_overlapping',
                      'non_coding',
                      'macro_lncRNA',
                      'antisense'],
           'pseudogene': ['unprocessed_pseudogene',
                          'translated_unprocessed_pseudogene',
                          'transcribed_unprocessed_pseudogene',
                          'processed_pseudogene',
                          'transcribed_processed_pseudogene',
                          'transcribed_unitary_pseudogene',
                          'unitary_pseudogene',
                          'polymorphic_pseudogene',
                          'pseudogene',
                          'translated_processed_pseudogene'],
           'miRNA': ['miRNA'],
           'other': ['snRNA', 'vault_RNA',
                     'misc_RNA', 'TEC',
                     'snoRNA', 'scaRNA',
                     'rRNA_pseudogene', 'rRNA',
                     'IG_V_pseudogene',
                     'scRNA', 'IG_V_gene',
                     'IG_C_gene', 'IG_J_gene',
                     'sRNA', 'ribozyme',
                     'vaultRNA', 'TR_C_gene',
                     'TR_J_gene', 'TR_V_gene',
                     'TR_V_pseudogene', 'TR_D_gene',
                     'IG_C_pseudogene', 'IG_D_gene',
                     'IG_pseudogene', 'Mt_tRNA',
                     'Mt_rRNA', 'TR_J_pseudogene',
                     'IG_J_pseudogene']}
    return map

def get_gene_info(gtf, o):
    df = pr.read_gtf(gtf, as_df=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # only gene
    df = df.loc[df.Feature == 'gene'].copy(deep=True)

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    # gene length
    df['length'] = (df.Start-df.End).abs()

    # add TF info
    df['tf'] = False
    d = os.path.dirname(__file__)
    tf_df = pd.read_csv('{}/../refs/biomart_tf_gids.tsv'.format(d), sep='\t')
    tf_gids = tf_df['Gene stable ID'].unique().tolist()
    df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
    df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
    df.drop('gid_stable', axis=1, inplace=True)

    # and save
    df = df[['gid', 'gname', 'length', 'biotype', 'biotype_category', 'tf']]
    df.to_csv(o, sep='\t', index=False)

def get_transcript_info(gtf, o):
    """
    Get a file with relevant information about each transcript in a gtf

    Parameters:
        gtf (str): Path to gtf
        o (str): Output file name
    """

    df = pr.read_gtf(gtf, as_df=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # only exons
    df = df.loc[df.Feature == 'exon'].copy(deep=True)

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    df['exon_len'] = (df.Start-df.End).abs()+1

    df = df[['gid', 'gname', 'tid', 'exon_len', 'biotype', 'biotype_category']]
    df_copy = df[['gid', 'gname', 'tid', 'biotype', 'biotype_category']].copy(deep=True)
    df_copy = df_copy.drop_duplicates(keep='first')

    df = df.groupby('tid').sum().reset_index()
    df.rename({'exon_len': 't_len'}, axis=1, inplace=True)
    df = df.merge(df_copy, on='tid', how='left')

    # add TF info
    df['tf'] = False
    d = os.path.dirname(__file__)
    tf_df = pd.read_csv('{}/../refs/biomart_tf_gids.tsv'.format(d), sep='\t')
    tf_gids = tf_df['Gene stable ID'].unique().tolist()
    df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
    df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
    df.drop('gid_stable', axis=1, inplace=True)

    # and save
    df.to_csv(o, sep='\t', index=False)

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
        sample (str): 'cell_line', 'tissue', 'mouse_match',
            'ljungman'

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
    elif sample == 'tissue' or sample == 'cell_line':
        datasets = df.loc[df.biosample_type == sample, 'hr'].tolist()
    elif sample == 'mouse_match':
        tissues = ['adrenal gland', 'brain', 'muscle', 'heart']
        tissue_df = get_tissue_metadata()
        df['biosample'] = df.hr.str.rsplit('_', n=2, expand=True)[0]
        df = df.merge(tissue_df, how='left', on='biosample')
        datasets = df.loc[df.tissue.isin(tissues), 'hr'].tolist()
    elif sample == 'heart':
        tissues = ['heart']
        tissue_df = get_tissue_metadata()
        df['biosample'] = df.hr.str.rsplit('_', n=2, expand=True)[0]
        df = df.merge(tissue_df, how='left', on='biosample')
        datasets = df.loc[df.tissue.isin(tissues), 'hr'].tolist()
    elif sample == 'ljungman':
        fname = '{}/../bru/ljungman_datasets.tsv'.format(d)
        sample_df = pd.read_csv(fname, sep='\t', header=None, names=['dataset'])
        df = df.merge(sample_df, how='inner', left_on='hr', right_on='dataset')
        datasets = df.hr.tolist()
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

def get_ad_metadata():
    """
    Get the human-readable dataset name <--> AD status mapping

    Returns:
        ad_df (pandas DataFrame): DataFrame containing dataset id
            and AD status
    """
    d = os.path.dirname(__file__)
    fname = '{}/../refs/ad_metadata.tsv'.format(d)
    ad_df = pd.read_csv(fname, sep='\t')
    ad_df = ad_df[['file_id', 'health_status']]
    fname = '{}/../lr_bulk/file_to_hr.tsv'.format(d)
    hr_df = pd.read_csv(fname, sep='\t', header=None, names=['file_id', 'hr'])

    ad_df = ad_df.merge(hr_df, on='file_id')
    return ad_df

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

def get_n_gencode_isos(subset=None, ver='v29'):
    """
    Get a DataFrame of the number of annotated isos / gene in GENCODE
    """

    df, _, _ = get_gtf_info(how='iso',
                            subset=subset, 
                            ver=ver,
                            add_stable_gid=True)
    df = df[['gid', 'tid']]
    df = df.groupby('gid').count().reset_index()
    df.rename({'tid': 'n_isos_gencode'}, axis=1, inplace=True)
    df.sort_values(by='n_isos_gencode', ascending=False, inplace=True)
    gene_df, _, _ = get_gtf_info(how='gene',
                                 subset=subset,
                                 ver=ver,
                                 add_stable_gid=True)
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

    if subset:
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
    ends = pr.pyranges.PyRanges(df=ends)

    return ends

def cluster_ends(ends,
                 slack,
                 cluster_start=1,
                 kind='tss',
                 first_sd=True):
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
    if kind == 'tss' and first_sd:
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

def get_subset_triplets(t_df,
                     df=None,
                     min_tpm=1,
                     sample='all',
                     groupby='sample',
                     source_name=None):
    """
    Compute the triplets on a subset of samples / libraries / isoforms

    Parameters:
        t_df (pandas DataFrame): t_df output from get_ref_triplets
        df (pandas DataFrame): Filtered TALON abundance or 90% set
            dataframe. Do not include if you want to compute for
            GENCODE + obs
        min_tpm (int): Min TPM to be considered detected
        sample (str): Choose 'cell_line', 'tissue', 'mouse_match'
        groupby (str): Choose 'library', 'sample', or 'all'

    Returns:
        counts (pandas DataFrame): DF w/ n tss, ic, tes, and
            splicing ratio for each gene in each sample / lib
    """

    # get table of which isoforms are detected in
    # which samples / libraries, or which isoforms
    # are part of the 90% set / sample
    if isinstance(df, pd.DataFrame):
        if 'transcript_novelty' in df.columns:
            df = get_det_table(df,
                               how='iso',
                               min_tpm=min_tpm,
                               sample=sample,
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
            counts = pd.concat([counts, temp])

    # if we don't have a df, we want to compute triplets for
    # all annotated + observed data
    else:
        counts = count_tss_ic_tes(t_df)

    # get gene info and add
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname', 'biotype',
                  'biotype_category', 'tf']]
    gene_df.drop_duplicates(inplace=True)
    counts = counts.merge(gene_df, how='left', left_index=True, right_on='gid')

    # if we were given a source name, replace all source names
    # with the given one
    if source_name:
        counts['source'] = source_name

    return counts

def get_ref_triplets(sg,
                     df,
                     first_sd=True,
                     annot_slack=200,
                     novel_slack=100,
                     verbose=False):
    """
    Get reference set of triplets for all the detected complete transcripts
    in our dataset as well as for GENCODE.

    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with both annotation
            and observed transcript data added
        df (pandas DataFrame): Filtered TALON abundance file
        annot_slack (int): Distance b/w which to merge annotated ends
        novel_slack (int): Distance b/w which to merge observed ends
        verbose (bool): Whether or not to print output

    Returns:
        t_df (pandas DataFrame): sg.t_df modified to include
            information about intron chain, tss, and tes
        regions (dict of pandas DataFrames): Indexed by
            'tss' and 'tes'. Bed regions for each end cluster
            as annotated in t_df
        counts (pandas DataFrame): DF of counts for intron
            chains, TSSs, TESs, and unique combinations of the three
            for annotated, observed, and both
    """

    all_df = add_tss_ic_tes(sg)

    # limit to those annotated or in list of detected tids that we allow
    # only ever allow known, nic, nnc
    _, inds = get_tpm_table(df,
                             how='iso',
                             sample='all',
                             min_tpm=1,
                             gene_subset='polya',
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
        reg = pr.pyranges.PyRanges(df=reg)
        ends = ends.join(reg, how='left',
                         slack=1,
                         strandedness=None, suffix='_annot')
        ends = ends.as_df()

        # limit to those w/ matching gid, first sd and add to cluster df
        if c == 'tss' and first_sd:
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
        if c == 'tss' and first_sd:
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
        nov_reg = pr.pyranges.PyRanges(df=nov_reg)
        temp = nov_reg.join(reg, how=None, strandedness=None, suffix='_annot')
        temp = temp.as_df()
        if verbose:
            if c == 'tss' and first_sd:
                temp = temp.loc[(temp.first_sd == temp.first_sd_annot)&(temp.gid == temp.gid_annot)]
            else:
                temp = temp.loc[temp.gid == temp.gid_annot]
            cols = ['gid']
            if c == 'tss' and first_sd:
                cols.append('first_sd')
            temp = temp.drop_duplicates(subset=cols)
            n = len(temp.index)
            print('{} new {} regions overlap annotated regions'.format(n,c))

        # add novel regions to clust and reg dfs
        reg = reg.as_df()
        nov_reg = nov_reg.as_df()
        clust = pd.concat([clust, nov_clust])
        reg = pd.concat([reg, nov_reg])

        # add strandedness to reg df
        temp2 = sg.t_df[['gid', 'path']]
        paths = temp2.path.values.tolist()
        first_edges = [path[0] for path in paths]
        temp2['first_edge'] = first_edges
        temp2 = temp2.merge(sg.edge_df[['strand']], how='left', left_on='first_edge', right_index=True)
        temp2 = temp2[['gid', 'first_edge', 'strand']].groupby(['gid', 'strand']).count().reset_index()
        temp2 = temp2[['gid', 'strand']]
        reg = reg.merge(temp2, how='left', on='gid')

        # some final formatting for these dfs
        cols = ['gid', 'gname', c,
                'Cluster', 'annotation']
        if c == 'tss' and first_sd:
            cols.append('first_sd')
        clust = clust[cols]
        clust.rename({'Cluster': '{}_cluster'.format(c),
                      'annotation': '{}_annotation'.format(c)},
                      axis=1, inplace=True)
        clust.drop_duplicates(inplace=True)
        end_regions[c] = reg

        # add cluster ids back into the original df
        cols = ['gid', 'gname', c]
        if c == 'tss' and first_sd:
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

    # get gene info and add
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname', 'biotype',
                  'biotype_category', 'tf']]
    gene_df.drop_duplicates(inplace=True)
    counts = counts.merge(gene_df, how='left', left_index=True, right_on='gid')

    t_df = all_df.copy(deep=True)

    # add tripletized transcript name to t_df
    t_df = get_gene_number(t_df, 'tss_cluster', 'tss')
    t_df = get_gene_number(t_df, 'tes_cluster', 'tes')
    t_df = get_gene_number(t_df, 'intron_chain', 'intron_chain')
    t_df['ttrip'] = t_df.gname +' ['+\
                t_df.tss_gene_num.astype('str')+','+\
                t_df.intron_chain_gene_num.astype('str')+','+\
                t_df.tes_gene_num.astype('str')+']'

    return t_df, end_regions, counts

def get_gene_number(df, col, pref):
    """
    Add number of tss / tes / ic occurrence to t_df
    from get_ic_tss_tes, based on number of unique
    values w/i the gene

    Parameters:
        df (pandas DataFrame): t_df from get_ic_tss_tes
        col (str): column with unique tss/ic/tes id
        pref (str): prefix to give new column

    Returns:
        df (pandas DataFrame): t_df with added numbers
            for tss/ic/tes id w/i gene
    """
    new_col = '{}_gene_num'.format(pref)
    temp = df[['gid', col, 'annotation']].copy(deep=True)
    temp.drop_duplicates(subset=['gid', col], inplace=True)
    temp[new_col] = temp.sort_values(['gid', col, 'annotation'],
                                 ascending=[True, True, False])\
                                 .groupby(['gid'])\
                                 .cumcount() + 1
    temp.drop('annotation', axis=1, inplace=True)
    df = df.merge(temp, how='left', on=['gid', col])
    return df

# def add_stable_gid(gtf):
#     """
#     Add stable gene id that accounts for PAR_X and PAR_Y
#     chromosomes to gtf df
#
#     Parameters:
#         gtf (pandas DataFrame): GTF dataframe
#
#     Returns:
#         gtf (pandas DataFrame): GTF dataframe with gene id turned into its
#             stable version
#     """
#     try:
#         gtf[['temp', 'par_region_1', 'par_region_2']] = gtf.gene_id.str.split('_', n=2, expand=True)
#         gtf['gene_id'] = gtf.gene_id.str.split('.', expand=True)[0]
#         gtf[['par_region_1', 'par_region_2']] = gtf[['par_region_1',
#                                                            'par_region_2']].fillna('')
#         gtf['gene_id'] = gtf.gene_id+gtf.par_region_1+gtf.par_region_2
#         gtf.drop(['temp', 'par_region_1', 'par_region_2'], axis=1, inplace=True)
#     except:
#         gtf['gene_id'] = gtf.gene_id.str.split('.', expand=True)[0]
#
#     return gtf

def add_feat(df, col, kind, as_index=False, as_number=False,
             drop_gid=True,
             drop_triplet=True):
    """
    Add ic, tss, tes info to a df from the transcript id
    
    Parameters:
        df (pandas DataFrame): Df w/ tid in there
        col (str): Which column to use
        kind (str): {'tss', 'ic', 'tes'}
        as_index (bool): Whether to replace current index with 
            new feat
        as_number (bool): Just add the number of element, not 
            geneid_#
    
    Returns:
        df (pandas DataFrame): DF w/ additional "kind" col
    """

    if col == 'index':
        df['temp_tid'] = df.index.tolist()
    else:
        df['temp_tid'] = df[col].tolist()
    df['triplet'] = df.temp_tid.str.split('[', expand=True)[1].str.split(']', expand=True)[0]
    df['temp_gid'] = df.temp_tid.str.split('[', expand=True)[0]
    if kind == 'tss':
        ind = 0
    elif kind == 'ic':
        ind = 1
    elif kind == 'tes':
        ind = 2
    df[kind] = df.triplet.str.split(',', expand=True)[ind]
    if not as_number:
        df[kind] = df.temp_gid+'_'+df[kind]
    drop_cols = list(['temp_tid'])
    if drop_gid:
        drop_cols += ['temp_gid']
    if drop_triplet:
        drop_cols += ['triplet']
    df.drop(drop_cols, axis=1, inplace=True)
    
    if as_index:
        df.reset_index(inplace=True, drop=True)
        df.index = df[kind]
        df.drop(kind, axis=1, inplace=True)
    return df
        
def get_gtf_info(how='gene',
                 subset=None,
                 add_stable_gid=False,
                 ver='v40_cerberus'):
    """
    Gets the info from the annotation about genes / transcripts

    Parameters:
        how (str): 'gene', 'iso', 'ic', 'tss', 'tes'
        subset (str): 'polya', 'tf', 'protein_coding' or None
        add_stable_gid (bool): Add stable gid (code from cerberus)
        ver (str): {'v29', 'v40_cerberus'}

    Returns:
        df (pandas DataFrame): DataFrame with info for gene / transcript
        biotype_counts (pandas DataFrame): DataFrame with the counts
            per biotype reported in gencode
        biotype_cat_counts (pandas DataFrame): DataFrame with the counts
            per meta biotype reported in gencode
    """
    iso_hows = ['iso', 'tss', 'tes', 'ic']
    d = os.path.dirname(__file__)
    if how == 'gene' and ver == 'v29':
        fname = '{}/../refs/gencode_v29_gene_metadata.tsv'.format(d)
    elif how in iso_hows and ver == 'v29':
        fname = '{}/../refs/gencode_v29_transcript_metadata.tsv'.format(d)
    elif how == 'gene' and ver == 'v40_cerberus':
        fname = '{}/../refs/cerberus/v40_gene_metadata.tsv'.format(d)
    elif how in iso_hows and ver == 'v40_cerberus':
        fname = '{}/../refs/cerberus/v40_transcript_metadata.tsv'.format(d)

    df = pd.read_csv(fname, sep='\t')
    
    # pdb.set_trace()

    if how == 'gene':
        id_col = 'gid'
    elif how == 'iso':
        id_col = 'tid'
    elif how == 'ic':
        iso_col = 'tid'
        id_col = 'ic'
    elif how == 'tss':
        iso_col = 'tid'
        id_col = 'tss'
    elif how == 'tes':
        iso_col = 'tid'
        id_col = 'tes'
    
    # if using cerberus features, drop duplicate  entries that come
    # from the different transcripts using the same features
    # also ignore the transcript length column
    if how in ['ic', 'tss', 'tes']:
        df = add_feat(df, kind=id_col, col=iso_col)
        df.drop([iso_col, 't_len'], axis=1, inplace=True)
        
        # double check
        n_feats = len(df[id_col].unique().tolist())
        n_drop_dupes = len(df.drop_duplicates().index)
        if n_feats != n_drop_dupes:
            print('Warning: number of unique {}s does not match length of deduped info table'.format(how))
        df.drop_duplicates(inplace=True)

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

    # add stable gid if requested
    if add_stable_gid:
        df['gid_stable'] = cerberus.get_stable_gid(df, 'gid')

    return df, biotype_counts, biotype_cat_counts

def get_det_table(df,
                  groupby='library',
                  min_tpm=0,
                  **kwargs):
    """
    Get a dataframe of True / False whether or not a gene / isoform
    was detected in a specific library or sample

    Parameters:
        df (pandas DataFrame): TALON abundance
        min_tpm (float): Min. TPM to be considered detected
        groupby (str): Either 'sample', 'library', or 'all'
            used to groupby datasets displayed

    Returns:
        df (pandas DataFrame): DataFrame with True / False entries
            for each isoform / gene per library / sample
    """
    if 'nov' not in kwargs:
        nov = df.transcript_novelty.unique().tolist()
    
    # add min_tpm to kwargs as it's needed for both 
    # functions
    kwargs['min_tpm'] = min_tpm

    # calc TPM per library on desired samples
    df, tids = get_tpm_table(df, **kwargs)

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

    if min_tpm != 0:
        df = (df >= min_tpm)
    # if min_tpm is 0, just take everything that has at least one read
    else:
        df = (df > min_tpm)
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

def get_ca_table(h5,
                 feat):
    ca = cerberus.read(h5)
    if feat == 'tss':
        df = ca.tss
    elif feat == 'tes':
        df = ca.tes
    elif feat == 'ic':
        df = ca.ic
    return df

def get_source_feats(h5,
                     feat, 
                     sources,
                     gene_subset):
    df = get_ca_table(h5, feat)
    df = filter_cerberus_sources(df, sources)

    if gene_subset: 
        gene_df, _, _ = get_gtf_info(how='gene', add_stable_gid=True, subset=gene_subset, ver='v40_cerberus')    
        df = df.loc[df.gene_id.isin(gene_df.gid_stable.tolist())]

    ids = df.Name.tolist()
    return ids
    
def get_det_feats(h5, 
                  filt_ab,
                  feat,
                  source,
                  **kwargs):
    
    """
    Get a list of ids corresponding to the queried 
    cerberus feature that are expressed
    
    Parameters:
        h5 (str): Path to cerberus annotation
        filt_ab (str): Path to filtered abundance file
        feat (str): {'tss', 'ic', 'tes'}
        source (str): Source to consider in cerberus obj
        
    Returns:
        ids (list of str): List of feature ids
    """
    
    # get these features from cerberus
    ca_df = get_ca_table(h5, feat)
    
    # get detected features
    df = pd.read_csv(filt_ab, sep='\t')
    df, ids = get_tpm_table(df, **kwargs)
    df = ca_df.loc[ca_df.Name.isin(ids)]
    
    return df.Name.tolist()

def get_tpm_table(df,
                    sample='all',
                    how='gene',
                    groupby='library',
                    nov=None,
                    min_tpm=0,
                    gene_subset=None,
                    save=False,
                    h5=None,
                    ic_nov=None,
                    tss_nov=None,
                    tes_nov=None,
                    **kwargs):
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
    df['gid_stable'] = cerberus.get_stable_gid(df, 'annot_gene_id')

    # merge with information about the gene
    gene_df, _, _ = get_gtf_info(how='gene', add_stable_gid=True, ver='v40_cerberus')
    gene_df = gene_df[['gid_stable', 'biotype_category', 'tf']]
    df = df.merge(gene_df, how='left', on='gid_stable')

    # get indices that we'll need to subset on
    if how == 'gene':
        id_col = 'gid_stable'
        nov_col = 'gene_novelty'
        nov = ['Known']
    elif how == 'iso':
        id_col = 'annot_transcript_id'
        nov_col = 'transcript_novelty'
    elif how == 'ic':
        id_col = 'ic'
        nov_col = 'transcript_novelty'
    elif how == 'tss':
        id_col = 'tss'
        nov_col = 'transcript_novelty'
    elif how == 'tes':
        id_col = 'tes'
        nov_col = 'transcript_novelty'
        
    # if we're looking at a cerberus feature, add that feature
    if how in ['tss', 'tes', 'ic']:
        df = add_feat(df, kind=how, col='annot_transcript_id')

    # filter on novelty
    if nov:
        print('Subsetting for novelty categories {}'.format(nov))
        nov_inds = df.loc[df[nov_col].isin(nov), id_col].tolist()
    else:
        nov_inds = df[id_col].tolist()
    
    # filter on ca feature novelties
    ca_inds = []
    for ca_feat, ca_novs in zip(['tss', 'ic', 'tes'], [tss_nov, ic_nov, tes_nov]):
        if h5 and ca_novs:
            print('Getting {} {}s'.format(ca_novs, ca_feat))
            if how != ca_feat:
                df = add_feat(df, col='annot_transcript_id', kind=ca_feat)
            ca_df = get_ca_table(h5, ca_feat)
            ca_df = ca_df.loc[ca_df.novelty.isin(ca_novs)]
            inds = df.loc[df[ca_feat].isin(ca_df.Name.tolist()), id_col].tolist()
            ca_inds.append(inds)
        else:
            ca_inds.append(df[id_col].tolist())
    feat_inds = list(set(ca_inds[0])&set(ca_inds[1])&set(ca_inds[2]))
            
    #     tss = get_ca_table(h5, 'tss')
    #     tss = tss.loc[tss.novelty.isin(tss_nov)]
    #     if how != 'tss':
    #         df = add_feat(df
    # elif h5 and tes_nov:
    #     tes = get_ca_table(h5, 'tes')
    #     tes = tes.loc[tes.novelty.isin(tes_nov)]
    # elif h5 and ic_nov:
    #     ic = get_ca_table(h5, 'ic')
    #     ic = ic.loc[ic.novelty.isin(ic_nov)]

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
    subset_inds = list(set(nov_inds)&set(gene_inds)&set(feat_inds))

    # sum up counts across the same gene, ic, tss, or tes
    sum_hows = ['gene', 'ic', 'tss', 'tes']
    if how in sum_hows:
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

def read_ends(h5,
              mode):
    """
    Read tss or tes bed file from cerberus

    Parameters:
        h5 (str): Path to file
        mode (str): {'tes', 'tss'}

    Returns:
        df
    """
    _, tss, tes, _, _, _ = read_h5(h5, as_pyranges=False)
    if mode == 'tss':
        df = tss
    elif mode == 'tes':
        df = tes

    return df

def read_h5(h5, as_pyranges=True):
    """
    Read h5 representation of a transcriptome

    Parameters:
        h5 (str): .h5 file to read from
        as_pyranges (bool): Convert bed representations to PyRanges objects

    Returns:
        ic (pandas DataFrame): Table detailing intron chains
        tss (pyranges PyRanges / pandas DataFrame): Bed representation of tss regions
        tes (pyranges PyRanges / pandas DataFrame): Bed represenation of tes regions
        tss_map (pyranges PyRanges / pandas DataFrame): Bed representation of
            all input tsss and which cerberus end they were assigned to
        tes_map (pyranges PyRanges / pandas DataFrame): Bed representation of
            all input tess and which cerberus end they were assigned to
        m (pandas DataFrame): Map of transcript id to tss / ic / tes
    """

    ic = pd.read_hdf(h5, key='ic')
    tss = pd.read_hdf(h5, key='tss')
    tes = pd.read_hdf(h5, key='tes')

    def read_empty_h5(h5, key, as_pyranges=False):
        try:
            df = pd.read_hdf(h5, key=key)
            if as_pyranges:
                df = pr.PyRanges(df)
        except:
            df = None
        return df

    m = read_empty_h5(h5, 'map', as_pyranges=False)
    tss_map = read_empty_h5(h5, 'tss_map', as_pyranges=as_pyranges)
    tes_map = read_empty_h5(h5, 'tes_map', as_pyranges=as_pyranges)

    if as_pyranges:
        tss = pr.PyRanges(tss)
        tes = pr.PyRanges(tes)

    return ic, tss, tes, tss_map, tes_map, m

def filter_cerberus_sources(df, sources):
    """
    Limit to entries that only come from given sources

    Parameters:
        df (pandas DataFrame): DF of ends or ICs from cerberus
        sources (list of str): Sources to retain
    """

    df['source'] = df.source.str.split(',')
    df = df.explode('source')
    df = df.loc[df.source.isin(sources)]
    df.drop('source', axis=1, inplace=True)
    df.drop_duplicates(inplace=True)
    return df

def filter_cerberus_genes(df, subset):
    """
    Limit to entries that only come from genes that map
    to specific biotypes

    Parameters:
        df (pandas DataFrame): DF of ends or ICs from cerberus
        subset (str): {'protein_coding', 'polya', 'tf', None}
    """

    # get gene info
    g_df, _, _ = get_gtf_info(how='gene', subset=subset, add_stable_gid=True)
    gids = g_df.gid_stable.tolist()

    # filter based on subset
    df = df.loc[df.gene_id.isin(gids)]

    return df

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

def count_gisx_region_genes(df, source, tss, tes, spl):
    df = df.loc[df.source == source].copy(deep=True)
    df['total'] = df.tss+df.tes+df.splicing_ratio
    df['tss_ratio'] = df.tss / df.total
    df['tes_ratio'] = df.tes / df.total
    df['spl_ratio'] = df.splicing_ratio / df.total

    t = len(df.index)
    print('{} genes are in {}'.format(t, source))

    # tss-high
    n = len(df.loc[df.tss_ratio > tss].index)
    print('{} ({:.2f}%) genes are TSS-high in {}'.format(n, (n/t)*100, source))

    # tes-high
    n = len(df.loc[df.tes_ratio > tes].index)
    print('{} ({:.2f}%) genes are TES-high in {}'.format(n, (n/t)*100, source))

    # splicing-high
    n = len(df.loc[df.spl_ratio > spl].index)
    print('{} ({:.2f}%) genes are splicing-high in {}'.format(n, (n/t)*100, source))

    # simple genes
    n = len(df.loc[(df.tss_ratio <= tss)&(df.tes_ratio <= tes)&(df.spl_ratio <= spl)].index)
    print('{} ({:.2f}%) genes are simple in {}'.format(n, (n/t)*100, source))

def count_gisx_region_genes(df, source, tss, tes, spl):
    df = df.loc[df.source == source].copy(deep=True)
    df['total'] = df.tss+df.tes+df.splicing_ratio
    df['tss_ratio'] = df.tss / df.total
    df['tes_ratio'] = df.tes / df.total
    df['spl_ratio'] = df.splicing_ratio / df.total

    t = len(df.index)
    print('{} genes are in {}'.format(t, source))

    # tss-high
    n = len(df.loc[df.tss_ratio > tss].index)
    print('{} ({:.2f}%) genes are TSS-high in {}'.format(n, (n/t)*100, source))

    # tes-high
    n = len(df.loc[df.tes_ratio > tes].index)
    print('{} ({:.2f}%) genes are TES-high in {}'.format(n, (n/t)*100, source))

    # splicing-high
    n = len(df.loc[df.spl_ratio > spl].index)
    print('{} ({:.2f}%) genes are splicing-high in {}'.format(n, (n/t)*100, source))

    # simple genes
    n = len(df.loc[(df.tss_ratio <= tss)&(df.tes_ratio <= tes)&(df.spl_ratio <= spl)].index)
    print('{} ({:.2f}%) genes are simple in {}'.format(n, (n/t)*100, source))

def assign_gisx_sector(df):
    if 'n_tss' in df.columns.tolist():
        tss = 'n_tss'
        tes = 'n_tes'
    else:
        tss = 'tss'
        tes = 'tes'
    spl = 'splicing_ratio'
    df['total'] = df[tss]+df[tes]+df[spl]
    df['tss_ratio'] = df[tss] / df.total
    df['tes_ratio'] = df[tes] / df.total
    df['spl_ratio'] = df[spl] / df.total

    df['sector'] = 'simple'

    df.loc[df.tss_ratio > 0.5, 'sector'] = 'tss'
    df.loc[df.tes_ratio > 0.5, 'sector'] = 'tes'
    df.loc[df.spl_ratio > 0.5, 'sector'] = 'splicing'
    
    # mixed genes
    df.loc[(df.sector=='simple')&(df.n_iso>1), 'sector'] = 'mixed'
    
    return df

def compare_species(h_counts, m_counts, source='obs'):

    d = os.path.dirname(__file__)
    fname = '{}/../refs/biomart_human_to_mouse.tsv'.format(d)
    conv = pd.read_csv(fname, sep='\t')

    # add sector
    h_counts = assign_gisx_sector(h_counts)
    m_counts = assign_gisx_sector(m_counts)

    h_counts = h_counts.loc[h_counts.source == source].copy(deep=True)
    m_counts = m_counts.loc[m_counts.source == source].copy(deep=True)
    conv = conv.drop_duplicates(subset=['Mouse gene stable ID', 'Mouse gene name'])

    # add non versioned gid
    m_counts['gid_2'] = m_counts.gid.str.rsplit('.', n=1, expand=True)[0]
    h_counts['gid_2'] = h_counts.gid.str.rsplit('.', n=1, expand=True)[0]

    # merge in with conversion table + mouse ID
    m_counts = m_counts.merge(conv, how='outer', left_on='gid_2', right_on='Mouse gene stable ID')

    # merge in with human counts
    h_counts = h_counts.merge(m_counts, how='outer', left_on='gid_2', right_on='Gene stable ID',
                              suffixes=('_human', '_mouse'))

    return h_counts
