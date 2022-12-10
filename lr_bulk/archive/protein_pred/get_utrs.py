import pandas as pd
import pyranges
import os
import sys
import seaborn as sns

p = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.append(p)

from scripts.utils import *

## gencode annotated transcripts
df = pr.read_gtf('../../refs/gencode_v29_sirv4_ercc.gtf')
df = df.df

# get rid of sirvs and erccs
df = df.loc[~df.Chromosome.str.contains('ERCC')]
df = df.loc[~df.Chromosome.str.contains('SIRV')]

# limit protein coding genes
gene_df = get_gtf_info(how='gene')[0]
gene_df = gene_df.loc[gene_df.biotype_category == 'protein_coding']
pc_genes = gene_df.gid.tolist()
print(len(pc_genes))

df = df.loc[df.gene_id.isin(pc_genes)]
len(df.loc[df.Feature == 'gene'])

t_df = df.loc[df.Feature == 'transcript'].copy(deep=True)
stop_df = df.loc[df.Feature == 'stop_codon'].copy(deep=True)
stop_df = stop_df[['transcript_id', 'Start', 'End']]

# add stop codon coords to transcript df
t_df = t_df.merge(stop_df, on='transcript_id',
                  suffixes=('','_stop_codon'))

# split into fwd and rev
fwd = t_df.loc[t_df.Strand == '+'].copy(deep=True)
rev = t_df.loc[t_df.Strand == '-'].copy(deep=True)

fwd['3_utr_start'] = fwd.End_stop_codon
fwd['3_utr_end'] = fwd.End

rev['3_utr_start'] = rev.Start
rev['3_utr_end'] = rev.Start_stop_codon

df = pd.concat([fwd, rev])
df = df[['Chromosome', '3_utr_start', '3_utr_end',
         'Strand', 'transcript_id']]
df.rename({'3_utr_start': 'Start',
           '3_utr_end': 'End',
           'transcript_id': 'Name'}, axis=1, inplace=True)

utr_df = df.copy(deep=True)

## novel transcripts
df = pd.read_csv('human_cds.bed', sep='\t',
                 header=None, usecols=[0,1,2,3,5,6,7])
df.columns = ['Chromosome', 'Start', 'Stop', 'fields', 'Strand',
              'CDS_Start', 'CDS_Stop']
df['tid'] = df.fields.str.split(';', expand=True)[1]

# get rid of sirvs and erccs
df = df.loc[~df.Chromosome.str.contains('ERCC')]
df = df.loc[~df.Chromosome.str.contains('SIRV')]

# check if transcript is novel
df['novel_transcript'] = df.tid.str.contains('ENCODE')

# check if transcript is nmd
df['nmd'] = ~df.fields.str.contains('prot_ok')

# limit to just novel transcripts that do not have NMD
df = df.loc[(df.novel_transcript==True)&(df.nmd==False)]

len(df.index)

# for transcripts without NMD, get the 3' UTR
fwd = df.loc[df.Strand == '+'].copy(deep=True)
rev = df.loc[df.Strand == '-'].copy(deep=True)

# fwd strand
fwd['3_utr_start'] = fwd.CDS_Stop+3
fwd['3_utr_end'] = fwd.Stop

# rev strand
rev['3_utr_start'] = rev.Start
rev['3_utr_end'] = rev.CDS_Start-3

df = pd.concat([fwd, rev])

df = df[['Chromosome', '3_utr_start', '3_utr_end', 'Strand', 'tid']]
df.rename({'3_utr_start': 'Start',
           '3_utr_end': 'End',
           'tid': 'Name'},
          axis=1, inplace=True)

utr_df = pd.concat([utr_df, df])
utr_df = pr.PyRanges(utr_df)

utr_df.to_bed('human_3_utr.bed')
