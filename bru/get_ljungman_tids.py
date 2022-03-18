import pandas as pd
import numpy as np
import scipy.stats as st
import seaborn as sns
import sys
import os
import gseapy as gp
import swan_vis as swan
from scipy import sparse
import sqlite3

p = os.path.dirname((os.getcwd()))
sys.path.append(p)

from scripts.utils import *
from scripts.plotting import *

# get file names that belong to ljungman cell lines
# exp_metadata.tsv is from the report query given in the readme
exp = pd.read_csv('exp_metadata.tsv', sep='\t')
exp['files'] = exp.Files.str.split(',')
exp.files = exp.apply(lambda x: [f.split('/')[-2] for f in x.files], axis=1)
exp = exp[['files']]
exp = exp.explode(column='files')

# merge with file ids from lr bulk to limit to
# get dataset ids
df = pd.read_csv('../lr_bulk/file_to_hr.tsv', sep='\t',
                 header=None, names=['files', 'dataset'])
df = df.merge(exp, how='inner', on='files')
df = df[['dataset']]

# save datasets so we can use them when obtaining samples
df.to_csv('ljungman_datasets.tsv', sep='\t', header=None, index=False)

# output detected transcript / gene ids to pass list
df = pd.read_csv('../lr_bulk/talon/human_talon_abundance_filtered.tsv', sep='\t')
df = get_det_table(df,
                   how='iso',
                   min_tpm=0,
                   gene_subset='polya',
                   sample='ljungman',
                   groupby='library',
                   nov=['Known', 'NIC', 'NNC'])

tids = df.columns.tolist()
df = pd.read_csv('../lr_bulk/talon/human_talon_abundance_filtered.tsv', sep='\t')
tids = list(set(tids)|set(df.loc[df.transcript_novelty=='Known', 'annot_transcript_id'].tolist()))
df = df.loc[df.annot_transcript_id.isin(tids)]
df = df[['gene_ID', 'transcript_ID']]

# also get all known gene / transcript IDs straight from the db
db = '../lr_bulk/talon/human.db'
with sqlite3.connect(db) as conn:
    query = """SELECT DISTINCT t.gene_ID, t.transcript_ID
                    FROM transcripts as t
                    LEFT JOIN transcript_annotations as ta
                        ON ta.ID = t.transcript_ID
                    WHERE (ta.attribute = 'transcript_status'
                        AND ta.value = 'KNOWN')
            """
    annot = pd.read_sql_query(query, conn)

df = pd.concat([df, annot])
df.to_csv('ljungman_complete_pass_list.csv', header=None, index=False)
