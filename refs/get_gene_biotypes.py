import pandas as pd

fname = '/dfs6/pub/freese/mortazavi_lab/data/rnawg/refs/gencode_v29_sirv4_ercc.gtf'

df = pd.read_csv(fname, sep='\t', usecols=[0,2,3,4,8], comment='#')
df.columns = ['chr', 'entry_type', 'start', 'stop', 'fields']
print(df.head())

# remove sirvs and erccs
print(len(df.index))
df = df.loc[(~df.chr.str.contains('SIRV'))&~(df.chr.str.contains('ERCC'))]
print(len(df.index))

# only genes
df = df.loc[df.entry_type == 'gene'].copy(deep=True)

df['gid'] = df.fields.str.split('gene_id "', n=1, expand=True)[1]
df['gid'] = df.gid.str.split('";', n=1, expand=True)[0]

df['gname'] = df.fields.str.split('gene_name "', n=1, expand=True)[1]
df['gname'] = df.gname.str.split('";', n=1, expand=True)[0]

df['biotype'] = df.fields.str.split('gene_type "', n=1, expand=True)[1]
df['biotype'] = df.biotype.str.split('";', n=1, expand=True)[0]

# biotype map
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
                      'transcribed_unprocessed_pseudogene',
                      'processed_pseudogene',
                      'transcribed_processed_pseudogene',
                      'transcribed_unitary_pseudogene',
                      'unitary_pseudogene',
                      'polymorphic_pseudogene',
                      'pseudogene',
                      'translated_processed_pseudogene',
                      'translated_unprocessed_pseudogene'],
       'miRNA': ['miRNA'],
       'other': ['vault_RNA', 'snRNA',
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
df['length'] = (df.start-df.stop).abs()

# add TF info
df['tf'] = False
tf_df = pd.read_csv('biomart_tf_gids.tsv', sep='\t')
tf_gids = tf_df['Gene stable ID'].unique().tolist()
df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
df.drop('gid_stable', axis=1, inplace=True)

# and save
df = df[['gid', 'gname', 'length', 'biotype', 'biotype_category', 'tf']]
fname = 'gencode_v29_gene_metadata.tsv'
df.to_csv(fname, sep='\t', index=False)

# # only transcripts
# df['tid'] = df.fields.str.split('transcript_id "', n=1, expand=True)[0]
# df['tid'] = df.tid.str.split('";', n=1, expand=True)[0]
#
# df['biotype'] = df.fields.str.split('gene_type "', n=1, expand=True)[0]
# df['biotype'] = df.fields.str.split('";', n=1, expand=True)[0]
