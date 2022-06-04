import pandas as pd
import pyranges as pr


df_gtf = pr.read_gtf(snakemake.input['gtf'], as_df=True)
chroms = [f'chr{i}' for i in snakemake.config['chroms']]
df_gtf = df_gtf[df_gtf['Chromosome'].isin(chroms)]

df_cds = pr.read_bed(snakemake.input['cds'], as_df=True) \
    .rename(columns={'Name': '_transcript_id'}) \
    .set_index('_transcript_id')

# Subset relevant genes, transcript and exons
df_gtf = df_gtf[
    (df_gtf.transcript_type == 'protein_coding') &
    (df_gtf.Feature == 'exon') &
    df_gtf.transcript_id.str.contains('#')
]
df_gtf['_transcript_id'] = df_gtf['transcript_id'].str.split('#').str.get(0)
df_gtf = df_gtf.set_index('_transcript_id')


df_pos = df_gtf[df_gtf['Strand'] == '+'] \
    .sort_values('Start', ascending=False).drop_duplicates('transcript_id')

df_neg = df_gtf[df_gtf['Strand'] == '-'] \
    .sort_values('End').drop_duplicates('transcript_id')

df_pos = df_pos.join(df_cds, rsuffix='_cds').reset_index()
df_pos = df_pos[
    (df_pos['Start'] < df_pos['End_cds']) &
    (df_pos['End_cds'] < df_pos['End'])
]
df_pos['Start'] = df_pos['End_cds'] + 3

df_neg = df_neg.join(df_cds, rsuffix='_cds').reset_index()
df_neg = df_neg[
    (df_neg['Start'] < df_neg['Start_cds']) &
    (df_neg['Start_cds'] < df_neg['End'])
]
df_neg['End'] = df_neg['Start_cds'] - 3

df = pd.concat([df_pos, df_neg])

core_cols = ['Chromosome', 'Start', 'End', 'Strand']

df = df[[*core_cols, 'transcript_id']]
df = df.groupby(core_cols) \
       .agg(lambda x: ';'.join(x)) \
       .reset_index() \
       .rename(columns={'transcript_id': 'Name'})

pr.PyRanges(df).to_bed(snakemake.output['bed'])
