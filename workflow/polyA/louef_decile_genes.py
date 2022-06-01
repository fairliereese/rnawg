import pandas as pd
import pyranges as pr


df_gtf = pr.read_gtf(snakemake.input['gtf'], as_df=True)
df_gtf = df_gtf[
    (df_gtf.Feature == 'gene') &
    (df_gtf.gene_type == 'protein_coding')
]
df_gtf['_gene_id'] = df_gtf['gene_id'].str.split('.').str.get(0)
df_gtf = df_gtf[['_gene_id', 'gene_id']]

df = pd.read_csv(snakemake.input['lof'], sep='\t')

df = df[~df['oe_lof_upper_bin'].isna()][['oe_lof_upper_bin', 'gene_id']] \
    .drop_duplicates('gene_id') \
    .rename(columns={'gene_id': '_gene_id'}) \
    .set_index('_gene_id') \
    .join(df_gtf.set_index('_gene_id'))

df = df[~df['gene_id'].isna()]
df['oe_lof_upper_bin'] = df['oe_lof_upper_bin'].astype(int)

df[['gene_id', 'oe_lof_upper_bin']].to_csv(
    snakemake.output['csv'], index=False)
