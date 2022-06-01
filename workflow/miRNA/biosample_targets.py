import numpy as np
import pandas as pd
import pyranges as pr


df_families = pd.read_csv(snakemake.input['families']).set_index('gene_id')

df_counts = pd.read_csv(snakemake.input['counts']).set_index('gene_id')
del df_counts['gene_name']

df_targetscan = pr.read_bed(snakemake.input['targetscan'], as_df=True)
df_targetscan['miR family'] = df_targetscan['Name'].str.split(':').str.get(1)

df = df_counts.join(df_families, how='inner') \
              .reset_index().groupby('miR family').sum()

df = df_targetscan.set_index('miR family') \
    .join(df)

gr_exp = pd.read_csv(snakemake.input['exp_bed'], sep='\t', header=None)
gr_exp = gr_exp.rename(columns={
    0: 'Chromosome',
    1: 'Start',
    2: 'End',
    5: 'Strand',
    9: 'cell_line',
    10: 'num_dataset'
})
gr_exp.columns = gr_exp.columns.astype(str)
gr_exp = pr.PyRanges(gr_exp, int64=True)
gr_exp = gr_exp.drop(['3', '4', '6', '7', '8', '11', '12'])

df = pr.PyRanges(df.reset_index(), int64=True).join(
    gr_exp, how='outer', strandedness='same').df

df = df[['Chromosome', 'Start', 'End', 'Strand', 'miR family',
         'cell_line', 'num_dataset', *df_counts]]
df = df[df['Start'] != -1].replace('-1', '.')

df.to_csv(snakemake.output['all_targets'], index=False)


for biosample in df_counts.columns:
    _df = df[[
        'Chromosome', 'Start', 'End', 'miR family', biosample,
        'Strand', 'cell_line', 'num_dataset'
    ]].rename(columns={biosample: 'Score', 'miR family': 'Name'})

    _df['Score'] = (_df['Score'] * 1_000_000) / _df['Score'].sum().round(2)
    _df = _df[_df['Score'] > snakemake.params['tpm_cutoff']]

    pr.PyRanges(_df).sort().to_bed(snakemake.config['miRNA']['targets']
                                   .format(biosample=biosample))
