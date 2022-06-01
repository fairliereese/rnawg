import pandas as pd


df = snakemake.params['df_short']

df = df['Biosample term name'].reset_index().rename(
    columns={'Biosample term name': 'biosample',
             'File accession': 'encode_id'})

df.to_csv(snakemake.output['csv'], index=False)
