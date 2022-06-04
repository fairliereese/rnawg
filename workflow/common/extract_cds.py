import pyranges as pr


df_gtf = pr.read_gtf(snakemake.input['gtf'], as_df=True)
df_gtf = df_gtf[df_gtf['Feature'] == 'CDS']

df = df_gtf.groupby('transcript_id').agg({
    'Chromosome': 'first',
    'Start': 'min',
    'End': 'max',
    'Strand': 'first',
    'transcript_id': 'first'
}).rename(columns={'transcript_id': 'Name'})

pr.PyRanges(df).to_bed(snakemake.output['bed'])
