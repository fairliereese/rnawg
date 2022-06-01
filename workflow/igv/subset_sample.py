import pandas as pd
import pyranges as pr


samples = snakemake.params['talon_samples']

df_abundance = pd.read_csv(snakemake.input['abundance'], sep='\t')
gr_gtf = pr.read_gtf(snakemake.input['gtf'])

df_abundance = df_abundance[
    df_abundance[samples].sum(axis=1) >= snakemake.params['cutoff']
]

expressed_transcripts = set(df_abundance['annot_transcript_id'])
expressed_genes = gr_gtf[
    gr_gtf.transcript_id.isin(expressed_transcripts) &
    (gr_gtf.Feature == 'transcript')
].gene_id

gr_gtf = gr_gtf[gr_gtf.gene_id.isin(expressed_genes)]
gr_gtf = gr_gtf[
    gr_gtf.transcript_id.isna() |
    gr_gtf.transcript_id.isin(expressed_transcripts)
]
gr_gtf.to_gtf(snakemake.output['gtf'])
