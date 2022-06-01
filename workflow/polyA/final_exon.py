import pyranges as pr


df_gtf = pr.read_gtf(snakemake.input['gtf'], as_df=True)

# Subset relevant genes, transcript and exons
df_gtf = df_gtf[
    (df_gtf.gene_type == 'protein_coding') &
    (df_gtf.Feature == 'exon') &
    (df_gtf.transcript_id.str.contains('#'))
]

# Get final exon from gtf
final_exons = df_gtf \
    .groupby('transcript_id')[['exon_number']] \
    .max() \
    .set_index('exon_number', append=True) \
    .index
df_final_exon = df_gtf \
    .set_index(['transcript_id', 'exon_number']) \
    .loc[final_exons] \
    .reset_index()

gr = pr.PyRanges(df_final_exon).merge(strand=True, by='gene_id')
gr.Name = gr.gene_id
gr = gr.drop('gene_id')
gr.to_bed(snakemake.output['bed'])
