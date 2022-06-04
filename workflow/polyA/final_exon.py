import pandas as pd
import pyranges as pr


df_gtf = pr.read_gtf(snakemake.input['gtf'], as_df=True)

chroms = [f'chr{i}' for i in snakemake.config['chroms']]
df_gtf = df_gtf[df_gtf['Chromosome'].isin(chroms)]

df_abundance = pd.read_csv(snakemake.input['abundance'], sep='\t')

# Subset major transcript isoform
cols = df_abundance.columns[11:]
df_counts = df_abundance \
    .groupby(['annot_gene_id', 'annot_transcript_id']) \
    .sum()[cols] \
    .sum(axis=1)

transcripts = set(df_counts.groupby(level=0).idxmax().str.get(1))

# Subset relevant genes, transcript and exons
df_gtf = df_gtf[
    (df_gtf.gene_type == 'protein_coding') &
    (df_gtf.Feature == 'exon') &
    (df_gtf.transcript_id.str.contains('#')) &
    (df_gtf.transcript_id.isin(transcripts))
]

# Get final exons (3'utr) from each major transcript
df_final_exon = pd.concat([
    df_gtf[df_gtf.Strand == '+']
    .sort_values('End', ascending=False)
    .drop_duplicates('gene_id'),

    df_gtf[df_gtf.Strand == '-']
    .sort_values('Start')
    .drop_duplicates('gene_id')
])

# Get final exon from gtf
# final_exons = df_gtf \
#     .groupby('transcript_id')[['exon_number']] \
#     .max() \
#     .set_index('exon_number', append=True) \
#     .index
# df_final_exon = df_gtf \
#     .set_index(['transcript_id', 'exon_number']) \
#     .loc[final_exons] \
#     .reset_index()
gr = pr.PyRanges(df_final_exon[['Chromosome', 'Start',
                                'End', 'Strand', 'gene_id', 'transcript_id']])
gr.Name = gr.transcript_id
gr = gr.drop('transcript_id')
gr.to_bed(snakemake.output['bed'])
