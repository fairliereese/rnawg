import pandas as pd
import pyranges as pr

gr_gencode = pr.read_gtf(snakemake.input['gencode_gtf'])
gr_gencode = gr_gencode[gr_gencode.Feature == 'gene']

gr_mirbase = pr.read_gff3(snakemake.input['mirbase'])
gr_mirbase_transcript = gr_mirbase[gr_mirbase.Feature ==
                                   'miRNA_primary_transcript']
df_mirbase_miRNA = gr_mirbase[gr_mirbase.Feature == 'miRNA'].df

gr = gr_gencode.join(gr_mirbase_transcript, strandedness='same')
exact_match = (gr.Start == gr.Start_b) & (gr.End == gr.End_b)
gr = gr[exact_match]

df_families = pd.read_csv(snakemake.input['families'], sep='\t')
df_families = df_families[df_families['Species ID']
                          == snakemake.params['ncbi_tax_id']]
df_families = df_families[~df_families['MiRBase Accession'].isna()]

df_families = df_families.set_index('MiRBase Accession').join(
    df_mirbase_miRNA.set_index('ID')[['Derives_from']], how='inner').reset_index()

df_families = df_families.set_index('Derives_from').join(
    gr.df.set_index('ID')[['gene_id']], how='inner')

df_families[['miR family', 'gene_id']] \
    .drop_duplicates() \
    .to_csv(snakemake.output['families'], index=False)
