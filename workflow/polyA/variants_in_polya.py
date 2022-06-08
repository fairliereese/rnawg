from tqdm import tqdm
import pandas as pd
import pyranges as pr
from lapa.result import LapaResult
from kipoiseq import Interval
from kipoiseq.extractors import MultiSampleVCF


df_3utr = pr.read_bed(snakemake.input['bed']).merge(by='ThickStart').df \
    .rename(columns={'ThickStart': 'gene_id'})
pos_strand = df_3utr['Strand'] == '+'
df_3utr.loc[pos_strand, 'Start'] = df_3utr.loc[pos_strand, 'End'] - 1
df_3utr.loc[~pos_strand, 'End'] = df_3utr.loc[~pos_strand, 'Start'] + 1
gr_3utr = pr.PyRanges(df_3utr)

df_cluster = LapaResult(snakemake.input['lapa_dir']).read_clusters()
df_cluster = df_cluster.drop_duplicates(
    ['Chromosome', 'Start', 'End', 'Strand'])
chrom = snakemake.wildcards['chrom']
df_cluster = df_cluster[df_cluster.Chromosome == 'chr' + chrom]
df_cluster = df_cluster[~(df_cluster.signal == 'None@None')]

df = pr.PyRanges(df_cluster).join(gr_3utr).df
df = df[df['gene_id'] == df['gene_id_b']]

vcf = MultiSampleVCF(snakemake.input['vcf'])

df['signal_seq'] = df['signal'].map(lambda x: x.split('@')[1])
df['signal_Start'] = df['signal'].map(lambda x: x.split('@')[0]).astype(int)
df['signal_End'] = df['signal_Start'] + 6

df_counts = list()

for row in tqdm(df.itertuples(), total=df.shape[0]):

    interval = Interval(row.Chromosome, row.signal_Start,
                        row.signal_End, strand=row.Strand)
    variants = vcf.fetch_variants(interval)

    df_counts.append({
        'Chromosome': row.Chromosome,
        'Start': row.signal_Start,
        'End': row.signal_End,
        'Strand': row.Strand,
        'gene_id': row.gene_id,
        'signal': row.signal_seq,
        'count': sum(
            (v.source.INFO.get('AC', 0) > 0)
            and (v.source.INFO.get('AF', 0) < 0.001)
            for v in variants
        )
    })

pd.DataFrame(df_counts).to_csv(snakemake.output['csv'], index=False)
