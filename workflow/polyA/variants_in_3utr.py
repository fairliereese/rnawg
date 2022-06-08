import pandas as pd
import pyranges as pr
from tqdm import tqdm
from kipoiseq.extractors import MultiSampleVCF
from kipoiseq import Interval

vcf = MultiSampleVCF(snakemake.input['vcf'])

df = pr.read_bed(snakemake.input['bed']).merge(by='ThickStart').df \
    .rename(columns={'ThickStart': 'gene_id'})


chrom = snakemake.wildcards['chrom']
df = df[df.Chromosome == 'chr' + chrom]


df_counts = list()

for row in tqdm(df.itertuples(), total=df.shape[0]):

    interval = Interval(row.Chromosome, row.Start,
                        row.End, strand=row.Strand)
    variants = vcf.fetch_variants(interval)

    num_3utr_variants = 0

    for variant in variants:

        # Rare variants
        if (variant.source.INFO.get('AC', 0) > 0) \
           and (variant.source.INFO.get('AF', 0) < 0.001):
            num_3utr_variants += 1

    df_counts.append({
        'Chromosome': row.Chromosome,
        'Start': row.Start,
        'End': row.End,
        'Strand': row.Strand,
        'gene_id': row.gene_id,
        'count': num_3utr_variants
    })

pd.DataFrame(df_counts).to_csv(snakemake.output['csv'], index=False)
