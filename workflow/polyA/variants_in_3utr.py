import pandas as pd
import pyranges as pr
from tqdm import tqdm
from kipoiseq.extractors import MultiSampleVCF
from kipoiseq import Interval

vcf = MultiSampleVCF(snakemake.input['vcf'])
df = pr.read_bed(snakemake.input['bed'], as_df=True)

chrom = snakemake.wildcards['chrom']
df = df[df.Chromosome == 'chr' + chrom]


df_counts = list()

for row in tqdm(df.itertuples(), total=df.shape[0]):

    interval = Interval(row.Chromosome, row.Start,
                        row.End, strand=row.Strand)
    variants = vcf.fetch_variants(interval)

    df_counts.append({
        'Chromosome': row.Chromosome,
        'Start': row.Start,
        'End': row.End,
        'Strand': row.Strand,
        'gene_id': row.Name,
        'count': sum(
            (v.source.INFO.get('AC', 0) > 0)
            and (v.source.INFO.get('AF', 0) < 0.001)
            for v in variants
        )
    })

pd.DataFrame(df_counts).to_csv(snakemake.output['csv'], index=False)
