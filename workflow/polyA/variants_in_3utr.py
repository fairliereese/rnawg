import pandas as pd
import pyranges as pr
from tqdm import tqdm
from kipoiseq.extractors import MultiSampleVCF
from kipoiseq import Interval

vcf = MultiSampleVCF(snakemake.input['vcf'])
df = pr.read_bed(snakemake.input['bed'], as_df=True) \
       .rename(columns={'Name': 'transcript_id', 'ThickStart': 'gene_id'})


chrom = snakemake.wildcards['chrom']
df = df[df.Chromosome == 'chr' + chrom]


df_counts = list()

for row in tqdm(df.itertuples(), total=df.shape[0]):

    interval = Interval(row.Chromosome, row.Start,
                        row.End, strand=row.Strand)
    variants = vcf.fetch_variants(interval)
    transcript_id = row.transcript_id.split('.')[0]

    num_3utr_variants = 0

    for variant in variants:

        # Rare variants
        if (variant.source.INFO.get('AC', 0) > 0) \
           and (variant.source.INFO.get('AF', 0) < 0.001):

            # parse consequences
            annotations = variant.source.INFO.get('vep') \
                                             .split(f',{variant.alt}|')
            annotations[0] = annotations[0].replace('{variant.alt}|', '')

            for annotation in annotations:
                annotation = annotation.split('|')

                vep_transcript_id = annotation[5]

                if transcript_id == vep_transcript_id:
                    csq = annotation[0]
                    if csq == '3_prime_UTR_variant':
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
