from kipoiseq.extractors import MultiSampleVCF
from lapa.results import LapaResult


vcf = MultiSampleVCF(snakemake.input['vcf'])


lapa = LapaResult(snakemake.input['lapa_dir'])

# TODO:
# - get polya seq site
# - calculate number of variants per gene
# - save numbers per gene
