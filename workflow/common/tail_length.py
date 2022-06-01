import matplotlib.pyplot as plt
from lapa.count import PolyaTailCounter


PolyaTailCounter(snakemake.input['bam'], min_tail_len=1) \
    .tail_len_dist() \
    .to_frame() \
    .reset_index() \
    .rename(columns={'index': 'tail_length', 0: 'count'}) \
    .to_csv(snakemake.output['counts'], index=False)
