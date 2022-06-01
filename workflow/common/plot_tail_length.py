import matplotlib.pyplot as plt
from lapa.count import PolyaTailCounter


ttc = PolyaTailCounter(snakemake.input['bam'], min_tail_len=1)
plt.figure(figsize=(5, 5), dpi=200)

ttc.plot_tail_len_dist()
plt.xscale('log')
plt.savefig(snakemake.output['plot'])
