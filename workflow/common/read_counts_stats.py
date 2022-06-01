from pathlib import Path
import pandas as pd


counts = dict()

for i in snakemake.input['counts']:
    counts[Path(i).stem] = int(open(i).read())

df_mapping = pd.read_csv(snakemake.input['mapping'], sep='\t', header=None)

df_read_annot_counts = pd.read_csv(snakemake.input['read_annot_counts'], sep='\t')

df = df_mapping.set_index(0).join(pd.Series(counts).to_frame())

import pdb
pdb.set_trace()


df.to_csv(snakemake.output['table'], sep='\t')