import dask.dataframe as dd
from dask.diagnostics import ProgressBar


with ProgressBar():
    df = dd.read_csv(snakemake.input['dfs'])
    df.to_csv(snakemake.output['df'], index=False, single_file=True)
