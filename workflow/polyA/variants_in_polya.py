from lapa.result import LapaResult


lapa = LapaResult(snakemake.input['lapa_dir'])

df_cluster = lapa.read_clusters()

__import__("pdb").set_trace()
