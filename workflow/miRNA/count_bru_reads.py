import pyranges as pr


gr_bam = pr.read_bam(snakemake.input['bam'])

gr_utr = pr.read_bed(snakemake.input['utr'])

gr = gr_utr.count_overlaps(gr_bam)
gr.Score = gr.NumberOverlaps
gr = gr.drop('NumberOverlaps')

gr.to_bed(snakemake.output['counts'])
