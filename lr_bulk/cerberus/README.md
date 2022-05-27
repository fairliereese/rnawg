Get CAGE ends
```bash
scp freese@hpc3.rcic.uci.edu:/data/resources/cage/merged.bed cage_tss.bed
```


Get ends from GTFs:
* GENCODE v39
* GENCODE v29
* TALON GTF


```bash
# gtf config
v40_gtf=~/mortazavi_lab/ref/gencode.v40/gencode.v40.annotation.gtf
v29_gtf=~/mortazavi_lab/data/rnawg/refs/gencode_v29_sirv4_ercc.gtf
# talon_gtf=~/mortazavi_lab/data/rnawg/lr_bulk/talon/human_known_nic_nnc_talon.gtf
lapa_gtf=~/mortazavi_lab/data/rnawg/lr_bulk/lapa/human_swan_talon.corrected.gtf
gtf_config=gtf_config.csv
printf "${v40_gtf},True,v40\n" > ${gtf_config}
printf "${v29_gtf},True,v29\n" >> ${gtf_config}
# printf "${talon_gtf},True,talon\n" >> ${gtf_config}
printf "${lapa_gtf},True,lapa\n" >> ${gtf_config}


# tss config - nothing yet

# tes config - nothing yet

h5=cerberus_ref.h5

cerberus gen_reference \
  --ref_gtf ${gtf_config} \
  -o ${h5} \
  --gtf_tss_dist 50 \
  --gtf_tss_slack 50 \
  --gtf_tes_dist 50 \
  --gtf_tes_slack 50 \
  --tss_slack 20 \
  --tes_slack 20 \
  --keep_tmp

# annotate transcriptome
h5=cerberus_ref.h5
talon_gtf=~/mortazavi_lab/data/rnawg/lr_bulk/talon/human_known_nic_nnc_talon.gtf
o=human_cerberus.h5
cerberus convert_transcriptome \
  --gtf ${talon_gtf} \
  --h5 ${h5} \
  -o ${o}

# replace ids in abundance
h5=human_cerberus.h5
ab=~/mortazavi_lab/data/rnawg/lr_bulk/talon/human_talon_abundance_filtered.tsv
o=human_cerberus.tsv
cerberus replace_ab_ids \
  --h5 ${h5} \
  --ab ${ab} \
  --collapse \
  -o ${o}
```

Other cmds
```bash
#Get tsss
cerberus gtf_to_bed \
  --gtf ${v39_gtf} \
  --mode tss \
  -o v39_tss.bed

cerberus gtf_to_bed \
  --gtf ${v29_gtf} \
  --mode tss \
  -o v29_tss.bed

cerberus gtf_to_bed \
  --gtf ${talon_gtf} \
  --mode tss \
  -o talon_tss.bed

#Aggregate TSSs
cfg=agg_tss_config_2.csv
cerberus agg_ends \
  --input ${cfg} \
  --mode tss \
  -o test_tss.bed

cfg=agg_tes_config_2.csv
cerberus agg_ends \
  --input ${cfg} \
  --mode tes \
  -o test_tes.bed
