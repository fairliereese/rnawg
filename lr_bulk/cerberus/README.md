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
# lapa, cage, rampage, ccre tss, ccre pels, ccre dels

# tes config - nothing yet
# lapa, pas seq

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

get intron chans
```bash
lapa_gtf=~/mortazavi_lab/data/rnawg/lr_bulk/lapa/human_lapa.gtf
out=~/mortazavi_lab/data/rnawg/lr_bulk/cerberus/test_ics.tsv
cerberus gtf_to_ics \
  --gtf ${lapa_gtf} \
  -o $out
```

# annotate transcriptome
h5=cerberus_ref.h5
# talon_gtf=~/mortazavi_lab/data/rnawg/lr_bulk/talon/human_known_nic_nnc_talon.gtf
lapa_gtf=~/mortazavi_lab/data/rnawg/lr_bulk/lapa/human_lapa.gtf
o=human_cerberus.h5
cerberus convert_transcriptome \
  --gtf ${lapa_gtf} \
  --h5 ${h5} \
  -o ${o}

# annotate transcriptome v29
h5=cerberus_ref.h5
gtf=~/mortazavi_lab/data/rnawg/refs/gencode_v29_sirv4_ercc.gtf
o=v29_cerberus.h5
cerberus convert_transcriptome \
  --gtf ${gtf} \
  --h5 ${h5} \
  -o ${o}


# annotate transcriptome v40
h5=cerberus_ref.h5
gtf=~/mortazavi_lab/ref/gencode.v40/gencode.v40.annotation.gtf
o=v40_cerberus.h5
cerberus convert_transcriptome \
  --gtf ${gtf} \
  --h5 ${h5} \
  -o ${o}

# replace ids in abundance
h5=cerberus_annot.h5
ab=../lapa/human_talon_abundance_filtered.corrected.tsv
o=human_cerberus_abundance.tsv
cerberus replace_ab_ids \
  --h5 ${h5} \
  --ab ${ab} \
  --collapse \
  -o ${o}

# replace ids in gtf
h5=cerberus_annot.h5
gtf=~/mortazavi_lab/data/rnawg/lr_bulk/lapa/human_lapa.gtf
o=human_cerberus.gtf
cerberus replace_gtf_ids \
  --h5 $h5 \
  --gtf $gtf \
  --update_ends \
  --collapse \
  -o $o

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
```
# Printed settings for agg_ends tes:
# ['/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/cerberus/temp/v40_tes.bed', '/Users/fairliereese/Documents/programming/mortazavi_
# lab/data/rnawg/lr_bulk/cerberus/temp/v29_tes.bed', '/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/cerberus/temp/lapa_tes.bed']
# ['v40', 'v29', 'lapa']
# [True, True, True]
# 20
# tes

## aggregate ends w/ lapa tss
```
# tss config - nothing yet
# lapa, cage, rampage, ccre tss, ccre pels, ccre dels

config=test_tss_config.csv
v40_bed=/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/cerberus/v40_tss.bed
v29_bed=/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/cerberus/v29_tss.bed
lapa_bed=/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/lapa/lapa_tss.bed

touch ${config}
printf "${v40_bed},True,v40\n" > $config
printf "${v29_bed},True,v29\n" >> $config
printf "${lapa_bed},True,lapa\n" >> $config

cerberus agg_ends \
  --input $config \
  --mode tss \
  -o test_agg_tss.bed
```

```
tss=agg_tss.bed
tes=agg_tes.bed
ics=agg_ics.tsv
ref=cerberus_ref.h5
cerberus write_reference \
        --tss ${tss} \
        --tes ${tes} \
        --ics ${ics} \
        -o ${ref}
```

```
# update gtf with cerberus
old_gtf=lr_bulk/lapa/human_lapa.gtf
new_gtf=lr_bulk/cerberus/cerberus.gtf
annot=lr_bulk/cerberus/cerberus_annot.h5
cerberus replace_gtf_ids \
    --h5 $annot \
    --gtf $old_gtf \
    --update_ends \
    --collapse \
    -o $new_gtf
```


Another test
```bash
config=test_tss_config_2.csv
v40_bed=/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/cerberus/v40_tss.bed
v29_bed=/Users/fairliereese/Documents/programming/mortazavi_lab/data/rnawg/lr_bulk/cerberus/v29_tss.bed

touch ${config}
printf "${v40_bed},True,v40\n" > $config
printf "${v29_bed},True,v29\n" >> $config

cerberus agg_ends \
  --input $config \
  --mode tss \
  -o test_agg_tss_2.bed
```
