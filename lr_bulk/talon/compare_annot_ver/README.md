Make a gtf for confusing transcript:

```bash
printf "72,222627" > 222627_pass_list.csv
db=/Users/fairliereese/mortazavi_lab/data/rnawg/lr_bulk/talon/human.db
annot=gencode_v29
build=hg38
opref=ljungman
talon_create_GTF \
    --db ${db} \
    -a ${annot} \
    -b ${build} \
    --whitelist 222627_pass_list.csv \
    --o confusing_transcript
```

Get GTF subset from original GTF for confusing transcript:
```bash
grep "ENCODEHT000222627" ../human_known_nic_nnc_talon.gtf > confusing_transcript_from_gtf.gtf
```