## Download LAPA-tized data from Hasan

```bash
lapa_gtf=/dfs3b/samlab/mcelik/rnawg/data/results/talon/human_talon.corrected.gtf
lapa_filt_ab=/dfs3b/samlab/mcelik/rnawg/data/results/talon/human_talon_abundance_filtered.corrected.tsv
scp freese@hpc3.rcic.uci.edu:${lapa_gtf} .
scp freese@hpc3.rcic.uci.edu:${lapa_filt_ab} .
```

## Generate version of GTF with just detected (>=1 TPM) known, nic, nnc
```bash
annot=gencode_v29
build=hg38
db=../talon/human.db

talon_create_GTF --db ${db} \
    -a ${annot} \
    -b ${build} \
    --whitelist human_talon_swan_pass_list.csv \
    --o human_swan
```
