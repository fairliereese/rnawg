## Create read annot file with everything in it
```bash
oprefix=human
sbatch ../../scripts/sbatch_talon_read_annot.sh $oprefix # running 4/14/22
```

## Get filtered abundance / GTF for everything
```bash
db=human.db
opref=human
sbatch ../../scripts/sbatch_talon_ab_gtf.sh $db $opref # running 4/14/22
```

## Make a version of the GTF w/o sirv/ercc
```bash
grep -v SIRV human_known_nic_nnc_talon.gtf > temp
grep -v chrEBV temp > human_known_nic_nnc_talon_ucsc.gtf # need to wait for above to run
rm temp
```
