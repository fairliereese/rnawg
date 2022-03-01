## Create read annot file with everything in it
```bash
oprefix=human
sbatch ../../scripts/sbatch_talon_read_annot.sh $oprefix # running 2/28/22
```

## Get filtered abundance / GTF for everything
```bash
db=human.db
opref=human
sbatch ../../scripts/sbatch_talon_ab_gtf.sh $db $opref # running 2/28/22
```
