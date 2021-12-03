## Create read annot file with everything in it
```bash
oprefix=human
sbatch ../../scripts/sbatch_talon_read_annot.sh $oprefix # running 11/29/21
```

## Get filtered abundance / GTF for everything
```bash
db=human.db
opref=human
sbatch ../../scripts/sbatch_talon_ab_gtf.sh $db $opref
```
