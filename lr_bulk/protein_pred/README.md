# Using TAMA to predict ORFs and NMD from transcripts

TAMA was forked from [GitHub](https://github.com/GenomeRIK/tama) on 10/21/21 into [my version](https://github.com/fairliereese/tama).

gffread was cloned from [GitHub](https://github.com/gpertea/gffread) on 10/21/21.

blast cli was downloaded from [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/) on 10/21/21.

I used Hmmer version 3.3 which is available on the cluster.

I used bedtools version 2.29.2 which is also available on the cluster.

For this, I'll use the GTF that only includes known, NIC and NNC filtered transcripts.

```bash
gtf=../talon/human_talon_known_nic_nnc.gtf
sample=human
sbatch ../../scripts/orf_nmd_pred.sh $gtf $sample
```
