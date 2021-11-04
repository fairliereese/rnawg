# FASTA

```bash
# sirv4
wget https://www.encodeproject.org/files/ENCFF473RJX/@@download/ENCFF473RJX.fasta.gz -O sirv4.fa.gz
gunzip sirv4.fa.gz

# ercc
wget https://www.encodeproject.org/files/ENCFF001RTP/@@download/ENCFF001RTP.fasta.gz -O ercc.fa.gz
gunzip ercc.fa.gz

# human reference
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O hg38.fa.gz
gunzip hg38.fa.gz

# cat them all together
touch hg38_sirv4_ercc.fa
cat hg38.fa >> hg38_sirv4_ercc.fa
cat sirv4.fa >> hg38_sirv4_ercc.fa
cat ercc.fa >> hg38_sirv4_ercc.fa
```

# GTF
Links to the ref that includes sirv and ercc
https://www.synapse.org/#!Synapse:syn25683628
https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/


```bash
# pull sirv and ercc from here
gunzip lrgasp.gtf.gz
awk 'index($1, "ERCC")' lrgasp.gtf > ercc.gtf
awk 'index($1, "SIRV")' lrgasp.gtf > sirv.gtf

# gencode
wget https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/@@download/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz -O gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

# cat them all together
touch gencode_v29_sirv4_ercc.gtf
cat gencode.v29.annotation.gtf >> gencode_v29_sirv4_ercc.gtf
cat ercc.gtf >> gencode_v29_sirv4_ercc.gtf
cat sirv.gtf >> gencode_v29_sirv4_ercc.gtf
```

How many transcripts are in GENCODE?
```bash
cat gencode.v29.annotation.gtf | cut -f3 | grep transcript | wc -l
```
206,761 transcripts

How many genes are in GENCODE?
```bash
cat gencode.v29.annotation.gtf | cut -f3 | grep gene | wc -l
```
58,780 genes


Add another column to tissue_metadata.csv
```python
import pandas as pd
df = pd.read_csv('tissue_metadata.csv')
df['biosample'] = df.biosample_name.str.replace(' ', '_')
df.biosample = df.biosample.str.lower()

df.to_csv('tissue_metadata.csv', index=False)
```
