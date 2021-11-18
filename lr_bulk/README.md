# Analysis of human ENCODE4 long read data

* [ENCODE cart of human cell line data](https://www.encodeproject.org/carts/723c6f14-e68b-4480-8a61-704a15ac5c7a/)
* [ENCODE cart of human tissue data](https://www.encodeproject.org/carts/26bd2879-329d-4168-98b9-6d132a1aad0f/)
* [ENCODE cart of all human data](https://www.encodeproject.org/carts/829d339c-913c-4773-8001-80130796a367/)
* [Query for human data](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=long+read+RNA-seq&lab.title=Ali+Mortazavi%2C+UCI&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.rfa=ENCODE4)

## Download processed files (post-TranscriptClean)
```bash
xargs -L 1 curl -O -J -L < files.txt
```

## Create human readable sample names and samples text file
Uses the manually created map of biosample term name + experiment ID to a more human-readable version for the cell line / in vitro differentiated cell data `exp_name_map.tsv`.
```bash
python format_metadata.py
```

## Move ENCODE accessioned files to human readable names
```bash
mkdir processing
while read sample
do
  file=`echo $sample | cut -f1 -d' '`.bam
  new=`echo $sample | cut -f2 -d' '`.bam
  mv $file $new
done < file_to_hr.tsv
mv *bam processing/
```
