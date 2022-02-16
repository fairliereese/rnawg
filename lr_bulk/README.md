# Analysis of human ENCODE4 long read data

* [ENCODE cart of human cell line data](https://www.encodeproject.org/carts/723c6f14-e68b-4480-8a61-704a15ac5c7a/)
* [ENCODE cart of human tissue data](https://www.encodeproject.org/carts/26bd2879-329d-4168-98b9-6d132a1aad0f/)
* [ENCODE cart of all human data](https://www.encodeproject.org/carts/829d339c-913c-4773-8001-80130796a367/)
* [Query for human data](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&perturbed=false&assay_title=long+read+RNA-seq&lab.title=Ali+Mortazavi%2C+UCI&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.rfa=ENCODE4&replicates.library.nucleic_acid_term_name=polyadenylated+mRNA&limit=all)

## Download processed files (post-TranscriptClean)
```bash
xargs -L 1 curl -O -J -L -n < files.txt
```

## Create human readable sample names and samples text file

Create manual map of biosample term name + experiment ID to a more human-readable version
```bash
python make_temp_biosamp_term_name_map.py
```

Manually edit missing shorthand entries with desired names.

Then rename file.
```bash
mv temp_biosamp_term_name_map.tsv biosamp_term_name_map.tsv
```

Uses the manually created map of biosample term name + experiment ID (`biosamp_term_name_map.tsv`) to a more human-readable version for the cell line / in vitro differentiated cell data `exp_name_map.tsv`.
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

## TALON label reads
```bash
opref=~/mortazavi_lab/data/rnawg/lr_bulk/processing/
samples=~/mortazavi_lab/data/rnawg/lr_bulk/samples  .txt
bash ../scripts/talon_label.sh $opref $samples
```

<!-- Check to see which things finished
```bashqs
tail processing/talon_label.o* | grep -B 8 "Run complete" | grep "talon_label.o"
``` -->

## Create TALON config file
```bash
mkdir talon
python make_config.py
```

## Run TALON
```bash
config=talon/talon_config.csv
oprefix=talon/human
sbatch ../scripts/sbatch_talon_bulk.sh $config $oprefix
```

<!-- Using more manageable chunks of data for TALON - 14 datasets at a time
```bash
oprefix=talon/human
sbatch ../scripts/sbatch_talon_bulk.sh talon/talon_config_1.csv $oprefix # finished 11/24/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/human.db talon/talon_config_2.csv $oprefix # done 11/25/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/human.db talon/talon_config_3.csv $oprefix # done 11/25/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/human.db talon/talon_config_4.csv $oprefix # done 11/27/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/human.db talon/talon_config_5.csv $oprefix # done 11/27/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/human.db talon/talon_config_6.csv $oprefix # done 11/28/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/human.db talon/talon_config_7.csv $oprefix # done 11/29/21
``` -->
