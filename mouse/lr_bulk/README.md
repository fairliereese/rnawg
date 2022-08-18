<!-- * [Experiment query](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&perturbed=false&assay_title=long+read+RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&lab.title=Ali+Mortazavi%2C+UCI&biosample_ontology.term_name=adrenal+gland&biosample_ontology.term_name=layer+of+hippocampus&biosample_ontology.term_name=left+cerebral+cortex&biosample_ontology.term_name=gastrocnemius&biosample_ontology.term_name=heart&status=released&replicates.library.biosample.donor.accession=ENCDO509HIY)
* [Report for experiment query](https://www.encodeproject.org/report.tsv?type=Experiment&control_type!=*&perturbed=false&assay_title=long+read+RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&lab.title=Ali+Mortazavi%2C+UCI&biosample_ontology.term_name=adrenal+gland&biosample_ontology.term_name=layer+of+hippocampus&biosample_ontology.term_name=left+cerebral+cortex&biosample_ontology.term_name=gastrocnemius&biosample_ontology.term_name=heart&status=released&replicates.library.biosample.donor.accession=ENCDO509HIY)
* [Cart with long-read mouse timecourse datasets](https://www.encodeproject.org/carts/a3ffb484-07f3-468c-8cac-ecfefcd93c86/) -->

* [Experiment query](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&perturbed=false&assay_title=long+read+RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&lab.title=Ali+Mortazavi%2C+UCI&status=released)
* [Report for experiment query](https://www.encodeproject.org/report.tsv?type=Experiment&control_type!=*&perturbed=false&assay_title=long+read+RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&lab.title=Ali+Mortazavi%2C+UCI&status=released)
* [Cart with all long-read mouse data](https://www.encodeproject.org/carts/55367842-f225-45cf-bfbe-5ba5e4182768/)

## Download processed files (post-TranscriptClean) from the portal

```bash
xargs -L 1 curl -O -J -L < lr_post_tc_files.txt
```

## Create human readable sample names and samples text file
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
opref=~/mortazavi_lab/data/rnawg/mouse/lr_bulk/processing/
samples=~/mortazavi_lab/data/rnawg/mouse/lr_bulk/samples.txt
bash ../scripts/talon_label.sh $opref $samples
```

## Create TALON config file
```bash
python make_config.py
```

## Run TALON
```bash
config=talon/talon_config.csv
oprefix=talon/mouse
sbatch ../scripts/sbatch_talon_bulk.sh $config $oprefix
```

```bash
snakemake -s workflow/mouse/lr_bulk/Snakefile -j 10 --latency-wait 120 --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END, --time=72:00:00" -n

snakemake -s workflow/mouse/lr_bulk/Snakefile -j 100 --latency-wait 120 -n
```


<!-- Using more manageable chunks of data for TALON - 14 datasets at a time
```bash
oprefix=talon/mouse
sbatch ../scripts/sbatch_talon_bulk.sh talon/talon_config_1.csv $oprefix # done 10/25/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_2.csv $oprefix # done 10/26/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_3.csv $oprefix # done 10/27/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_4.csv $oprefix # done 10/27/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_5.csv $oprefix # done 10/27/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_6.csv $oprefix # done 3/29/22
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_7.csv $oprefix # done? 3/31/22 ran out of memory but db finished updating at 9:06 and OOM exception thrown at 9:25
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_8.csv $oprefix # done 4/1/22 NOT april fools joke
``` -->

<!--
Update with most recently released data 12/16/21
```bash
# first store old bois
mv talon/talon_config_1.csv talon/211216_talon_config_1.csv
mv talon/talon_config_2.csv talon/211216_talon_config_2.csv
mv talon/talon_config_3.csv talon/211216_talon_config_3.csv
mv talon/talon_config_4.csv talon/211216_talon_config_4.csv
mv talon/talon_config_5.csv talon/211216_talon_config_5.csv

oprefix=talon/mouse
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_1.csv $oprefix # done 12/17/21
sbatch ../scripts/sbatch_talon_bulk.sh -d talon/mouse.db talon/talon_config_2.csv $oprefix # running 12/17/21

```
-->

<!-- ## Download end-corrected GTF
Hasan corrected the ends of the transcripts and provided a GTF (2022/2/2)
[Link](https://drive.google.com/file/d/1ymdbWtNkpO_jnqSOg7HoDj0yDRUPJpu1/view?usp=sharing)

Download it using the browser and the move it to the reference folder.
```bash
mv ~/Downloads/mouse_talon.corrected.gtf ../refs/mouse_lr_bulk_corrected.gtf
``` -->
