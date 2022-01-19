# Analysis of human ENCODE3 and 4 short read RNA-seq data
* [ENCODE cart of all human data](https://www.encodeproject.org/carts/789031a1-e247-4e0a-b699-7a885e22ea47/)
* [Query for human data](https://www.encodeproject.org/matrix/?type=Experiment&control_type!=*&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=total+RNA-seq&assay_title=polyA+plus+RNA-seq&assembly=GRCh38&award.project=ENCODE&lab.title=Thomas+Gingeras%2C+CSHL&lab.title=Barbara+Wold%2C+Caltech&files.file_type=tsv&replicates.library.biosample.subcellular_fraction_term_name!=cytosol&replicates.library.biosample.subcellular_fraction_term_name!=nucleus&replicates.library.biosample.subcellular_fraction_term_name!=nucleolus&replicates.library.biosample.subcellular_fraction_term_name!=chromatin&replicates.library.biosample.subcellular_fraction_term_name!=nucleoplasm&award.rfa=ENCODE3&award.rfa=ENCODE4)

## Download processed files (TSV gene quantifications)
```bash
xargs -L 1 curl -O -J -L < files.txt
```
