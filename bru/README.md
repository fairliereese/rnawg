# Analysis of human ENCODE4 Bru data
* [Query for data](https://www.encodeproject.org/report/?type=Experiment&control_type!=*&status=released&replicates.library.biosample.biosample_ontology.term_id=EFO:0002106&replicates.library.biosample.biosample_ontology.term_id=EFO:0001203&replicates.library.biosample.biosample_ontology.term_id=EFO:0006711&replicates.library.biosample.biosample_ontology.term_id=EFO:0002713&replicates.library.biosample.biosample_ontology.term_id=EFO:0002847&replicates.library.biosample.biosample_ontology.term_id=EFO:0002074&replicates.library.biosample.biosample_ontology.term_id=EFO:0001200&replicates.library.biosample.biosample_ontology.term_id=EFO:0009747&replicates.library.biosample.biosample_ontology.term_id=EFO:0002824&replicates.library.biosample.biosample_ontology.term_id=CL:0002327&replicates.library.biosample.biosample_ontology.term_id=CL:0002618&replicates.library.biosample.biosample_ontology.term_id=EFO:0002784&replicates.library.biosample.biosample_ontology.term_id=EFO:0001196&replicates.library.biosample.biosample_ontology.term_id=EFO:0001187&replicates.library.biosample.biosample_ontology.term_id=EFO:0002067&replicates.library.biosample.biosample_ontology.term_id=EFO:0001099&replicates.library.biosample.biosample_ontology.term_id=EFO:0002819&replicates.library.biosample.biosample_ontology.term_id=EFO:0009318&replicates.library.biosample.biosample_ontology.term_id=EFO:0001086&replicates.library.biosample.biosample_ontology.term_id=EFO:0007950&replicates.library.biosample.biosample_ontology.term_id=EFO:0003045&replicates.library.biosample.biosample_ontology.term_id=EFO:0003042&replicates.library.biosample.internal_tags=Deeply+Profiled&assay_title=long+read+RNA-seq)

## Get list of dataset names and detected transcript IDs from the LR bulk Ljungman data
```bash
python get_ljungman_tids.py
```

## Generate all annotated + detected novel GTF from just the LR bulk Ljungman data
```bash
db=../lr_bulk/talon/human.db
annot=gencode_v29
build=hg38
opref=ljungman
talon_create_GTF \
    --db ${db} \
    -a ${annot} \
    -b ${build} \
    --whitelist ${opref}_complete_pass_list.csv \
    --o ${opref}_known_nic_nnc
```

## Generate detected novel GTF from just the LR bulk Ljungman data
```bash
db=../lr_bulk/talon/human.db
annot=gencode_v29
build=hg38
opref=ljungman
talon_create_GTF \
    --db ${db} \
    -a ${annot} \
    -b ${build} \
    --whitelist ${opref}_novel_pass_list.csv \
    --o ${opref}_nic_nnc
```
