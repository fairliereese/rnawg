```bash
annot=gencode_v29
build=hg38
db=../talon/human.db

talon_create_GTF --db ${db} \
    -a ${annot} \
    -b ${build} \
    --whitelist human_talon_swan_pass_list.csv \
    --o human_swan
```