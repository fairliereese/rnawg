```bash
conda activate encode_submissions
```

## Documents
```bash
eu_register.py -m dev -p document -i documents.tsv
# eu_register.py -m prod -p document -i documents.tsv
```

## Reference files
```bash
eu_register.py -m dev -p file -i file_refs.tsv
# eu_register.py -m prod -p file -i file_refs.tsv
```

## CAGE, RAMPAGE, PAS
```bash
eu_register.py -m prod -p file -i cage.tsv
eu_register.py -m prod -p file -i rampage.tsv
eu_register.py -m prod -p file -i pas.tsv
```
