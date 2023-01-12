# Analysis of human ENCODE3 and 4 short read RNA-seq data
* [ENCODE cart of all human data](https://www.encodeproject.org/carts/4ea7a43f-e564-4656-a0de-b09c92215e52/)
* [Query for human data](https://www.encodeproject.org/cart-report/?type=Experiment&cart=/carts/4ea7a43f-e564-4656-a0de-b09c92215e52/)

## Download processed files (TSV gene quantifications)
```bash
xargs -L 1 curl -O -J -L < files.txt
```
