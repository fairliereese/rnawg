```bash
conda activate encode_submissions
```

## Documents
```bash
eu_register.py -m prod -p document -i documents.tsv # done
```

## Reference files
```bash
eu_register.py -m prod -p file -i file_refs.tsv # done
```

## CAGE, RAMPAGE, PAS
```bash
eu_register.py -m prod -p file -i cage.tsv # done
eu_register.py -m prod -p file -i rampage.tsv # done
eu_register.py -m prod -p file -i pas.tsv # done
```

## Human LR data
```bash
eu_register.py -m prod -p file -i lr_human.tsv # done
```

## Mouse LR data
```bash
eu_register.py -m prod -p file -i lr_human.tsv # done
```

## Document patches
```bash
eu_register.py -m prod -p reference -i ref_doc_patch.tsv --patch -w # done
eu_register.py -m prod -p annotation -i annot_doc_patch.tsv --patch -w # done
```

## File patches
```bash
eu_register.py -m prod -p file -i file_patch.tsv --patch -w # done
eu_register.py -m prod -p file -i file_patch_2.tsv --patch -w # done
eu_register.py -m prod -p file -i file_reupload.tsv --patch -w # donehpe

```

```python
import encode_utils as eu
from encode_utils.connection import Connection
conn = Connection("www.encodeproject.org")
file_ids = ['ali-mortazavi:rnawg_mouse_lapa_tes_minus_norm_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tes_plus_norm_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tes_minus_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tes_plus_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tss_minus_norm_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tss_plus_norm_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tss_minus_signal',
            'ali-mortazavi:rnawg_mouse_lapa_tss_plus_signal']
paths = ['/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tes/ratio/all_polyA_ratio_neg.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tes/ratio/all_polyA_ratio_pos.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tes/counts/all_polyA_counts_neg.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tes/counts/all_polyA_counts_pos.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tss/ratio/all_tss_ratio_neg.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tss/ratio/all_tss_ratio_pos.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tss/counts/all_tss_counts_neg.bw',
         '/share/crsp/lab/seyedam/share/mousewg/data/results/lapa/tss/counts/all_tss_counts_pos.bw']      
for fid, path in zip(file_ids, paths):
  if '_tes_' in fid:
    assert 'polyA' in path
  if '_tss_' in fid:
    assert '_tss_' in path
  if '_plus_' in fid:
    assert '_pos' in path
  if '_minus_' in fid:
    assert '_neg' in path
  if not '_norm_' in fid:
    assert '_counts_' in path
  if '_norm_' in fid:
    assert '_ratio_' in path

for fid, path in zip(file_ids, paths):
  conn.upload_file(file_id=fid, file_path=path)
```

Remove the dumb chromosomes from the mouse polyA atlas
```python
import pyranges as pr
fname = '/dfs7/samlab/mcelik/mousewg/data/resources/polyasite_atlas/atlas.clusters_formatted.bed'
ofile = '/dfs7/samlab/mcelik/mousewg/data/resources/polyasite_atlas/atlas.clusters_formatted_canon_only.bed'

df = pr.read_bed(fname).as_df()
df = df.loc[~df.Chromosome.str.contains('JH')]

df = pr.PyRanges(df)
df.to_bed(ofile)
```

```bash
gzip -c /dfs7/samlab/mcelik/mousewg/data/resources/polyasite_atlas/atlas.clusters_formatted_canon_only.bed > /dfs7/samlab/mcelik/mousewg/data/resources/polyasite_atlas/atlas.clusters_formatted_canon_only.bed.gz
```
