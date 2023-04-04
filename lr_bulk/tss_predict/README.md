## Machine learning models predict the support for long-read TSS peaks by other TSS-annotating assays and in a cross-cell type manner

This section contains the code to train/test logistic regression models to predict the overlap (support) for a peak in different TSS assays
The overall order for preprocessing, training and testing the logit models:

```
./combine_all_chrs_jamboree_files.sh
Rscript combine_cerberus_replicates.R
./preproc_files.sh
./make_labels.sh
./bigwigavg_DHS.sh
./split_train_test.sh

cd training/
Rscript train_diff_models_AIC.R
Rscript train_glm.R
Rscript predict_glm.R
Rscript predict_glm_diff_cell_line.R

```
