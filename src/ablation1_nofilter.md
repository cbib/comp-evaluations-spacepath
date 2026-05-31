# Ablation without filtering features: "no filter"

The script in `../1_ablation_nofilter.sh` calls SpacePath
and uses the parameters of `../config/no_filter.tsv` that
are tuned to avoid the filtering of the values, and that for all the tissues
(datasets). These are the thresholds:
* offSample_cutoff : 1.1
* msm_cutoff: -1.1
* fdr_cutoff: 1.1
* Moran_I_cutoff: -1.1
* Moran_padj_cutoff: 1.1

Expected values and logics of the thresholds in the pipeline:

* offSample : interval (0,1), the values below cutoff are retained
* msm : interval (0,1), the values above cutoff are retained
* fdr : interval (0,1), the values below cutoff are retained
* Moran_I : interval (-1,1), the values above cutoff are retained
* Moran_padj : interval (0,1), the values below cutoff are retained
