# Ablation without filtering features: "no filter"

The script in `../1_ablation_nofilter.sh` calls SpacePath
and uses the parameters of `../config/no_filter.tsv` that
are tuned to avoid the filtering of the ions ($m/z$ 's), and that for all the tissues
(datasets). These are the thresholds fixed for the _"no filter"_ ablation:
* `offSample_cutoff` : 1.1
* `msm_cutoff`: -1.1
* `fdr_cutoff`: 1.1
* `Moran_I_cutoff`: -1.1
* `Moran_padj_cutoff`: 1.1

And these values are the same for all the datasets (tissues).

### How these cutoffs work
Explanation of typical intervals of values for each metric,
and logics how thresholds are used by the pipeline:

* _offSample_ : interval (0,1), the ions with values below `cutoff` are retained
* _msm_ : interval (0,1), the ions with values above `cutoff` are retained
* _fdr_ : interval (0,1), the ions with values below `cutoff` are retained
* _Moran_I_ : interval (-1,1), the ions with values above `cutoff` are retained
* _Moran_padj_ : interval (0,1), the ions with values below `cutoff` are retained
