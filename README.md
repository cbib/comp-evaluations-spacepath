# Computational validations of the SpacePath pipeline

## How to obtain summarized reporting of the results of the ablations vs the baseline SpacePath pipeline

Before starting the computational validations, the input folder must be populated. 
Specific instructions for the setup of the input folder are in `src/setup_input` folder.

The computational validations are performed by the scripts:
* `0_baseline.sh`
* `1_ablation_nofilter.sh`
* `2_ablation_random_id.sh`


Finally, the scripts that are located in `src/reporting` produce the 
tables and figures that summarize all the results for the Supplementary Methods
section of the paper.


------------------

## Appendix

Structure of this work:

```text
├── 0_baseline.sh
├── 1_ablation_nofilter.sh
├── 2_ablation_random_id.sh
├── Reporting.md
│
├── config
│   ├── baseline.tsv
│   ├── no_filter.tsv
│   ├── README.md
│   └── tmp
│       ├── no_filter_RESCUE.tsv
│       ├── rescue_random.tsv
│       └── the_data_list.tsv
│
├── input
│   ├── Bladder-cancer-hum
│   │   ├── Bladder-cancer-hum.h5ad
│   │   ├── Bladder-cancer-hum-HMDB.csv
│   │   └── Bladder-cancer-hum-LipidMaps.csv
│   ├── ...
│   └── UBcalibr4
│       ├── UBcalibr4.h5ad
│       ├── UBcalibr4-HMDB.csv
│       └── UBcalibr4-LipidMaps.csv
├── logs
├── output
│   ├── 0_baseline
│   │   ├── Bladder-cancer-hum
│   ├── ...
│   │   └── UBcalibr4
│   ├── 1_ablation_nofilter
│   │   ├── Bladder-cancer-hum
│   │   ├── Brain-center-rat
│   │   ├── ...
│   │   └── UBcalibr4
│   ├── 2_ablation_random_id
│   │   ├── Bladder-cancer-hum
│   │   ├── ...
│   │   └── UBcalibr4
│   └── output.md 
└── src
    ├── ablation1_nofilter.md
    ├── ablation2_random_id.py
    ├── bmhelpers.py
    ├── reporting
    │   ├── ambiguity_propag_score.py
    │   ├── numbers_sets_significant.pdf
    │   ├── redo_differential_keep_all.py
    │   └── run_visualizations.py
    └── setup_input
        ├── annex
        ├── setup_input.md
        └── setup_input.py

```
------------------
Author: 
Johanna Galvis


