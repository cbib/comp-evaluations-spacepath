import os
import SpacePath.core.data as spd
from SpacePath import api as sap
from SpacePath.src.orchestrator import DatasetsProcessor as DP
import numpy as np
import pandas as pd
from ambiguity_propag_score import compute_differential_keep_all
import SpacePath.src.helpers as s_helpers

"""
Run this script in the first time to generate the "extended-differential-SETS.tsv"
files per dataset, as we require unfiltered differential results
"""


np.random.seed(42)


out_dir ="../../output"
#tissues = ["Brain-center-rat", "Stomach-rat"]
mok = pd.read_csv("../../config/no_filter.tsv", sep="\t", index_col=None, header=0)
tissues = mok['dataset_name'].unique().tolist()

ablations_out_dirs = {'baseline': '0_baseline',
                      'no filter': '1_ablation_nofilter',
                      'random id': '2_ablation_random_id'}

dfs_metrics_I_lfc = list()
dfs_metric_new = list()
for tissue in tissues:
    for ablation_type in ablations_out_dirs.keys():
        cl_pr = mok.loc[
            mok['dataset_name'] == tissue, 'cluster_params'].values[0].split(':')
        clustering_parameters = {'choice': cl_pr[0].lower(), 'numeric_parameter': cl_pr[1]}
        print(clustering_parameters)

        ablation_dir = ablations_out_dirs[ablation_type]

        pathways_types = ['KEGG', 'WIKI', 'FAMILIES']
        if ablation_type == 'pathways_only':
            pathways_types = ['KEGG', 'WIKI']

        nest_dir = os.path.join(out_dir, ablation_dir)
        df_sets_res = compute_differential_keep_all(
            pathways_types,
            clustering_parameters, tissue, tissue, nest_dir)
        df_sets_res.to_csv(os.path.join(out_dir, ablation_dir, tissue,
                                  f"{tissue}-extended-differential-SETS.tsv"),
                       sep="\t", index=False, header=True)