import os
import smpath.data as spd
from SpacePath import api as sap
from SpacePath.src.orchestrator import DatasetsProcessor as DP
import numpy as np
import pandas as pd
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from ambiguity_propag_score import (compute_ambig_propag_across_families,
                                    compute_differential_keep_all)
import SpacePath.src.helpers as s_helpers


plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8})

np.random.seed(42)


out_dir ="../../output"
#tissues = ["Brain-center-rat", "Stomach-rat"]
mok = pd.read_csv("../../config/no_filter.tsv", sep="\t", index_col=None, header=0)
tissues = mok['dataset_name'].unique().tolist()

ablations_out_dirs = {'baseline': '0_baseline',
                      'no filter': '1_ablation_nofilter',
                      'random id': '2_ablation_random_id',
                      'pathways only': '3_ablation_pathways_alone'}

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