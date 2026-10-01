#!/usr/bin/env python
# johaGL
"""
## Ablation 2 ##
Runs SpacePath with a special modification in this code:
 selecting a random compound ID for pathways and families computation
"""

import os
import sys
import SpacePath.core.data as spd
from SpacePath import api as sap
from SpacePath.src.orchestrator import DatasetsProcessor as DP
import numpy as np
import pandas as pd
import bmhelpers as bmhelp


np.random.seed(42)

params_file = sys.argv[1]
centralized_data_dir = sys.argv[2]
plot_folder = sys.argv[3]

params_df = pd.read_csv(os.path.join(
    ".", params_file
), sep="\t", index_col=None, header=0)

print(params_file, '\n', centralized_data_dir, '\n', plot_folder)


for _, row in params_df.iterrows():
    (sample, annotations_dict, offSample_cutoff, msm_cutoff,
     fdr_cutoff, moran_I_cutoff,
     moran_padj_cutoff, clustering_parameters, pathways_types,
     min_pathway_size, cutoff_tree,
     min_family_size) = DP.extract_params_per_dataset(row)

    subfolder = sample

    annots_obj = spd.MetabolomeAnnotations(
        name=sample, data_dir=centralized_data_dir,
        subfolder=subfolder,
        annotation_files=annotations_dict)
    annots_obj.load()  # load annotations with library

    sap.normalize_classify_and_save_anndata_mz(
        annots_obj,
        offSample_cutoff=offSample_cutoff,
        msm_cutoff=msm_cutoff, fdr_cutoff=fdr_cutoff,
        moran_I_cutoff=moran_I_cutoff,
        moran_padj_cutoff=moran_padj_cutoff,
        SAMPLE=sample,
        subfolder=subfolder,
        data_centralized_folder=centralized_data_dir,
        plot_folder=os.path.join(centralized_data_dir,
                                 subfolder)
    )

    sap.clean_annotations_and_save_tsv(
        annots_obj, pathways_types,
        sample,
        subfolder=subfolder,
        data_centralized_folder=centralized_data_dir)


    sap.compute_clustering_of_spots(
        pathways_types, clustering_parameters,
        sample,  subfolder=subfolder,
        data_centralized_folder=centralized_data_dir,
        plot_folder=os.path.join(centralized_data_dir,subfolder))


    sap.compute_differential_analysis_of_mzs(pathways_types,
                                         clustering_parameters,
                                         sample,
        subfolder=subfolder,
        data_centralized_folder=centralized_data_dir
                                         )

    bmhelp.compute_pathway_scores_random_id(pathways_types,
                                            min_pathway_size,
                                            sample,
                                            subfolder=subfolder,
                                            data_centralized_folder=centralized_data_dir)

    try:
        bmhelp.compute_molecules_families_random_id(cutoff_tree,
                                                sample,
                                                subfolder=subfolder,
                                                data_centralized_folder=centralized_data_dir,
                                                plot_folder=os.path.join(centralized_data_dir,
                                 subfolder) )
    except ValueError:
        pass

    try:
        sap.compute_families_scores(min_family_size,
                            sample,
                            subfolder=subfolder,
                               data_centralized_folder=centralized_data_dir )
    except Exception as e: # no IDs satisfying condition for families for scoring
        bmhelp.handle_error_families_scores( sample,
                            subfolder=subfolder,
                               data_centralized_folder=centralized_data_dir)

    sap.compute_differential_families_and_pathways(
        pathways_types,
        clustering_parameters, sample,
        subfolder=subfolder,
        data_centralized_folder=centralized_data_dir)
