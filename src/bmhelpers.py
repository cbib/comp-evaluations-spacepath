import os
import numpy as np
from typing import List, Tuple, Dict, Union
from scipy.sparse import csr_matrix
import pandas as pd
import anndata as ad
import pickle
import smpath.data as spd
import smpath.method as spm
import smpath.processing.measures_methods as measures_methods
import SpacePath.src.helpers as s_helpers



pathways_files_location_info = {
        "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG'},
        "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                          'subfolder': "metabolites-mappings",
                          'id_type': 'moleculeIds'}}


##

def compute_diff_analysis_for_sets(
        pathways_types: List[str],
        clustering_parameters: Dict[str, str | int | float],
        SAMPLE: str, subfolder: str, data_centralized_folder: str,
       ) -> None:  # before: n_top_per_cluster 25
    if not "FAMILIES" in pathways_types:
        pathways_types.append("FAMILIES")
    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    scores_merged = s_helpers.join_scores_adatas(pathways_types, SAMPLE, subfolder,
                                       data_centralized_folder)

    clusters_df = pd.read_csv(os.path.join(data_centralized_folder, subfolder,
                                           f'{SAMPLE}-clusters-{str_out}.txt'),
                              header=0)

    scores_merged =  s_helpers.safe_insert_clusters_column(scores_merged,
                                                clusters_df,
                                                str_out)
    #  shift to avoid negative values
    scores_merged.X = csr_matrix(
        scores_merged.X.todense() + abs(scores_merged.X.min().min()))

    diff_df = s_helpers.squidpy_differential_analysis(
        scores_merged, group_str=str_out)

    # Moran : compute for the sets
    scores_merged = measures_methods.compute_moran_i(scores_merged)
    scores_merged.var['I'] = scores_merged.uns['moranI'].loc[
        scores_merged.var.index.tolist(), 'I'].astype(float).tolist()

    diff_df = pd.merge(
        diff_df,
        scores_merged.var[['pathway_id', 'pathway_name', 'I']].copy(),
        how='left', left_on="names", right_on="pathway_id")

    diff_df.to_csv(os.path.join(data_centralized_folder, subfolder,
                                         f'{SAMPLE}-diff-results-SETS.tsv'),
                            sep='\t', header=True, index=False)
##


def random_adapt_half_per_mz_selector(annots_df, usable_id):
    """
    yields an annotations df with a half of the IDs (per m/z),
    randomly selected whenever it is possible

    """
    ius_list = []
    annots_df = annots_df.assign(
        ius=annots_df['mz'].astype(str) + "ww" + annots_df[usable_id])
    mzs_list = list(set(annots_df['mz'].astype(str).tolist()))
    for imz in mzs_list:
        here_ids_list = annots_df.loc[annots_df['mz'].astype(str) == imz,
        usable_id].unique().tolist()

        local_size = int(len(here_ids_list) / 2) if len(
            here_ids_list) >= 2 else 1
        tmps = np.random.choice(here_ids_list,
                                size=local_size, replace=False)
        for j in tmps:
            ius_list.append(f'{imz}ww{j}')

    annots_df = annots_df.loc[annots_df['ius'].astype(str).isin(ius_list),:].copy()
    annots_df = annots_df.drop(columns=['ius'])

    return annots_df


def random_id_per_mz_selector(annots_df, usable_id):
    ius_list = []
    annots_df = annots_df.assign(
        ius=annots_df['mz'].astype(str) + "ww" + annots_df[usable_id])
    mzs_list = list(set(annots_df['mz'].astype(str).tolist()))
    seen = []
    for imz in mzs_list:
        here_ids_list = annots_df.loc[annots_df['mz'].astype(str) == imz,
        usable_id].unique().tolist()
        here_ids_list = list(set(here_ids_list) - set(seen))
        try:
            tmp = np.random.choice(here_ids_list, size=1, replace=False)
        except ValueError:
            tmp = ['']
        ius_list.append(f'{imz}ww{tmp[0]}')   # always is only one element
        seen.append(tmp[0])

    annots_df = annots_df.loc[annots_df['ius'].astype(str).isin(ius_list),:].copy()
    annots_df = annots_df.drop(columns=['ius'])

    return annots_df


def ids_in_good_order(current_mzs, annots_df, usable_id):
    annots_df = annots_df.drop_duplicates(subset=['mz', usable_id])
    annots_df.index = annots_df['mz'].astype(str).tolist()
    tmp = annots_df.T
    ids_list_ordered = tmp.loc[usable_id, current_mzs].tolist()
    ko = 00

    return ids_list_ordered


def compute_pathway_scores_random_id(pathways_types,
                           min_pathway_size,
                           SAMPLE, subfolder,data_centralized_folder) -> None:
    """
        Compute pathways scores exclusively,
        saves each type of scores (e.g. KEGG or WIKI) to distinct .h5ad files
        """
    for pw_type in pathways_types:
        adata_mz, clean_annots_dict = s_helpers.safe_load_data_and_clean_annotations(
            SAMPLE, pathways_types, subfolder, data_centralized_folder
        )
        annots_df = clean_annots_dict[pw_type]
        pathways_obj = spd.MetabolicPathways(
            data_dir=data_centralized_folder,
            subfolder=pathways_files_location_info[pw_type]['subfolder'],
            file=pathways_files_location_info[pw_type]['file']
        )
        pathways_obj.molecule_id_type = pathways_files_location_info[pw_type][
            'id_type']

        pathways_obj.load()

        # -- the random id selection part; see the function  --
        annots_df = random_id_per_mz_selector(
            annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type']
        )
        pathways_obj, annots_df = spm.SyncPathwaysAndAnnots.run(
            pathways_obj.model_copy(),
            annots_df,
            usable_id=pathways_files_location_info[pw_type]['id_type'],
            min_set_size=min_pathway_size)

        #   -- ablation : no weights --
        # metabolites_computed_df, coords = spm.ComputeAssignmentsWeights.run(
        #     adata_mz, annots_df,
        #     usable_id=pathways_files_location_info[pw_type]['id_type'])

        coords = adata_mz.obsm['spatial']

        metabolites_computed_df = adata_mz[:, annots_df['mz'].astype(str).tolist()].to_df()
        del adata_mz  # free some space
        usable_id = pathways_files_location_info[pw_type]['id_type']
        metabolites_computed_df.columns = annots_df[usable_id].tolist()
        print(metabolites_computed_df.shape)
        if metabolites_computed_df.shape[1] == 0:
            metabolites_computed_df['nope'] = np.repeat(0, metabolites_computed_df.shape[0])

        scores_adata, time_elapsed = s_helpers.run_pathway_scores_steps_produce_anndata(
            metabolites_computed_df, coords, pathways_obj,
            min_set_size=min_pathway_size,
            type_abbr=pw_type)


        out_file = os.path.join(data_centralized_folder, subfolder,
                                f'{SAMPLE}-scores-{pw_type}.h5ad')
        scores_adata.write_h5ad(out_file)
        print("scores saved to file %s", out_file)


def compute_molecules_families_random_id(cutoff_tree,
                               SAMPLE,
                               subfolder,
                               data_centralized_folder,
                               plot_folder
                               )->None:
    np.random.seed(42)
    top_differential_mzs = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-top-differential-mzs.tsv'), sep='\t',
        header=0)

    # to ensure consistency of families,  refilter for 25 top mz per cluster
    top_differential = s_helpers.differential_results_filter(top_differential_mzs,
                                                   n_top_per_cluster=25)

    clean_annots_nokegg_ids = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    # -- ablation : random id
    clean_annots_nokegg_ids = random_id_per_mz_selector(
        clean_annots_nokegg_ids, usable_id='moleculeIds')

    molecules_ids_df = s_helpers.prepare_molecular_ids_df(
        top_differential, clean_ann_df=clean_annots_nokegg_ids)

    mol_similarities_result = spm.ComputeMoleculesSimilarityMatrix.run(
        molecules_ids_list=molecules_ids_df["moleculeIds"])

    with open(os.path.join(data_centralized_folder, subfolder,
                           f"{SAMPLE}-mol-similarities-matrix.pkl"),
              "wb") as f:
        pickle.dump(mol_similarities_result, f)

    updated_config = spm.ComputeMoleculesFamiliesConfig(
        cutoff_tree=cutoff_tree)
    computed_fam_results = spm.ComputeMoleculesFamilies.run(
        mol_similarities_result,
        molecules_ids_df,
        SAMPLE,
        subfolder,
        data_centralized_folder,
        plot_folder,
        config=updated_config)  # all results: maybe of relevance for the user

    with open(os.path.join(
            data_centralized_folder, subfolder, f"{SAMPLE}-mol-fam-dict.pkl"),
            "wb") as f:
        pickle.dump(computed_fam_results, f)


def handle_error_families_scores(
                            sample,
                            subfolder,
                               data_centralized_folder ):
    """
    """
    adata_mz = ad.read_h5ad(os.path.join(data_centralized_folder, subfolder,
                                         f'{sample}.h5ad' ))
    print(adata_mz.obsm)
    coords = adata_mz.obsm['spatial']
    sets_ids_df = pd.DataFrame({'pathway_id': ['none'], 'pathway_name': ['none'] })
    piled = pd.DataFrame({'pathway_id': ['none'], 'mz': ['none'],
                          'moleculeIds': ['none'] })

    scores_df = pd.DataFrame({
        'none': np.repeat(np.nan, adata_mz.n_obs)
    }, index=adata_mz.obs.index)

    scores_adata = s_helpers.score_as_adata(scores_df,
                                  sets_ids_df,
                                  piled,
                                  coords)
    scores_adata.write_h5ad(os.path.join(data_centralized_folder,
                                         subfolder, f'{sample}-scores-FAMILIES.h5ad'))

