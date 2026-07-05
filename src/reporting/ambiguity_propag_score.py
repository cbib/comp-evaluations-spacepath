# johaGL2026
# ambiguity propagation score


import os
from typing import Dict, List

import smpath.processing as smp
from SpacePath import api as sap
from SpacePath.src.orchestrator import DatasetsProcessor as DP
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import SpacePath.src.helpers as s_helpers
import squidpy as sq
import seaborn as sns
import re


def compute_differential_keep_all(pathways_types: List[str],
            clustering_parameters: Dict[str, str | int | float],
            SAMPLE: str, subfolder: str, data_centralized_folder: str):
    """modified version of differential to use only in ablations analysis"""

    str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"

    scores_merged = s_helpers.join_scores_adatas(pathways_types, SAMPLE, subfolder,
                                       data_centralized_folder)

    clusters_df = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-clusters-{str_out}.txt'),
        header=0)

    scores_merged = s_helpers.safe_insert_clusters_column(scores_merged,
                                                clusters_df,
                                                str_out)
    #  shift to avoid negative values
    scores_merged.X = csr_matrix(
        scores_merged.X.todense() + abs(scores_merged.X.min().min()))

    differential_results_by_cluster = s_helpers.squidpy_differential_analysis(
        scores_merged, group_str=str_out)

    scores_merged = smp.measures_methods.compute_moran_i(scores_merged,
                                                         n_perms=2) # fast comput, pvalues not required for this viz
    scores_merged.var['I'] = scores_merged.uns['moranI'].loc[
        scores_merged.var.index.tolist(), 'I'].astype(float).tolist()


    out_diff_all = pd.merge(
        differential_results_by_cluster,
        scores_merged.var[['pathway_id', 'pathway_name', 'I']].copy(),
        how='left', left_on="names", right_on="pathway_id")
    return out_diff_all



def extract_contracted_name(lipid_name):
    """
    Transforms full name into contracted species expression

    :param lipid_name: example: "PE(O-16:1(1Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))"
    :returns: contracted species name, example "PE_38:7"
    """
    chains = re.findall(r'(\d+:\d+)', lipid_name)
    head_match =  re.match(r'^([^(]+)', lipid_name)
    try:
        polar_head = head_match.group(1)
    except:
        polar_head = "UNDETERMINED"

    try:
        carbon_n_list = []
        insaturations_n_list = []
        for chain in chains:
            carbon_n_list.append(chain.split(':')[0])
            insaturations_n_list.append(chain.split(':')[1])

        carbon_sum = np.sum(np.array(carbon_n_list, dtype=np.int32))
        insaturations_sum = np.sum(
            np.array(insaturations_n_list, dtype=np.int32))
    except:
        carbon_sum = ''
        insaturations_sum = ''

    return f'{polar_head}_{carbon_sum}:{insaturations_sum}'


def compute_lipid_species(names_list: List[str]):
    """transforms a list of full lipid names into a list of unique species"""
    lipid_species = set()
    if len(names_list) > 0:
        for name in names_list:
            lipid_species = lipid_species.union(
                {extract_contracted_name(name)})

    return list(lipid_species)


# def extract_lipid_class(full_lipid_name:str):
#     # Extract the substring before the first '('
#     match = re.match(r'^([^(]+)', full_lipid_name)
#
#     return match.group(1) if match else "UNDETERMINED"


# def compute_lipid_classes(names_list: List[str]):
#     lipid_classes = set()
#     if len(names_list) > 0:
#         for name in names_list:
#             lipid_classes = lipid_classes.union({extract_lipid_class(name)})
#
#     return list(lipid_classes)


def mappings_only_families(mapps_df: pd.DataFrame)->pd.DataFrame:
    mapps_df = mapps_df.loc[
        ~mapps_df['pathway_id'].str.startswith("hsa"), :].copy()
    mapps_df = mapps_df.loc[
        ~mapps_df['pathway_id'].str.startswith("WP"), :].copy()
    mapps_df = mapps_df.dropna(subset=["moleculeNames"])

    return mapps_df


def family_molecule_mz_ambiguity_propagation_score(mapps_df: pd.DataFrame, i_fam: str, remaining_fams):
    remain_df = mapps_df.loc[mapps_df["pathway_id"].isin(remaining_fams),:].copy()
    i_df = mapps_df.loc[mapps_df["pathway_id"] == i_fam,:].copy()

    ambigs_scores_mzs = []

    for i_mz in i_df["mz"].unique().tolist():

        remain_f_mz = remain_df.loc[
            remain_df["mz"].astype(str) == str(i_mz), :].copy()
        i_f_mz = i_df.loc[i_df["mz"].astype(str) == str(i_mz), :].copy()

        # all_ids__i_mz = list(
        #     set(remain_f_mz["moleculeNames"].unique().tolist()).union(
        #         set(i_f_mz["moleculeNames"].unique().tolist())
        #     )
        # )

        # species_all__i_mz = compute_lipid_species(all_ids__i_mz)
        species_in = compute_lipid_species(i_f_mz["moleculeNames"].unique().tolist())
        species_out = compute_lipid_species(remain_f_mz["moleculeNames"].unique().tolist())
        # species_union = set(species_in).union(set(species_out))
        species_common = set(species_in).intersection(set(species_out)) # this penalizes the score

        ambigs_scores_mzs.append(
            len(species_common) + (len(species_out)/len(species_in))
        )
        
    return np.array(ambigs_scores_mzs).mean()


def compute_ambig_propag_across_families( SAMPLE,
                                         subfolder, data_centralized_folder):

    scores_merged = s_helpers.join_scores_adatas(['FAMILIES'],
                                                 SAMPLE,
                                                 subfolder,
                                                 data_centralized_folder)

    mapps_df = scores_merged.uns['current_pathways_df'].copy()

    clean_annots_nokegg_ids = pd.read_csv(
        os.path.join(data_centralized_folder, subfolder,
                     f'{SAMPLE}-WIKI-CleanAnnotations.tsv'),
        sep='\t', header=0, index_col=None)

    mapps_df['mz'] = mapps_df['mz'].astype(str)
    clean_annots_nokegg_ids['mz'] = clean_annots_nokegg_ids['mz'].astype(str)
    mapps_df = pd.merge(mapps_df, clean_annots_nokegg_ids, how="left",
                        on=["moleculeIds", "mz"])

    # discard all rows with  pathways (leads to errors):
    mapps_df = mappings_only_families(mapps_df)

    # mapps_df = mapps_df.loc[mapps_df['database'] == "HMDB", :].copy()
    # mzs_uniq = mapps_df['mz'].astype(str).unique().tolist()

    fams_uniq = [h for h in mapps_df['pathway_id'].astype(str).unique().tolist()
                 if h.startswith("MolFam")]

    ambig_propag_score_per_fam = []
    for i_fam in fams_uniq:
        remaining_fams = list(set(fams_uniq) - {i_fam})
        i_ambig_propag_score = family_molecule_mz_ambiguity_propagation_score(
            mapps_df, i_fam, remaining_fams)
        # print(i_fam, "<<<<", i_ambig_propag_score)
        ambig_propag_score_per_fam.append(i_ambig_propag_score)

    ambig_propag_score_df = pd.DataFrame(
        {"pathway_id": fams_uniq,
         "ambig_propag_score": ambig_propag_score_per_fam})

    return ambig_propag_score_df