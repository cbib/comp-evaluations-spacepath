import os
import smpath.data as spd
from SpacePath import api as sap
from SpacePath.src.orchestrator import DatasetsProcessor as DP
import numpy as np
import pandas as pd
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from ambiguity_propag_score import compute_ambiguity_propag_score
import SpacePath.src.helpers as s_helpers


plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8})


np.random.seed(42)


def save_plot_metric_tissue_specific(df_piled: pd.DataFrame,
                                     metric_to_plot="logfoldchanges", out_file=''):
    g = sns.catplot(
        data=df_piled,
        x="eval_type",
        y=metric_to_plot,
        hue="set_type",
        row="tissue_name",
        kind="box",
        sharey=True,
        height=3,
        aspect=2,
        palette="Pastel2"
    )

    for tissue, ax in zip(g.row_names, g.axes.flat):
        sns.stripplot(
            data=df_piled[df_piled["tissue_name"] == tissue],
            x="eval_type",
            y=metric_to_plot,
            hue="set_type",
            dodge=True,
            size=2,
            palette="Set2",
            ax=ax
        )
        # Remove duplicate legends from stripplot
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    for ax in g.axes.flat:
        xticks = ax.get_xticks()

        for x in xticks[:-1]:
            ax.axvline(
                x + 0.5,
                color="0.8",
                linestyle=":",
                linewidth=0.8,
                zorder=0
            )

    g.figure.set_size_inches(4, 5)  # width=8", height=10"
    if out_file != '':
        g.figure.savefig(
            out_file,
            dpi=150,
            bbox_inches="tight"
        )
    else:
        print("!! missing out file for figure !!")


def save_plot_metric_global(df_piled: pd.DataFrame,
                            metric_to_plot="logfoldchanges",
                            out_file=''):
    """all values together all tissues mixed"""
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 3))
    sns.boxplot(data=df_piled, x="eval_type", y=metric_to_plot,
                     hue="set_type", fill=True, gap=.1, palette="Pastel2",
                ax=ax)
    sns.stripplot(data=df_piled, x="eval_type", y=metric_to_plot,
                  hue="set_type", s=2, palette="Set2", dodge=True,
                  legend=False,
                  ax=ax)
    if out_file != '':
        plt.savefig(
            out_file,
            dpi=150,
            bbox_inches="tight"
        )
    else:
        print("!! missing out file for figure !!")


# --
out_dir ="../../output"
tissues = ["Brain-center-rat", "Stomach-rat"]
ablations_out_dirs = {'baseline': '0_baseline',
                      'no filter': '1_ablation_nofilter',
                      'random id': '2_ablation_random_id',
                      'pathways only': '3_ablation_pathways_alone'}

dfs_metrics_I_lfc = list()
dfs_metric_new = list()
for tissue in tissues:
    for ablation_type in ablations_out_dirs.keys():
        ablation_dir = ablations_out_dirs[ablation_type]
        df_top_sets = pd.read_csv(os.path.join(out_dir, ablation_dir, tissue,
                                    f"{tissue}-top-differential-SETS.tsv"),
                       sep="\t", index_col=None, header=0)

        df_top_sets['clus_set'] = df_top_sets['group'].astype(str) + "_" + df_top_sets[
            'pathway_id']
        df_top_sets = df_top_sets[
            ['clus_set', 'pathway_id', 'logfoldchanges', 'pvals', 'pvals_adj', 'I']]
        df_top_sets['eval_type'] = ablation_type
        df_top_sets['set_type'] = ['ghost' if x.startswith('MolFam')
                            else 'pathway' for x in
                            df_top_sets['pathway_id'].tolist()]

        if ablation_dir == '3_ablation_pathways_alone':
            df_top_sets = df_top_sets.loc[
                df_top_sets['set_type'] == 'pathway',:].copy()

        df_top_sets['tissue_name'] = tissue
        dfs_metrics_I_lfc.append(df_top_sets)
        del df_top_sets

        # * ambiguity propagation score
        if ablation_dir != '3_ablation_pathways_alone':
            pathways_types = [ 'FAMILIES']
            df_aps= compute_ambiguity_propag_score(
                pathways_types, tissue,
                subfolder=os.path.join(ablation_dir, tissue),
                data_centralized_folder=out_dir)

            df_aps['eval_type'] = ablation_type
            df_aps['set_type'] = pd.Categorical(np.repeat("ghost", df_aps.shape[0]),
                                                categories=["ghost", "pathway"])
            df_aps['tissue_name'] = tissue
            dfs_metric_new.append(df_aps)


df_piled = pd.concat(dfs_metrics_I_lfc)
print( df_piled.columns)
df_ambig_propag_score = pd.concat(dfs_metric_new)
# ------ plot metric 1

# save_plot_metric_global(df_piled,  metric_to_plot="logfoldchanges",
#                         out_file="metric1-mix.png")
# save_plot_metric_tissue_specific(df_piled, metric_to_plot="logfoldchanges",
#                                  out_file="metric1.png")
#
# # ------ plot metric 2
# save_plot_metric_global(df_piled,  metric_to_plot="I",
#                         out_file="metric2-mix.png")
# save_plot_metric_tissue_specific(df_piled, metric_to_plot="I",
#                                  out_file="metric2.png")


# ------- plot metric 3

save_plot_metric_global(df_ambig_propag_score,
                        metric_to_plot="ambig_propag_score",
                        out_file="metric3-mix.png")
save_plot_metric_tissue_specific(df_ambig_propag_score,
                                 metric_to_plot="ambig_propag_score",
                                out_file="metric3.png")

