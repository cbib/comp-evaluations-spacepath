import os
import numpy as np
import pandas as pd
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from ambiguity_propag_score import compute_ambig_propag_across_families

np.random.seed(42)

"""
Run after redo_differential_keep_all.py
"""


plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9})

palette_set_type_dark = {"ghost" : "#fc8d62",   "pathway" : "#66c2a5"}
palette_set_type_pastel = {"ghost" : "#FDCDAC",   # pastel orange/coral
                           "pathway": "#B3E2CD"}  # pastel teal

height_subplots = 38  # 38 if all  # 28 if selection
width_subplots = 4  #
dpi_png_plots = 220


def save_plot_metric_tissue_specific(df_piled: pd.DataFrame,
                                     metric_to_plot="logfoldchanges", out_file='',
                                     sharey_bool=True, dodge_bool=True):
    g = sns.catplot(
        data=df_piled,
        x="eval_type",
        y=metric_to_plot,
        hue="set_type",
        row="dataset_name",
        kind="box",
        sharey=sharey_bool,
        sharex=True,
        height=1.5,
        aspect=2.5,
        palette=palette_set_type_pastel
    )

    for tissue, ax in zip(g.row_names, g.axes.flat):
        sns.stripplot(
            data=df_piled[df_piled["dataset_name"] == tissue],
            x="eval_type",
            y=metric_to_plot,
            hue="set_type",
            dodge=dodge_bool,  # stripplot: ambiguity score only has 1 class
            size=1.5,
            palette=palette_set_type_dark,
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

    g.figure.set_size_inches(w=width_subplots, h=height_subplots)  # width=8", height=10"
    if out_file != '':
        g.figure.savefig(
            out_file,
            dpi=dpi_png_plots,
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
                     hue="set_type", fill=True, gap=.1, palette=palette_set_type_pastel,
                ax=ax)
    sns.stripplot(data=df_piled, x="eval_type", y=metric_to_plot,
                  hue="set_type", s=1, palette=palette_set_type_dark,
                  dodge=True,
                  legend=False,
                  alpha=0.6,
                  ax=ax)
    if out_file != '':
        plt.savefig(
            out_file,
            dpi=dpi_png_plots,
            bbox_inches="tight"
        )
    else:
        print("!! missing out file for figure !!")


def categorize_specific_columns(df, datasetlist):
    df = df.assign(
        eval_type=pd.Categorical(
            df['eval_type'].tolist(),
            categories=['baseline', 'no filter', 'random id']))
    df = df.assign(
        set_type=pd.Categorical(df['set_type'].tolist(),
                                categories=['ghost', 'pathway']))
    df = df.assign(
        dataset_name=pd.Categorical(df['dataset_name'].tolist(),
                                categories=datasetlist))

    return df


# --
out_dir ="../../output"
#tissues = ["Brain-center-rat", "Stomach-rat"]
mok = pd.read_csv("../../config/no_filter.tsv", sep="\t", index_col=None, header=0)

#mok = mok.sample(12).copy()
tissues = mok['dataset_name'].unique().tolist()

only_selection = False  # True only if small demo   (e.g. slides paste)
if only_selection:
    tissues = ['Brain-center-rat','Brain-trv-mouse','Kidney-xenogr-2',
               'Lung-mouse', 'Liver-brca-mouse', 'Bladder-cancer-hum',
               'UBcalibr4',  'UBcalibr1', 'Brain-healthy-mouse',
               'Brain-demyel-mouse', 'Stomach-rat']


ablations_out_dirs = {'baseline': '0_baseline',
                      'no filter': '1_ablation_nofilter',
                      'random id': '2_ablation_random_id'
                      }

dfs_metrics_I_lfc = list()
dfs_metric_new = list()
for tissue in tissues:
    for ablation_type in ablations_out_dirs.keys():
        cl_pr = mok.loc[
            mok['dataset_name'] == tissue, 'cluster_params'].values[0].split(':')
        clustering_parameters = {'choice': cl_pr[0], 'numeric_parameter': cl_pr[1]}
        ablation_dir = ablations_out_dirs[ablation_type]

        df_sets_res = pd.read_csv(os.path.join(out_dir, ablation_dir, tissue,
                                 f"{tissue}-extended-differential-SETS.tsv"),
                       sep="\t", index_col=None, header=0)

        if df_sets_res.shape[0] == 0:
            df_sets_res = pd.DataFrame({'pathway_id': ['WPX'],
                                        'group' : ['1'],
                                        'logfoldchanges': [0],
                                        'pvals': [1], 'pvals_adj': [1],
                                        'I': [0]})

        df_sets_res['clus_set'] = df_sets_res['group'].astype(str) + "_" + df_sets_res[
            'pathway_id']
        df_sets_res = df_sets_res[
            ['clus_set', 'pathway_id', 'logfoldchanges', 'pvals', 'pvals_adj', 'I']]
        df_sets_res['eval_type'] = ablation_type
        df_sets_res['set_type'] = ['ghost' if x.startswith('MolFam')
                            else 'pathway' for x in
                            df_sets_res['pathway_id'].tolist()]

        # if ablation_dir == '3_ablation_pathways_alone':
        #     df_sets_res = df_sets_res.loc[
        #         df_sets_res['set_type'] == 'pathway',:].copy()

        df_sets_res['dataset_name'] = tissue
        dfs_metrics_I_lfc.append(df_sets_res)
        del df_sets_res

        # * ambiguity propagation score exclusively for families
        #if ablation_dir != '3_ablation_pathways_alone':
        try:
            df_aps= compute_ambig_propag_across_families(
                 tissue,
                subfolder=os.path.join(ablation_dir, tissue),
                data_centralized_folder=out_dir)
            df_aps['eval_type'] = ablation_type
            df_aps['set_type'] = "ghost"
            df_aps['dataset_name'] = tissue
        except Exception as e:
            print(e, "\n » impossible ghost scores, adding nan")
            df_aps = pd.DataFrame({'pathway_id': ['MolFamImpossible'],
                               'ambig_propag_score' : [np.nan],
                               'eval_type': [ablation_type],
                               'set_type': ['ghost'],
                               'dataset_name': [tissue] })

        # # add a fake pathway for consistent plot aesthetics
        df_aps = pd.concat([df_aps, pd.DataFrame({'pathway_id': ['WPZ'],
                               'ambig_propag_score' : [np.nan],
                              'eval_type': [ablation_type],
                               'set_type': ['pathway'],
                               'dataset_name': [tissue] })],
            ignore_index=True )
        dfs_metric_new.append(df_aps)


df_piled = pd.concat(dfs_metrics_I_lfc, ignore_index=True)
df_ambig_propag_score = pd.concat(dfs_metric_new, ignore_index=True)


df_piled = categorize_specific_columns(df_piled, tissues)
print( df_piled.columns)


df_ambig_propag_score = categorize_specific_columns(df_ambig_propag_score, tissues)

# ------ plot metric 1
save_plot_metric_tissue_specific(df_piled, metric_to_plot="logfoldchanges",
                                 out_file="metric1.pdf")

# ------ plot metric 2
df_I_piled = df_piled.copy()
df_I_piled = df_I_piled.drop(columns=['clus_set', 'logfoldchanges', 'pvals', 'pvals_adj' ])
df_I_piled = df_I_piled.drop_duplicates()
df_I_piled = df_I_piled.sort_values(by=[ 'eval_type' ,'set_type', 'pathway_id'], ascending=[True, True,True])

save_plot_metric_tissue_specific(df_I_piled, metric_to_plot="I",
                                 out_file="metric2.pdf")


# ------- plot metric 3
save_plot_metric_tissue_specific(df_ambig_propag_score,
                                 metric_to_plot="ambig_propag_score",
                                out_file="metric3.pdf",
                                 sharey_bool=False)

# ------- plot metric 4: number of significant sets (ghost and pathways)
df_piled['abslfc'] = np.abs(df_piled['logfoldchanges'].to_numpy())
# identical filter as  in SpacePath
signif_sets = df_piled.loc[(df_piled['pvals_adj'].astype(float) <= 0.05) & (
    df_piled['abslfc'].astype(float) >= 0.8) & (
    df_piled['I'] >= 0.2), 'pathway_id'].unique().tolist()

df_cc = df_I_piled.copy()
df_cc = df_cc.loc[df_cc['pathway_id'].isin(signif_sets), :].copy()
df_cc = df_cc.drop_duplicates()
df_cc['cc_var'] = df_cc['dataset_name'].astype(str) + '|' + df_cc['set_type'].astype(str) + '|' + df_cc['eval_type'].astype(str)
cool = df_cc.groupby(['cc_var']).count()
cool.rename(columns={'pathway_id':'n_sets'}, inplace=True)
cool.reset_index(inplace=True)
cool[['dataset_name', 'set_type', 'eval_type']] = cool['cc_var'].str.split("|", expand=True)

cool = categorize_specific_columns(cool, tissues)



def save_n_sets(df: pd.DataFrame, metric_to_plot='n_sets',
                  out_file=''):
    g = sns.catplot(
        data=df,
        x="eval_type",
        y=metric_to_plot,
        hue="set_type",
        dodge=True,
        gap=0.89,
        row="dataset_name",
        kind="bar",
        sharey=False,
        sharex=True,
        height=1.5,
        aspect=2.5,
        width=0.8,
        palette=palette_set_type_pastel
    )

    for tissue, ax in zip(g.row_names, g.axes.flat):
        sns.stripplot(
            data=df[df["dataset_name"] == tissue],
            x="eval_type",
            y=metric_to_plot,
            hue="set_type",
            dodge=True,
            size=7,
            palette=palette_set_type_dark,
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

    g.figure.set_size_inches(w=width_subplots, h=height_subplots)  # width=8", height=10"
    if out_file != '':
        g.figure.savefig(
            out_file,
            bbox_inches="tight"
        )
    else:
        print("!! missing out file for figure !!")


save_n_sets(cool, metric_to_plot='n_sets', out_file='numbers_sets_significant.pdf')


### end
