import sys
import os
import pandas as pd
import anndata as ad

params_dir = "../../config/"
orig_data = "../../../data/"            # "../../../dataLOCAL/"
input_dir =  "../../input/"
output_dir = "../../output/"

params_f = pd.read_csv(os.path.join(params_dir, 'baseline.tsv'),
                       sep="\t", index_col=None, header=0)

initialized_files_and_folders = True  # set to False if input and output are empty

if not initialized_files_and_folders:

    #  a) copy from original location to the input of this benchmark
    for tissue in params_f['dataset_name'].tolist():

        if not os.path.exists(os.path.join(input_dir, tissue)):
            os.makedirs(os.path.join(input_dir, tissue))

        os.system(f"cp {orig_data}{tissue}/{tissue + '.h5ad'} {input_dir}{tissue}/{tissue + '.h5ad'}")

        foo = params_f.loc[
            params_f['dataset_name'] == tissue, 'annotations_to_use'].values[0]

        for type_ann in foo.split():
            print(type_ann)
            os.system(
                f"cp {orig_data}{tissue}/{tissue + '-' + type_ann + '.csv'} \
                  {input_dir}{tissue}/{tissue + '-' + type_ann + '.csv'}")

        print("ok, populated:  ", input_dir, tissue)


    # b)clear objects from any computed values
    for tissue in params_f['dataset_name'].tolist():

        adata = ad.read_h5ad(os.path.join(input_dir, tissue, tissue + '.h5ad')
                             )
        print(adata.var.columns)
        print(adata.obs.columns)
        print(adata.uns.keys())

        adata.var = adata.var.drop(columns=['quality_msm_off_fdr', 'I', 'pval_sim_fdr_bh', 'var_sim',
           'quality_autocorr', 'final_quality', 'msm', 'rawOffSampleProb', 'fdr', 'formula'], errors='ignore')

        adata.obs = adata.obs.drop(
            columns=['sum_in_tissue_quants', 'sum_out_tissue_quants', 'ratio_in_out',
           'spot_quality'], errors='ignore')

        if 'moranI' in list(adata.uns.keys()):
            del adata.uns['moranI']

        if 'is_X_TIC_normalized' in list(adata.uns.keys()):
            del adata.uns['is_X_TIC_normalized']

        if 'spot_quality_colors' in list(adata.uns.keys()):
            del adata.uns['spot_quality_colors']

        print("==========================")
        print(adata.var.columns)
        print(adata.obs.columns)
        print(adata.uns.keys())
        # overwrite cleared object
        adata.write_h5ad(os.path.join(input_dir, tissue, tissue + '.h5ad')
                             )


    #  c) now distribute all these clear inputs to the output to "intialize"
    for tissue in params_f['dataset_name'].tolist():
        for out_subdir in ['0_baseline', '1_ablation_nofilter',
                           '2_ablation_random_id']:
            if not os.path.exists(os.path.join(output_dir, out_subdir, tissue)):
                os.makedirs(os.path.join(output_dir, out_subdir, tissue))
            # fi

            os.system(
                f"cp {input_dir}{tissue}/*.* {output_dir}{out_subdir}/{tissue}/")

    # d) add the metabolites mappings to each parent folder
    for out_subdir in ['0_baseline', '1_ablation_nofilter',
                           '2_ablation_random_id']:
        os.system(
            f"cp -r {orig_data}metabolites-mappings/ {output_dir}{out_subdir}/")

    print("Done")
else:
    print("nothing to do")