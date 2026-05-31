import sys
import os
import pandas as pd
import anndata as ad

params_dir = "../../../config/"
orig_data = "../../../dataLOCAL/"
x00_dir = "../../../output/0_baseline/"


params_f = pd.read_csv(os.path.join(params_dir, 'baseline.tsv'),
                       sep="\t", index_col=None, header=0)

for tissue in params_f['dataset_name'].tolist():

    if not os.path.exists(os.path.join(x00_dir, tissue)):
        os.makedirs(os.path.join(x00_dir, tissue))

    os.system(f"cp {orig_data}{tissue}/* {x00_dir}{tissue}/")


    print("ok, populated:  ", x00_dir, tissue)