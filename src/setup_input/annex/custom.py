import sys
import os
import pandas as pd

params_dir = "../../../config/"
input_dir = "../../../input/"

params_f = pd.read_csv(os.path.join(params_dir, 'baseline.tsv'),
                       sep="\t", index_col=None, header=0)


f1 = pd.read_csv(os.path.join(
    params_dir, "the_data_list.tsv"
), sep="\t", index_col=None, header=0)

f2 = pd.read_csv(os.path.join(
    params_dir, "public-datasets-parameters.tsv"
), sep="\t", index_col=None, header=0)

f3 = pd.read_csv(os.path.join(
    params_dir, "private-datasets-parameters.tsv"
), sep="\t", index_col=None, header=0)

f4 = pd.concat([f2, f3])

some_already_omit = list(set(
    f4["dataset_name"].tolist()) - set(f1['dataset_name'].tolist()))

more_omit = ['UBcalibr2', 'MALDI_275_T', 'MALDI_248_T', 'Heart-sgt-mouse', 'Kidney-xenogr-1',
            'Brain-sgt-mouse', 'Brain-wt1-mouse']

all_omit = set(some_already_omit).union(
    set(more_omit)
)
omit_df = f4.loc[f4["dataset_name"].isin(all_omit), :]
include_df = f4.loc[~f4["dataset_name"].isin(all_omit), :]

include_df.to_csv(os.path.join(params_dir, 'baseline.tsv'),
                            sep='\t', header=True, index=False)
omit_df.to_csv(os.path.join(params_dir, 'for_later.tsv'),
                            sep='\t', header=True, index=False)

