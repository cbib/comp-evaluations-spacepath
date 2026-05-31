#!/bin/bash
##############################################################################
##  Runs SpacePath, selecting random compound ID for pathways and families  ##
##############################################################################

# works with smpath env activated

CONFIG="config/baseline.tsv"  # ok as the changement is in the .py script

OUT_DIR="output/2_ablation_random_id"
FIGS_DIR=$OUT_DIR

python src/ablation2_random_id.py $CONFIG $OUT_DIR $FIGS_DIR
