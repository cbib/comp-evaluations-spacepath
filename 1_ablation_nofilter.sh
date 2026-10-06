#!/bin/bash
##############################################################################
##          Runs SpacePath without filtering the m/z features               ##
##############################################################################
# works with 'spacepath' environment activated


PARAMS_FILE="config/no_filter.tsv"   # config for retaining all features
# same parent folder for both figures and results
OUT_DIR="output/1_ablation_nofilter/"
FIGS_DIR="output/1_ablation_nofilter/"

SpacePath $PARAMS_FILE $OUT_DIR $FIGS_DIR

