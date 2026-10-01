#!/bin/bash
##############################################################################
##               Runs SpacePath with parameters as published                ##
##############################################################################
# with the recommended parameters as stated in the config/baseline.tsv  file
# works with 'spacepath' environment activated
echo "First make sure that src/setup_input/input.py populated ../input folder"
echo ""
echo "Calling SpacePath"

CONFIG="config/baseline.tsv"   # config/minimal_foo.tsv

# same parent folder for both figures and results
OUT_DIR="output/0_baseline/"
FIGS_DIR="output/0_baseline/"

SpacePath $CONFIG $OUT_DIR $FIGS_DIR


