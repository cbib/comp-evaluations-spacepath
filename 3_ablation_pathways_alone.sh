#!/bin/bash
##############################################################################
## Pathways alone: copy the sets results of baseline into the output folder  ##
##############################################################################

# Does not re-run SpacePath as it is already available in the baseline output

BL_DIR="output/0_baseline"
OUT_DIR="output/3_ablation_pathways_alone"

list_tissues=$(cat config/baseline.tsv | cut -f 1)

#rm -r $OUT_DIR/*/   # clean it

for tiss in ${list_tissues[@]}; do
    echo $tiss

    mkdir -p "$OUT_DIR/${tiss}"

    cp "${BL_DIR}/${tiss}/"*scores*.h5ad "${OUT_DIR}/$tiss/"
    cp "${BL_DIR}/${tiss}/"*differential-SETS.tsv "${OUT_DIR}/$tiss/"
    cp "${BL_DIR}/${tiss}/"*.mmt "${OUT_DIR}/${tiss}/"
    cp "${BL_DIR}/${tiss}/"*CleanAnnotations*.* "${OUT_DIR}/${tiss}/"
done

echo "OK"
