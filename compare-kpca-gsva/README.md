# Compare k-PCA vs GSVA
**Across spatially resolved metabolomics/lipidomics datasets**

This comparison was performed using the datasets specified in:

```
the_data_list.tsv
```
which is provided within this current location (`compare-kpca-gsva`).

The steps were performed in order, always in the current location `compare-kpca-gsva`:
* `cd scripts0` then `./copy_files.sh`: to copy the base .h5ad files and their annotations, plus the **KEGG** pathways
* `cd scripts1` then `module load smpath` (in our server. Adapt to yours) then `python3 run.py`: to obtain compound resolved matrix per dataset
* `cd scripts2` then `sbatch slurm_launch.sh`: to run **GSVA** (it includes `module load r_gsva`), can take up to 1 hour
* `cd scripts3` then  `module load smpath` (in our sever. Adapt to yours) and finally: 

  
```
python3 compare_outs__gsva_kpca.py
python3 plots_compare_outs.py
python3 viz_runtime.py
```

which ran **k-PCA** (via SpacePath which in turn uses the sspa implementation of k-PCA).

The input matrices (compound resolved) as well as the reference pathways were the same for both methods.

---

Virtual environments:
- the `smpath` environment in our server has the same dependencies as `spacepath` which is explained in SpacePath's README.
- the `r_gsva.yml` file to create the `r_gsva` environment is located in `compare-kpca-gsva/scripts2`




______________

author: Johanna Galvis
