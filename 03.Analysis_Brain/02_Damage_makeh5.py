import scanpy as sc 
import pandas as pd
import numpy as np
import snapatac2 as snap
import anndata
import os
import re

sample_list = pd.read_table("/gpfs/commons/groups/zhu_lab/czhu/01.project_data/03.DNA_damage/11.Brain_dataset/sample_list.xls", 
                            sep="\t", header=0)

#read the fragment files

for i in range(len(sample_list)):
    RNA_lib = sample_list.iloc[i]['RNA']
    DNA_lib = sample_list.iloc[i]['DNA']

    # Define file paths
    frag_path = "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/fragments_mod/Modified_SnapATAC2_{}_fragment.tsv.gz".format(DNA_lib)
    out_file = "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/data/MouseBrain/snapATAC2_h5ad/MouseBrain_damage_{}.h5ad".format(DNA_lib)
    tsse_figs = "/gpfs/commons/groups/zhu_lab/zcao/Paired_damage/analysis/MouseBrain/DNA/Individual/{}_tsse.png".format(DNA_lib)
    # Import data
    adata = snap.pp.import_data(frag_path, 
                                chrom_sizes=snap.genome.mm10, 
                                file=out_file, 
                                n_jobs = -1)
    # Compute and plot TSSe metrics
    snap.metrics.tsse(adata, snap.genome.mm10)
    snap.pl.tsse(adata, interactive=False, out_file = tsse_figs)

    # Add tile matrix
    snap.pp.add_tile_matrix(adata, bin_size=50000, n_jobs=-1)
    adata.close()
    print("Finished:", DNA_lib)