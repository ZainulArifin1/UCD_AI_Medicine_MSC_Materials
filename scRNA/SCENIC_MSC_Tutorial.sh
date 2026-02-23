#!/bin/bash
conda activate scenic_docker
        
# PAN ALL
docker run --rm \
    -v /mnt/8TB/Projects/POIAZ/Zainul/scRNA_PD1:/data \
    aertslab/pyscenic:0.12.1 pyscenic grn \
        --num_workers 5 \
        --transpose \
        -o /data/SCENIC_out/MSC/MSC_adj.csv \
        /data/Data_sc/MSC/Input_Scenic.tsv \
        /data/cistarget/TF_hg38.txt

docker run --rm \
    -v /mnt/8TB/Projects/POIAZ/Zainul/scRNA_PD1:/data \
    aertslab/pyscenic:0.12.1 pyscenic ctx \
        /data/SCENIC_out/MSC/MSC_adj.csv \
        /data/cistarget/hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
        /data/cistarget/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
        --annotations_fname /data/cistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname /data/Data_sc/MSC/Input_Scenic.tsv \
        --mode "custom_multiprocessing" \
        --output /data/SCENIC_out/MSC/MSC_regulon.csv \
        --transpose \
        --num_workers 5

docker run --rm \
    -v /mnt/8TB/Projects/POIAZ/Zainul/scRNA_PD1:/data \
    aertslab/pyscenic:0.12.1 pyscenic aucell \
        /data/Data_sc/MSC/Input_Scenic.tsv \
        /data/SCENIC_out/MSC/MSC_regulon.csv  \
        -o /data/SCENIC_out/MSC/MSC_auc_mtx.csv \
        --transpose \
        --num_workers 5