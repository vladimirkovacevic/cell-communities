#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate gcn

# 1. Whole brain (MOSTA) - 1 slice 
python main.py -f /goofys/Samples/Stereo_seq/E16.5_E1S3_cell_bin_whole_brain_noborderct.h5ad -o results/whole_brain -a "sim anno" --win_sizes 150 --sliding_steps 50 --scatter_thres 0.12 --downsample_rate 80 -c spectral --n_clusters 16 --resolution 0.25 --hide_plots --plotting 4 &

# 2. Adult brain - 3 slices
python main.py -f /goofys/Samples/Stereo_seq/Mouse_brain/SS200000141TL_B5.h5ad,/goofys/Samples/Stereo_seq/Mouse_brain/SS200000141TL_A4.h5ad,/goofys/Samples/Stereo_seq/Mouse_brain/SS200000128TR_E2.h5ad -o results/adult_brain -a celltype_pred --win_sizes 200 --sliding_steps 50 -c agglomerative --n_clusters 16 --resolution 0.2 --spot_size 20 --hide_plots --plotting 4 &

# 3. Mouse testis - 2 slices
# https://www.sciencedirect.com/science/article/pii/S2211124721013887
python main.py -f /goofys/Samples/slide_seq/mouse_testis/wt/WT3_ct.h5ad,/goofys/Samples/slide_seq/mouse_testis/diabetes/Diabetes2_ct.h5ad -o results/testis/ -a annotation --win_sizes 150 --sliding_steps 50 -c leiden --n_clusters 8 --resolution 0.25 --hide_plots --plotting 4 &

# 4. Mouse kidney - Figure 2
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8971939/
# https://cellxgene.cziscience.com/collections/8e880741-bf9a-4c8e-9227-934204631d2a
python main.py -f /goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_191204_15.h5ad,/goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_191204_22.h5ad -o results/kidney_fig2 -a author_cell_type --win_sizes 150 --sliding_steps 50 -c leiden --n_clusters 8 --resolution 0.3 --spot_size 12 --min_cluster_size 100 --hide_plots --plotting 4 &

# 5. Mouse kidney - Figure 3
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8971939/
# https://cellxgene.cziscience.com/collections/8e880741-bf9a-4c8e-9227-934204631d2a
python main.py -f /goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_191223_19.h5ad,/goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_200104_07.h5ad -o results/kidney_fig3 -a author_cell_type --win_sizes 300,150 --sliding_steps 150,50 -c leiden --n_clusters 7 --resolution 0.2 --spot_size 12 --hide_plots --plotting 4 &

# 6. Carcinoma - tumor, normal, border
# https://www.biorxiv.org/content/10.1101/2021.10.21.465135v1.full
python main.py -f /goofys/projects/ctc/LC6-M_DW1_web1_bin50/LC5-P_FE2.h5ad,/goofys/projects/ctc/LC6-M_DW1_web1_bin50/LC5-M_DU3.h5ad,/goofys/projects/ctc/LC6-M_DW1_web1_bin50/LC5-T_FD4.h5ad -o results/carcinoma -a cell_type --win_sizes 12 --sliding_steps 6 --downsample_rate 1 -c agglomerative --n_clusters 6 --resolution 0.1 --spot_size 4 --min_cells_coeff 1 --hide_plots --plotting 4 &

wait
