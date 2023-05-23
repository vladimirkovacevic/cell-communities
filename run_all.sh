#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate gcn

# 1. Whole brain (MOSTA) - 1 slice 
python main.py -f /goofys/Samples/Stereo_seq/E16.5_E1S3_cell_bin_whole_brain_spagft_4_noborderct.h5ad -o /home/ubuntu/results/whole_brain/ --win_sizes 150 --sliding_steps 50 -a "sim anno" --scatter_thres 0.12 --downsample_rate 80 --spot_size 30 --resolution 0.25 --plotting 3 --min_cells_coeff 1.5 --min_num_celltype 1 --min_perc_celltype 10 --min_perc_to_show 8 --color_plot_system rgb &

# 2. Adult brain - 3 slices
python main.py -f /goofys/Samples/Stereo_seq/Mouse_brain/SS200000128TR_E2.h5ad,/goofys/Samples/Stereo_seq/Mouse_brain/SS200000141TL_A4.h5ad,/goofys/Samples/Stereo_seq/Mouse_brain/SS200000141TL_B5.h5ad -o /home/ubuntu/results/adult_brain/ -a 'celltype_pred' --win_sizes 300 --sliding_steps 50 --resolution 0.2 --spot_size 50 --plotting 3 &

# 3. Mouse testis - 2 slices
# https://www.sciencedirect.com/science/article/pii/S2211124721013887
python main.py -f /goofys/Samples/slide_seq/mouse_testis/wt/WT1_ct.h5ad,/goofys/Samples/slide_seq/mouse_testis/diabetes/Diabetes1_ct.h5ad -o /home/ubuntu/results/testis/ -a annotation --win_sizes 150 --sliding_steps 50 --resolution 0.2 --downsample_rate 40 --entropy_thres 1.0 --scatter_thres 1.0 --min_cells_coeff 1.5 --plotting 3 &

# 4. Mouse kidney - Figure 2
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8971939/
# https://cellxgene.cziscience.com/collections/8e880741-bf9a-4c8e-9227-934204631d2a
python main.py -f /goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_191204_15.h5ad,/goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_191204_22.h5ad -o /home/ubuntu/results/kidney -a author_cell_type --win_sizes 150 --sliding_steps 50 --resolution 0.3 --spot_size 20 --plotting 3 &

# 5. Mouse kidney - Figure 3
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8971939/
# https://cellxgene.cziscience.com/collections/8e880741-bf9a-4c8e-9227-934204631d2a
python main.py -f /goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_191223_19.h5ad,/goofys/Samples/slide_seq/cellxgene_kidney_slide_seq_v2/Puck_200104_07.h5ad -o /home/ubuntu/results/kidney -a author_cell_type --win_sizes 300,150 --sliding_steps 150,50 --resolution 0.2 --spot_size 20 --plotting 3 &

# 6. Carcinoma - tumor, normal, border
# https://www.biorxiv.org/content/10.1101/2021.10.21.465135v1.full
python main.py -f /goofys/projects/ctc/LC6-M_DW1_web1_bin50/LC5-M_DU3_web3_bin35_CelltypeTrans_MaligComb_Spotlight_First.h5ad,/goofys/projects/ctc/LC6-M_DW1_web1_bin50/LC5-P_FE2_web0_bin50_CelltypeTrans_MaligComb_Spotlight_First.h5ad,/goofys/projects/ctc/LC6-M_DW1_web1_bin50/LC5-T_FD4_web2_bin50_CelltypeTrans_MaligComb_Spotlight_First.h5ad -o /home/ubuntu/results/LC5 -a cell_type --win_sizes 9 --sliding_steps 3 --resolution 0.05 --spot_size 1.5 --plotting 3 &


wait
