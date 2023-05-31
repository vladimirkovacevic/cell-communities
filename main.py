import os
import logging
import time

import argparse as ap
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

from core import *


if __name__ == '__main__':
    start_time = time.perf_counter()

    parser = ap.ArgumentParser(description='A script that performs cell communities clusterization on single and multiple slices of ST data.')
    parser.add_argument('-f', '--files', help='csv list of file paths to file that contain data to be analyzed/clustered', type=str, required=True)
    parser.add_argument('-t', '--tfile', help='File path to Anndata object with calculated cell mixtures for data windows, output of calc_feature_matrix', type=str, required=False, default=None)
    parser.add_argument('-a', '--annotation', help='Annotation label for cell types', type=str, required=True)
    parser.add_argument('-o', '--out_path', help='Absolute path to store outputs', type=str, default='results')
    parser.add_argument('-r', '--resolution', help='All: Resolution of the clustering algorithm', type=float, required=False, default=0.2)
    parser.add_argument('-s', '--spot_size', help='Size of the spot on plot', type=float, required=False, default=30)
    parser.add_argument('-v', '--verbose', help='Show logging messages. 0 - Show warrnings, >0 show info, <0 no output generated.', type=int, default=0)
    parser.add_argument('-p', '--plotting', help='Save plots flag. 0 - No plotting/saving, 1 - save clustering plot, 2 - save all plots (cell type images, statisctics and cell mixture plots)', type=int, required=False, default=2)
    parser.add_argument('--project_name', help='Project name that is used to name a directory containing all the slices used', type=str, required=False, default="Project")
    parser.add_argument('--skip_stats', help='Skip statistics calculation on cell community clustering result. A table of cell mixtures and comparative spatial plots of cell types and mixtures will not be created.', type=bool, required=False, default=False)
    parser.add_argument('--total_cell_norm', help='Total number of cells per window mixture after normalization', type=int, required=False, default=10000)
    parser.add_argument('--downsample_rate', help='Rate by which the binary image of cells is downsampled before calculating the entropy and scatteredness metrics', type=int, required=False, default=80)
    parser.add_argument('--num_threads', help='Number of threads that will be used to speed up community calling', type=int, required=False, default=5)
    parser.add_argument('--entropy_thres', help='Threshold value for spatial cell type entropy for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('--scatter_thres', help='Threshold value for spatial cell type scatteredness for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('-w', '--win_sizes', help='Comma separated list of window sizes for analyzing the cell community', type=str, required=False, default='150')
    parser.add_argument('--sliding_steps', help='Comma separated list of sliding steps for sliding window', type=str, required=False, default='50')
    parser.add_argument('--min_cluster_size', help='Minumum number of cell for cluster to be plotted in plot_stats()', type=int, required=False, default=500)
    parser.add_argument('--min_perc_to_show', help='Minumum percentage of cell type in cluster for cell type to be plotted in plot_stats()', type=int, required=False, default=5)
    parser.add_argument('--min_num_celltype', help='Minimum number of cell types that have more than --min_perc_celltype in a cluster, for a cluster to be shown in plot_celltype_table()', type=int, required=False, default=2)
    parser.add_argument('--min_perc_celltype', help='Minimum percentage of cells of a cell type which at least min_num_celltype cell types need to have to show a cluster in plot_celltype_table()', type=int, required=False, default=15)
    parser.add_argument('--min_cells_coeff', help='Multiple od standard deviations from mean values where the cutoff for m', type=float, required=False, default=1.5)
    parser.add_argument('--color_plot_system', help='Color system for display of cluster specific windows.', type=str, required=False, default='rgb', choices={'hsv', 'rgb'})
    parser.add_argument('--save_adata', help='Save adata file with resulting .obs column of cell community labels', type=bool, required=False, default=False)
    parser.add_argument('--min_count_per_type', help='Minimum number of cells per cell type needed to use the cell type for cell communities extraction (in percentages)', type=float, required=False, default=0.1)

    args = parser.parse_args()

    if args.verbose == 0:
        logging.basicConfig(level=logging.WARNING, force=True)
    elif args.verbose > 0:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.NOTSET)
        
    # Create directory to store outputs
    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)

    algo_list = []
    tissue_list = []
    win_sizes = "_".join([i for i in args.win_sizes.split(',')])
    args.project_name_orig = args.project_name
    args.out_path_orig = args.out_path
    args.project_name += f"_r{args.resolution}_ws{win_sizes}_en{args.entropy_thres}_sct{args.scatter_thres}_dwr{args.downsample_rate}_mcc{args.min_cells_coeff}"
    args.out_path = os.path.join(args.out_path, args.project_name)
    if not os.path.exists(args.out_path):
            os.mkdir(args.out_path)
    # FOR all slices
    for slice_id, file in enumerate(args.files.split(',')):
        # READ CELL TYPE ADATA
        if file.endswith('.h5ad'):
            adata = sc.read(file)
            adata.uns['slice_id'] = slice_id
            if 'X_spatial' in adata.obsm:
                adata.obsm['spatial'] = adata.obsm['X_spatial'].copy()
            elif 'spatial_stereoseq' in adata.obsm:
                adata.obsm['spatial'] = np.array(adata.obsm['spatial_stereoseq'].copy())
        else:
            # TODO: Consider adding GEF support
            raise AttributeError(f"File '{file}' extension is not .h5ad")  # or .gef
        # FEATURE EXTRACTION (SLIDING_WINDOW)
        algo = SlidingWindowMultipleSizes(adata, slice_id, file, **vars(args))
        # plot original annotation
        if args.plotting > 1:
            algo.plot_annotation()
        # run algorithm for feature extraction and cell type filtering based on entropy and scatteredness
        algo.run()
        if args.plotting > 1:
            algo.plot_histogram_cell_sum_window()
        # CELL TYPE FILTERING
        # [NOTE] This is not valid for multislice. A consensus on removing a cell type must
        # be made for all slices before removing it from any slice.
        # here I have tissue, I want to calculate entropy and scatteredness for each cell type in adata
        # and based on this information remove certain cell types
        entropy, scatteredness, cell_type_images = \
            calculate_spatial_metrics(algo.adata, algo.unique_cell_type, algo.downsample_rate, algo.annotation)
        # init var layer of tissue anndata object
        algo.tissue.var = algo.tissue.var.copy()
        algo.tissue.var.loc[:, 'entropy'] = entropy.loc[algo.tissue.var.index]
        algo.tissue.var.loc[:, 'scatteredness'] = scatteredness.loc[algo.tissue.var.index]
        algo.tissue.uns['cell_t_images'] = cell_type_images
        # save a .csv file with metrics per cell type
        algo.save_metrics()
        # plot binary images of cell types spatial positions
        if args.plotting > 3:
            algo.plot_celltype_images()
        # filter the cell types which are not localized using calculated metrics (entropy and scatteredness)
        algo.cell_type_filtering()

        # add algo object for each slice to a list
        algo_list.append(algo)
    
    if args.plotting > 0 and len(algo_list)>1:
        plot_all_annotation(args.out_path, algo_list)

    # MERGE TISSUE ANNDATA
    # each tissue has slice_id as 3rd coordinate in tissue.obsm['spatial']
    merged_tissue = ad.concat([a.get_tissue() for a in algo_list], axis=0, join='outer')
    # if tissues have different sets of cell those columns are filled with NANs
    # this is corrected by writing 0s
    merged_tissue.X[np.isnan(merged_tissue.X)] = 0.0

    # CLUSTERING (WINDOW_LABELS)
    sc.pp.neighbors(merged_tissue, use_rep='X')
    sc.tl.leiden(merged_tissue, resolution=args.resolution)

    for slice_id, algo in enumerate(algo_list):
        # extract clustering data from merged_tissue
        algo.set_clustering_labels(
            merged_tissue.obs.loc[merged_tissue.obsm['spatial'][:, 2] == slice_id, 'leiden'])

        # COMMUNITY CALLING (MAJORITY VOTING)
        algo.community_calling()

        # save anndata objects for further use
        if args.save_adata:
            algo.save_anndata()
        algo.save_community_labels()
        algo.save_tissue()

        # PLOT COMMUNITIES & STATISTICS
        # plot cell communities clustering result
        if args.plotting > 0:
            algo.plot_clustering()

        # if flag skip_stats is active, skip cell mixture statistics analysis
        if not args.skip_stats:
            algo.calculate_cell_mixture_stats()
            algo.save_mixture_stats()
            if args.plotting > 1:
                algo.plot_stats()
                algo.plot_celltype_table()
            if args.plotting > 2:
                algo.plot_cluster_mixtures()
                algo.boxplot_stats()
                algo.colorplot_stats(color_system=args.color_plot_system)
            if args.plotting > 3:
                algo.colorplot_stats_per_cell_types()
            # save final tissue with stats
            algo.save_tissue(suffix='_stats')
    
    if args.plotting > 0 and len(algo_list)>1:
        plot_all_clustering(args.out_path, algo_list)
    if args.plotting > 2:
        plot_celltype_mixtures_total([algo.get_cell_mixtures().to_dict() for algo in algo_list], args.out_path)
        plot_cell_abundance_total(algo_list, args.out_path)
        plot_cluster_abundance_total(algo_list, args.out_path)
    if args.plotting > 3:
        plot_cell_perc_in_community_per_slice(algo_list, args.out_path)
        plot_cell_abundance_per_slice(algo_list, args.out_path)
        plot_cluster_abundance_per_slice(algo_list, args.out_path)

    generate_report(args)
    end_time = time.perf_counter()
    total_time = end_time - start_time
    print(f'main.py took {total_time:.4f}s')

    logging.info(f'main.py took {total_time:.4f}s')
    logging.warning('END')
