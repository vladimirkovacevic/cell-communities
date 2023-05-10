import os
import logging

import argparse as ap
import anndata as ad
import scanpy as sc

from core import *


if __name__ == '__main__':

    sc.settings.verbosity = 3      
    sc.settings.set_figure_params(dpi=300, facecolor='white')
    parser = ap.ArgumentParser(description='A script that performs cell communities clusterization on single and multiple slices of ST data.')
    parser.add_argument('-f', '--files', help='csv list of file paths to file that contain data to be analyzed/clustered', type=str, required=True)
    parser.add_argument('-t', '--tfile', help='File path to Anndata object with calculated cell mixtures for data windows, output of calc_feature_matrix', type=str, required=False, default=None)
    parser.add_argument('-a', '--annotation', help='Annotation label for cell types', type=str, required=True)
    parser.add_argument('-m', '--methods', help='Comma separated list of methods to perform. Available: sliding_window', type=str, default='sliding_window')
    parser.add_argument('-o', '--out_path', help='Absolute path to store outputs', type=str, default='results')
    parser.add_argument('-r', '--resolution', help='All: Resolution of the clustering algorithm', type=float, required=False, default=0.2)
    parser.add_argument('-s', '--spot_size', help='Size of the spot on plot', type=float, required=False, default=30)
    parser.add_argument('-v', '--verbose', help='Show logging messages. 0 - Show warrnings, >0 show info, <0 no output generated.', type=int, default=0)
    parser.add_argument('-p', '--plotting', help='Save plots flag. 0 - No plotting/saving, 1 - save clustering plot, 2 - save all plots (cell type images, statisctics and cell mixture plots)', type=int, required=False, default=2)
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
    parser.add_argument('--min_cells_coeff', help='Multiple od standard deviations from mean values where the cutoff for m', type=float, required=False, default=1.5)
    parser.add_argument('--save_adata', help='Save adata file with resulting .obs column of cell community labels', type=bool, required=False, default=False)

    args = parser.parse_args()

    if args.verbose == 0:
        logging.basicConfig(level=logging.WARNING, force=True)
    elif args.verbose > 0:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.NOTSET)

    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)

    # Parse requested and installed methods to make sure that requested methods are installed
    available_methods = ['sliding_window']

    chosen_methods = args.methods.split(',')
    assert set(chosen_methods).issubset(set(available_methods)), \
        "The requested methods could not be executed because your environment lacks needed libraries."
    all_methods = {}  # Only one method for now
    if 'sliding_window' in chosen_methods:
        all_methods['sliding_window'] = SlidingWindowMultipleSizes

    # Process requested methods
    for method in all_methods:
        algo_list = []
        tissue_list = []
        # FOR
        for slice_id, file in enumerate(args.files.split(',')):
            # READ CELL TYPE ADATA
            if file.endswith('.h5ad'):
                adata = sc.read(file)
                adata.uns['slice_id'] = slice_id
            else:
                # TODO: Consider adding GEF support
                raise AttributeError(f"File '{file}' extension is not .h5ad") # or .gef
            # FEATURE EXTRACTION (SLIDING_WINDOW)
            algo = all_methods[method](adata, slice_id, file, **vars(args))
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
            if args.plotting > 1:
                algo.plot_celltype_images()
            # filter the cell types which are not localized using calculated metrics (entropy and scatteredness)
            algo.cell_type_filtering()
            
            # add algo object for each slice to a list
            algo_list.append(algo)
        
        # MERGE TISSUE ANNDATA
        # each tissue has slice_id as 3rd coordinate in tissue.obsm['spatial']
        merged_tissue = ad.concat([a.get_tissue() for a in algo_list], axis=0, join='outer')

        # CLUSTERING (WINDOW_LABELS)
        sc.pp.neighbors(merged_tissue, use_rep='X')
        sc.tl.leiden(merged_tissue, resolution=args.resolution)

        for slice_id, _ in enumerate(args.files.split(',')):
            # extract clustering data from merged_tissue
            algo_list[slice_id].set_clustering_labels(
                merged_tissue.obs.loc[merged_tissue.obsm['spatial'][:, 2] == slice_id, 'leiden'])
           
            # COMMUNITY CALLING (MAJORITY VOTING)
            algo_list[slice_id].community_calling()

            # save anndata objects for further use
            if args.save_adata:
                algo_list[slice_id].save_adata()
            algo_list[slice_id].save_tissue()

            # PLOT COMMUNITIES & STATISTICS
            # plot cell communities clustering result
            if args.plotting > 0:
                algo_list[slice_id].plot_clustering()

            # if flag skip_stats is active, skip cell mixture statistics analysis
            if not args.skip_stats:
                algo_list[slice_id].calculate_cell_mixture_stats()
                algo_list[slice_id].save_mixture_stats()
                if args.plotting > 1:
                    algo_list[slice_id].plot_stats()
                # save final tissue with stats
                algo_list[slice_id].save_tissue(suffix='_stats')

    logging.warning('END')
