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
    parser.add_argument('--total_cell_norm', help='Total number of cells per window mixture after normalization', type=int,required=False, default=10000)
    parser.add_argument('--downsample_rate', help='Rate by which the binary image of cells is downsampled before calculating the entropy and scatteredness metrics', type=int,required=False, default=80)
    parser.add_argument('--entropy_thres', help='Threshold value for spatial cell type entropy for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('--scatter_thres', help='Threshold value for spatial cell type scatteredness for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('-w', '--win_sizes', help='Comma separated list of windows sizes', type=str,required=False, default='150')
    parser.add_argument('--sliding_step', help='Slide step for sliding window method', type=int, required=False, default=50)

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
        all_methods['sliding_window'] = SlidingWindow

    win_sizes_split = args.win_sizes.split(',')
    win_size_list = [int(w) for w in win_sizes_split]

    # Process requested methods
    for method in all_methods:
        algo_dict = {}
        tissue_list = []
        
        for slice_id, file in enumerate(args.files.split(',')):
            algo_list = []

            # READ CELL TYPE ADATA
            if file.endswith('.h5ad'):
                adata = sc.read(file)
                adata.uns['slice_id'] = slice_id
            else:
                # TODO: Consider adding GEF support
                raise AttributeError(f"File '{file}' extension is not .h5ad") # or .gef

            # FOR
            for win_size in win_size_list:
                args.win_size = win_size
                # FEATURE EXTRACTION (SLIDING_WINDOW)
                algo = all_methods[method](adata, slice_id, file, **vars(args))
                # run algorithm for feature extraction and cell type filtering based on entropy and scatteredness
                algo.run()

                # CELL TYPE FILTERING
                # [NOTE] This is not valid for multislice. A consensus on removing a cell type must
                # be made for all slices before removing it from any slice.
                # here I have tissue, I want to calculate entropy and scatteredness for each cell type in adata
                # and based on this information remove certain cell types
                algo.tissue.var['entropy'], algo.tissue.var['scatteredness'], algo.tissue.uns['cell_t_images'] = \
                    calculate_spatial_metrics(algo.adata, algo.unique_cell_type, algo.downsample_rate, algo.annotation)
                # save a .csv file with metrics per cell type
                algo.save_metrics()
                # plot binary images of cell types spatial positions
                if args.plotting > 1: algo.plot_celltype_images()
                # filter the cell types which are not localized using calculated metrics (entropy and scatteredness)
                algo.cell_type_filtering()

                # add algo object for each slice and each window size to a list
                algo_list.append(algo)
            
            algo_dict[slice_id] = algo_list
        
        algos = []
        for algo_list in algo_dict.values():
            algos = algos + algo_list
        # MERGE TISSUE ANNDATA
        # each tissue has slice_id as 3rd coordinate and windows_size as 4rd coordinate in tissue.obsm['spatial']
        merged_tissue = ad.concat([a.get_tissue() for a in algos], axis=0, join='outer')

        # CLUSTERING (WINDOW_LABELS)
        sc.pp.neighbors(merged_tissue, use_rep='X')
        sc.tl.leiden(merged_tissue, resolution=args.resolution)




        # COMMUNITY CALLING NE ZA SUBWINDOW NEGO ZA CELIJU

        for slice_id, algo_list in algo_dict.items():
            # extract clustering data from merged_tissue

            # community calling (majority voting) se pokrece za iste za iste slice-ove
            # ovaj korak ostaje
            algo_list[slice_id].set_clustering_labels(
                merged_tissue.obs.loc[merged_tissue.obsm['spatial'][:, 2] == slice_id, 'leiden'])
            
            # COMMUNITY CALLING (MAJORITY VOTING)
            algo_list[slice_id].community_calling()

            # save anndata object for further use
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
