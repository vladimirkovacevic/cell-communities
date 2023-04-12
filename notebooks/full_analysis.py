import argparse as ap
import anndata as ad
import logging
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import scanpy as sc
from sliding_window import SlidingWindow

def community_calling(adata, tissue, win_size, sliding_step, method_key):
    bin_slide_ratio = int(win_size/sliding_step)
    x_min = adata.obs['Centroid_X'].min()
    y_min = adata.obs['Centroid_Y'].min()
    # max voting on cluster labels
    # init the new obs column
    tissue.obs['leiden_max_vote'] = list('x' for x in range(len(tissue.obs.index)))
    for x_curr, y_curr, _ in tissue.obsm['spatial']:
        # index of subwindow is in the top left corner of the whole window
        subwindow_labels = {}
        for slide_x in range(0, np.min([bin_slide_ratio, x_curr - x_min + 1])):
            for slide_y in range(0, np.min([bin_slide_ratio, y_curr - y_min + 1])):
                # check if location exist (spatial area is not complete)
                if (f'{x_curr - slide_x}_{y_curr - slide_y}') in tissue.obs.index:
                    new_value = tissue.obs.loc[f'{x_curr - slide_x}_{y_curr - slide_y}', 'leiden']
                    subwindow_labels[new_value] = subwindow_labels[new_value] + 1 if new_value in subwindow_labels.keys() else 1
        
        # max vote
        # max vote should be saved in a new obs column so that it does not have diagonal effect on
        # other labels during refinment
        tissue.obs.loc[f'{x_curr}_{y_curr}', 'leiden_max_vote'] = max(subwindow_labels, key=subwindow_labels.get)

    adata.obs[f'tissue_{method_key}'] = list(tissue.obs.loc[adata.obs['x_y'], 'leiden_max_vote'])

    logging.info(r"Sliding window cell mixture calculation done. Added results to adata.obs['sliding_window']")


if __name__ == '__main__':

    sc.settings.verbosity = 3      
    sc.settings.set_figure_params(dpi=300, facecolor='white')
    parser = ap.ArgumentParser(description='A script that performs cell communities clusterization on single and multiple slices of ST data.')
    parser.add_argument('-f', '--files', help='csv list of file paths to file that contain data to be analyzed/clustered', type=str, required=True)
    parser.add_argument('-t', '--tfile', help='File path to Anndata object with calculated cell mixtures for data windows, output of calc_feature_matrix', type=str, required=False, default=None)
    parser.add_argument('-a', '--annotation', help='Annotation label for cell types', type=str, required=True)
    parser.add_argument('-m', '--methods', help='Comma separated list of methods to perform. Available: sliding_window', type=str, required=True, default='sliding_window')
    parser.add_argument('-o', '--out_path', help='Absolute path to store outputs', type=str, required=True)
    parser.add_argument('-r', '--resolution', help='All: Resolution of the clustering algorithm', type=float, required=False, default=0.2)
    parser.add_argument('-s', '--spot_size', help='Size of the spot on plot', type=float, required=False, default=30)
    parser.add_argument('-v', '--verbose', help='Show logging messages. 0 - Show warrnings, >0 show info, <0 no output generated.', type=int, default=0)
    
    parser.add_argument('--total_cell_norm', help='Total number of cells per window mixture after normalization', type=int,required=False, default=10000)
    parser.add_argument('--downsample_rate', help='Rate by which the binary image of cells is downsampled before calculating the entropy and scatteredness metrics', type=int,required=False, default=80)
    parser.add_argument('--entropy_thres', help='Threshold value for spatial cell type entropy for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('--scatter_thres', help='Threshold value for spatial cell type scatteredness for filtering out overdispersed cell types', type=float, required=False, default=1.0)

    parser.add_argument('-w', '--win_size', help='Window size for anlyzing the cell community', type=int, required=False, default=150)
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

    # # Most algorithms demand sparse cell gene matrix
    # if not scipy.sparse.issparse(adata.X):
    #     adata.X = scipy.sparse.csr_matrix(adata.X)

    # Parse requested and installed methods to make sure that requested methods are installed
    # available_methods = [module.__name__ for module in sys.modules.values() if re.search('^core.+', module.__name__)]
    # available_methods = [m.split('.')[1] for m in available_methods]
    available_methods = ['sliding_window']

    chosen_methods = args.methods.split(',')
    assert set(chosen_methods).issubset(set(available_methods)), "The requested methods could not be executed because your environment lacks needed libraries."
        
    all_methods = {}
    if 'sliding_window' in chosen_methods:
        all_methods['sliding_window'] = SlidingWindow

    
    # Process requested methods
    for method in all_methods:
        adata_list = []
        tissue_list = []
        # FOR
        for slice_id, file in enumerate(args.files.split(',')):
            # READ CELL TYPE ADATA
            if file.endswith('.h5ad'):
                adata = sc.read(file)
                adata.uns['slice_id'] = slice_id
                adata_list.append(adata)
            # elif args.file.endswith('.gef'):
            #     data = st.io.read_gef(file_path=args.file, bin_type='cell_bins')
            #     adata = st.io.stereo_to_anndata(data)
            else:
                raise AttributeError(f"File '{file}' extension is not .h5ad") # or .gef
            # FEATURE EXTRACTION (SLIDING_WINDOW) & CELL TYPE FILTERING
            algo = all_methods[method](adata, file, **vars(args))
            # run algorithm for feature extraction and cell type filtering based on entropy and scatteredness
            algo.run()
            # algo.plot_clustering(color=[f'tissue_{algo.method_key}'], sample_name=f'clusters_cellspots_{algo.params_suffix}.png')
            # algo.calculate_cell_mixture_stats()
            # algo.plot_stats()
            # algo.save_results()

            tissue_list.append(algo.get_tissue())
        
        # MERGE TISSUE ANNDATA
        # each tissue has slice_id as 3rd coordinate in tissue.obsm['spatial']
        merged_tissue = ad.concat(tissue_list, axis=0, join='outer')
        # CLUSTERING (WINDOW_LABELS)
        sc.pp.neighbors(merged_tissue, use_rep='X')
        sc.tl.leiden(merged_tissue, resolution=args.resolution)

        for slice_id, _ in enumerate(args.files.split(',')):
            # extract clustering data from merged_tissue
            tissue_list[slice_id].obs = tissue_list[slice_id].obs.copy()
            tissue_list[slice_id].obs['leiden'] = pd.Series(merged_tissue.obs.loc[merged_tissue.obsm['spatial'][:,2]==slice_id, 'leiden'], index=tissue_list[slice_id].obs.index)
            # COMMUNITY CALLING (MAJORITY VOTING)
            community_calling(adata=adata_list[slice_id], tissue=tissue_list[slice_id], win_size=args.win_size, sliding_step=args.sliding_step, method_key=args.methods)
            
            # PLOT COMMUNITIES & STATISTICS
            figure, ax = plt.subplots(nrows=1, ncols=1)
            sc.pl.spatial(adata_list[slice_id], color=f'tissue_{args.methods}', palette=None, spot_size=args.spot_size, ax=ax, show=False)
            ax.axis('off')
            figure.savefig(os.path.join(algo.dir_path, f'clusters_cellspots_slice{slice_id}.png'), dpi=300, bbox_inches='tight')
            plt.close()

        print('END')
        

