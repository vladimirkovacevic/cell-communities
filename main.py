import logging
import time

import argparse as ap
import scanpy as sc

from community_detection import CommunityDetection

if __name__ == '__main__':
    start_time = time.perf_counter()

    parser = ap.ArgumentParser(description='A script that performs cell communities clusterization on single and multiple slices of ST data.')
    parser.add_argument('-f', '--files', help='csv list of file paths to file that contain data to be analyzed/clustered', type=str, required=True)
    parser.add_argument('-t', '--tfile', help='File path to Anndata object with calculated cell mixtures for data windows, output of calc_feature_matrix', type=str, required=False, default=None)
    parser.add_argument('-a', '--annotation', help='Annotation label for cell types', type=str, required=True)
    parser.add_argument('-o', '--out_path', help='Absolute path to store outputs', type=str, default='results')
    parser.add_argument('-c', '--cluster_algo', help='Clustering algorithm', type=str, required=False, default='leiden', choices={'leiden', 'spectral', 'agglomerative'})
    parser.add_argument('-r', '--resolution', help='Resolution of leiden clustering algorithm. Ignored for spectral and agglomerative.', type=float, required=False, default=0.2)
    parser.add_argument('-s', '--spot_size', help='Size of the spot on plot', type=float, required=False, default=30)
    parser.add_argument('-v', '--verbose', help='Show logging messages. 0 - Show warrnings, >0 show info, <0 no output generated.', type=int, default=0)
    parser.add_argument('-p', '--plotting', help='Save plots flag. 0 - No plotting/saving, 1 - save clustering plot, 2 - save all plots (cell type images, statisctics and cell mixture plots)', type=int, required=False, default=3)
    parser.add_argument('--project_name', help='Project name that is used to name a directory containing all the slices used', type=str, required=False, default="Project")
    parser.add_argument('--skip_stats', help='Skip statistics calculation on cell community clustering result. A table of cell mixtures and comparative spatial plots of cell types and mixtures will not be created.', type=bool, required=False, default=False)
    parser.add_argument('--total_cell_norm', help='Total number of cells per window mixture after normalization', type=int, required=False, default=10000)
    parser.add_argument('--downsample_rate', help='Rate by which the binary image of cells is downsampled before calculating the entropy and scatteredness metrics', type=int, required=False, default=50)
    parser.add_argument('--dpi', help='DPI (dots per inch) used for plotting figures', type=int, required=False, default=100)
    parser.add_argument('--num_threads', help='Number of threads that will be used to speed up community calling', type=int, required=False, default=5)
    parser.add_argument('--entropy_thres', help='Threshold value for spatial cell type entropy for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('--scatter_thres', help='Threshold value for spatial cell type scatteredness for filtering out overdispersed cell types', type=float, required=False, default=1.0)
    parser.add_argument('-w', '--win_sizes', help='Comma separated list of window sizes for analyzing the cell community', type=str, required=False, default='NA')
    parser.add_argument('--sliding_steps', help='Comma separated list of sliding steps for sliding window', type=str, required=False, default='NA')
    parser.add_argument('--min_cluster_size', help='Minumum number of cell for cluster to be plotted in plot_stats()', type=int, required=False, default=200)
    parser.add_argument('--min_perc_to_show', help='Minumum percentage of cell type in cluster for cell type to be plotted in plot_stats()', type=int, required=False, default=4)
    parser.add_argument('--min_num_celltype', help='Minimum number of cell types that have more than --min_perc_celltype in a cluster, for a cluster to be shown in plot_celltype_table()', type=int, required=False, default=1)
    parser.add_argument('--min_perc_celltype', help='Minimum percentage of cells of a cell type which at least min_num_celltype cell types need to have to show a cluster in plot_celltype_table()', type=int, required=False, default=10)
    parser.add_argument('--min_cells_coeff', help='Multiple od standard deviations from mean values where the cutoff for m', type=float, required=False, default=1.5)
    parser.add_argument('--color_plot_system', help='Color system for display of cluster specific windows.', type=str, required=False, default='rgb', choices={'hsv', 'rgb'})
    parser.add_argument('--save_adata', help='Save adata file with resulting .obs column of cell community labels', type=bool, required=False, default=False)
    parser.add_argument('--min_count_per_type', help='Minimum number of cells per cell type needed to use the cell type for cell communities extraction (in percentages)', type=float, required=False, default=0.1)
    parser.add_argument('--n_clusters', help='Number of clusters for spectral and agglomerative clustering. Ignored for leiden', type=int, required=False, default=10)

    args = parser.parse_args()

    if args.verbose == 0:
        logging.basicConfig(level=logging.WARNING, force=True, format='%(message)s')
    elif args.verbose > 0:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    else:
        logging.basicConfig(level=logging.NOTSET, format='%(message)s')
    
    slices = []
    for file in args.files.split(','):
        # READ CELL TYPE ADATA
        if file.endswith('.h5ad'):
            adata = sc.read(file)
            slices.append(adata)
        else:
            # TODO: Consider adding GEF support
            raise AttributeError(f"File '{file}' extension is not .h5ad")  # or .gef
    cd = CommunityDetection(slices, **vars(args))
    cd.run()

    end_time = time.perf_counter()
    total_time = end_time - start_time

    logging.info(f'main.py took {total_time:.4f}s')
    logging.warning('END')
