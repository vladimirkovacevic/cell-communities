import argparse as ap
import logging
import os
import re
import sys

import scanpy as sc
from sliding_window import SlidingWindow

if __name__ == '__main__':

    sc.settings.verbosity = 3      
    sc.settings.set_figure_params(dpi=300, facecolor='white')
    parser = ap.ArgumentParser(description='A script that performs cell communities clusterization on single and multiple slices of ST data.')
    parser.add_argument('-f', '--file', help='File path to file that contain data to be analyzed/clustered', type=str, required=True)
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


    parser.add_argument('-w', '--win_size', help='Window size for anlyzing the cell community', type=int,required=False, default=150)
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

    if args.file.endswith('.h5ad'):
        adata = sc.read(args.file)
    # elif args.file.endswith('.gef'):
    #     data = st.io.read_gef(file_path=args.file, bin_type='cell_bins')
    #     adata = st.io.stereo_to_anndata(data)
    else:
        raise AttributeError(f"File '{args.file}' extension is not .h5ad") # or .gef

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
        algo = all_methods[method](adata, **vars(args))
        algo.run()

