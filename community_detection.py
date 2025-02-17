import os
import time
import logging
from collections import defaultdict
from functools import reduce
from itertools import cycle
from typing import List

import anndata as ad
import matplotlib.ticker as mtick
import seaborn as sns
import scanpy as sc
from matplotlib import pyplot as plt
from sklearn.cluster import SpectralClustering, AgglomerativeClustering

from anndata import AnnData
from ccd import *
from ccd.community_clustering_algorithm import cluster_palette

class CommunityDetection():
    """
    Class for performing community detection on a set of slices.
    """

    def __init__(
            self,
            slices: List[AnnData],
            annotation: str,
            **kwargs) -> None:
        """
        Initialize the CommunityDetection object.

        Example:
        import scanpy as sc
        from community_detection import CommunityDetection
        files = ... # List of input .h5ad file paths
        slices = []
        for file in files:
            adata = sc.read(file)
            slices.append(adata)
        cd = CommunityDetection(slices, **vars(args))
        cd.run()

        Parameters:
        - slices (List[AnnData]): A list of AnnData objects representing the slices of a tissue.
        - annotation (str): The annotation string.
        - **kwargs: Additional keyword arguments.
        """
        self.params = { **COMMUNITY_DETECTION_DEFAULTS, **kwargs }
        init_logger(level=self.params['verbose'])
        self.params['annotation'] = annotation
        self.slices = slices
        self.cell_types = set([])
        self.annotation_palette = None
        missing_cell_type_palette = False
        for slice in self.slices:
            if 'X_spatial' in slice.obsm:
                slice.obsm['spatial'] = slice.obsm['X_spatial'].copy()
            elif 'spatial_stereoseq' in slice.obsm:
                slice.obsm['spatial'] = np.array(slice.obsm['spatial_stereoseq'].copy())
            # annotation data must be of string type
            slice.obs[annotation] = slice.obs[annotation].astype('str')
            # create a set of existing cell types in all slices
            self.cell_types = self.cell_types.union(set(slice.obs[annotation].unique()))
            # if any of the samples lacks the cell type palette, set the flag
            if f'{annotation}_colors' not in slice.uns_keys():
                missing_cell_type_palette = True

        self.cell_types = list(sorted(self.cell_types))
                
        self.file_names = [fname for fname in self.params['files'].split(',')] if 'files' in self.params else [f"Slice_{id}" for id in range(len(slices))]
        if self.params['win_sizes'] == 'NA' or self.params['sliding_steps'] == 'NA':
            logging.info("Window sizes and/or sliding steps not provided by user - proceeding to calculate optimal values")
            self.params['win_sizes'], self.params['sliding_steps'] = self.calc_optimal_win_size_and_slide_step()
        else:
            self.log_win_size_full_info()

        # if downsample_rate is not provided, use the dimension of the smallest window size as reference
        if self.params['downsample_rate'] == None:
            logging.info("Downsample rate is not provided by user - proceeding to calculate one based on minimal window size.")
            min_win_size = np.min([int(i) for i in self.params['win_sizes'].split(',')])
            self.params['downsample_rate'] = min_win_size // 2
            logging.info(f"donwsample_rate = {self.params['downsample_rate']}")

        # if cell type palette is not available, create and add to each slice anndata object
        if missing_cell_type_palette:
            self.annotation_palette = {cellt: cluster_palette[-id-1] for id, cellt in enumerate(self.cell_types)}
            for slice in self.slices:
                slice.uns[f'{annotation}_colors'] = [self.annotation_palette[cellt] for cellt in list(sorted(slice.obs[annotation].unique()))]
        

    def run(self):
        """
        Executes the community detection algorithm.

        This method performs community detection using the specified parameters and data slices. It follows the following steps:

        1. Data Processing Loop:
        - Iterates over each data slice, identified by `slice_id` and the corresponding file name.
        - Initializes a `SlidingWindowMultipleSizes` algorithm object, `algo`, with the slice and other parameters.
        - Optionally plots the original annotation if the plotting level is greater than 0.
        - Runs the algorithm for feature extraction and cell type filtering based on entropy and scatteredness.
        - Optionally plots the histogram of cell sum per window if the plotting level is greater than 1.
        - Calculates entropy, scatteredness, and cell type images using the `calculate_spatial_metrics` function.
        - Optionally plots binary images of cell types' spatial positions if the plotting level is greater than 3.
        - Filters out cell types that are not localized based on calculated metrics.
        - Appends the `algo` object to the `algo_list`.

        2. Annotation Plotting:
        - If the plotting level is greater than 0 and there are multiple algorithm objects in `algo_list`, plots the annotation for all slices together.

        3. Tissue Merging:
        - Merges the tissue Anndata objects from all algorithm objects into a single Anndata object, `merged_tissue`.

        4. Clustering:
        - Performs clustering on the merged tissue using the specified clustering algorithm.

        5. Algorithm Execution:
        - Performs community calling using the majority voting method on each algorithm object.
        - Saves the Anndata objects, community labels, and tissue data for further use.
        - Optionally plots the clustering results if the plotting level is greater than 0.
        - If the `skip_stats` flag is not active, calculates cell mixture statistics, saves them, and generates corresponding plots.
        - Optionally saves the final tissue with statistics.
    
        6. Clustering Plotting:
        - If the plotting level is greater than 0 and there are multiple algorithm objects in `algo_list`, plots the clustering results for all slices together.

        7. Additional Plots:
        - Generates additional plots based on the plotting level and the data from `algo_list`. These include cell type mixtures, cell abundance, and cluster abundance plots.

        8. Report Generation:
        - Generates a HTML report o the results.
        """

        start_time = time.perf_counter()

        if not os.path.exists(self.params['out_path']):
            os.makedirs(self.params['out_path'])

        self.algo_list = []
        win_sizes = "_".join([i for i in self.params['win_sizes'].split(',')])
        sliding_steps = "_".join([i for i in self.params['sliding_steps'].split(',')])
        self.params['project_name_orig'] = self.params['project_name']
        self.params['out_path_orig'] = self.params['out_path']
        cluster_string = f"_r{self.params['resolution']}" if self.params['cluster_algo'] == 'leiden' else f"_nc{self.params['n_clusters']}"
        self.params['project_name'] += f"_c{self.params['cluster_algo']}{cluster_string}_ws{win_sizes}_ss{sliding_steps}_sct{self.params['scatter_thres']}_dwr{self.params['downsample_rate']}_mcc{self.params['min_cells_coeff']}"
        self.params['out_path'] = os.path.join(self.params['out_path'], self.params['project_name'])

        if not os.path.exists(self.params['out_path']):
            os.makedirs(self.params['out_path'])

        for slice_id, (slice, file) in enumerate(zip(self.slices, self.file_names)):
            slice.uns['slice_id'] = slice_id

            algo = SlidingWindowMultipleSizes(slice, slice_id, file, **self.params)
            # plot original annotation
            if self.params['plotting'] > 0:
                algo.plot_annotation()
            # run algorithm for feature extraction and cell type filtering based on entropy and scatteredness
            algo.run()
            if self.params['plotting'] > 1:
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
            if self.params['plotting'] > 3:
                algo.plot_celltype_images()
            # filter the cell types which are not localized using calculated metrics (entropy and scatteredness)
            algo.cell_type_filtering()
            # check if cell type filtering removes all cell types and raise an Error
            if algo.tissue.n_vars == 0:
                raise ValueError(r'Empty algo.tissue object. All cell types removed. \
                                 Adjust scatter_thres and entropy_thres to preserve useful cell types.')

            # add algo object for each slice to a list
            self.algo_list.append(algo)
        
        if self.params['plotting'] > 0 and len(self.algo_list) > 1:
            self.plot_all_annotation()

        # MERGE TISSUE ANNDATA
        # each tissue has slice_id as 3rd coordinate in tissue.obsm['spatial']
        merged_tissue = ad.concat([a.get_tissue() for a in self.algo_list], axis=0, join='outer')
        # if tissues have different sets of cell those columns are filled with NANs
        # this is corrected by writing 0s
        merged_tissue.X[np.isnan(merged_tissue.X)] = 0.0

        # CLUSTERING (WINDOW_LABELS)
        self.cluster(merged_tissue)

        for slice_id, algo in enumerate(self.algo_list):
            # extract clustering data from merged_tissue
            algo.set_clustering_labels(
                merged_tissue.obs.loc[merged_tissue.obsm['spatial'][:, 2] == slice_id, self.params['cluster_algo']])

            # COMMUNITY CALLING (MAJORITY VOTING)
            algo.community_calling()
            # copy final Cell Community Detection (CCD) result to original slices
            self.slices[slice_id].obs.loc[algo.adata.obs[f'tissue_{algo.method_key}'].index, 'cell_communities'] = algo.adata.obs[f'tissue_{algo.method_key}']

            # save anndata objects for further use
            if self.params['save_adata']:
                algo.save_anndata()
            algo.save_community_labels()
            algo.save_tissue()

            # PLOT COMMUNITIES & STATISTICS
            # plot cell communities clustering result
            if self.params['plotting'] > 0:
                algo.plot_clustering()

            # if flag skip_stats is active, skip cell mixture statistics analysis
            if not self.params['skip_stats']:
                algo.calculate_cell_mixture_stats()
                algo.save_mixture_stats()
                if self.params['plotting'] > 1:
                    algo.plot_stats()
                    algo.plot_celltype_table()
                if self.params['plotting'] > 2:
                    algo.plot_cluster_mixtures()
                    algo.boxplot_stats()
                if self.params['plotting'] > 4:
                    algo.colorplot_stats(color_system=self.params['color_plot_system'])
                    algo.colorplot_stats_per_cell_types()
                # save final tissue with stats
                algo.save_tissue(suffix='_stats')
        
        if self.params['plotting'] > 0 and len(self.algo_list) > 1:
            self.plot_all_clustering()
        if self.params['plotting'] > 2:
            self.plot_celltype_mixtures_total([algo.get_cell_mixtures().to_dict() for algo in self.algo_list])
            self.plot_cell_abundance_total()
            self.plot_cluster_abundance_total()
        if self.params['plotting'] > 3:
            self.plot_cell_abundance_per_slice()
            self.plot_cluster_abundance_per_slice()
            self.plot_cell_perc_in_community_per_slice()

        end_time = time.perf_counter()

        self.params['execution_time'] = end_time - start_time
        generate_report(self.params)
   
    @timeit
    def cluster(self, merged_tissue):
        """
        Perform clustering on merged tissue data from all slices.
        Supported clustering algorithms are:
        'leiden' - Leiden (scanpy) with neighbors similarity metric,
        'spectral' - Spectral (skimage) with neighbors similarity metric, and
        'agglomerative' - Agglomerative (skimage) with 'ward' linkage type
        and 'euclidian' distance metric.
        Cluster labels are stored in merged_tissue.obs[cluster_algo]
        and updated inplace.

        Parameters:
        - merged_tissue (AnnData): AnnData object containin features of all slices

        """
        if self.params['cluster_algo'] == 'leiden':
            sc.pp.neighbors(merged_tissue, use_rep='X')
            sc.tl.leiden(merged_tissue, resolution=self.params['resolution'])
        elif self.params['cluster_algo'] == 'spectral':
            sc.pp.neighbors(merged_tissue, use_rep='X')
            spcl = SpectralClustering(n_clusters=self.params['n_clusters'], eigen_solver='arpack', random_state=0, affinity='precomputed', n_jobs=5)
            merged_tissue.obs[self.params['cluster_algo']] = (spcl.fit_predict(merged_tissue.obsp['connectivities'])).astype('str')
        elif self.params['cluster_algo'] == 'agglomerative':
            ac = AgglomerativeClustering(n_clusters=self.params['n_clusters'], affinity='euclidean', compute_full_tree=False, linkage='ward', distance_threshold=None)
            merged_tissue.obs[self.params['cluster_algo']] = (ac.fit_predict(merged_tissue.X)).astype('str')
        else:
            logging.error('Unsupported clustering algorithm')
            raise ValueError("Unsupported clustering algorithm")
    
    def log_win_size_full_info(self):
        for slice, fname in zip(self.slices, self.file_names):
            x_min, x_max = np.min(slice.obsm['spatial'][:, 0]), np.max(slice.obsm['spatial'][:, 0])
            y_min, y_max = np.min(slice.obsm['spatial'][:, 1]), np.max(slice.obsm['spatial'][:, 1])
            x_range, y_range = abs(abs(x_max) - abs(x_min)), abs(abs(y_max) - abs(y_min))
            for win_size, slide_step in zip([int(w) for w in self.params['win_sizes'].split(',')], [int(s) for s in self.params['sliding_steps'].split(',')]):
                self.log_win_size_info_per_slice(slice, fname, win_size, slide_step, x_range, y_range)


    def log_win_size_info_per_slice(self, slice, fname, win_size, slide_step, x_range, y_range):
        """
        Logs window size information for a given slice.

        Parameters:
        - slice: The slice.
        - fname: The filename of the slice.
        - win_size: The size of the window.
        - slide_step: The sliding step for the window.
        - x_range: The range of x-coordinates.
        - y_range: The range of y-coordinates.

        """
        cell_to_loc = defaultdict(int)
        for x, y in slice.obsm['spatial']:
            cell_to_loc[(x // win_size, y // win_size)] += 1
        
        logging.info(f"""Window size info for slice: {fname}     
                     window size: {win_size}
                     sliding step: {slide_step}
                     cells mean: {np.mean(list(cell_to_loc.values())):.2f}
                     cells median: {np.median(list(cell_to_loc.values()))}
                     num horizontal windows: {int(x_range // win_size)}
                     num vertical windows: {int(y_range // win_size)}\n
                     """)

    @timeit
    def calc_optimal_win_size_and_slide_step(self):
        """
        Method for calculating the optimal window size and sliding step.
        Window size is calculated such that it covers between MIN_COVERED and MAX_COVERED cells.
        Sliding step is set to the half of the window size.

        """
        MAX_ITERS = 10
        MIN_COVERED = 30
        MAX_COVERED = 60
        AVG_COVERED_GOAL = (MAX_COVERED + MIN_COVERED) // 2
        
        x_min, x_max = np.min(self.slices[0].obsm['spatial'][:, 0]), np.max(self.slices[0].obsm['spatial'][:, 0])
        y_min, y_max = np.min(self.slices[0].obsm['spatial'][:, 1]), np.max(self.slices[0].obsm['spatial'][:, 1])
        x_range, y_range = abs(abs(x_max) - abs(x_min)), abs(abs(y_max) - abs(y_min))

        win_size = int(x_range // 50 if x_range < y_range else y_range // 50)
        delta_multiplier = win_size * 0.15

        avg_covered = -1
        delta = -1
        iters = 0
        while iters < MAX_ITERS:
            cell_to_loc = defaultdict(int)
            for x, y in self.slices[0].obsm['spatial']:
                cell_to_loc[(x // win_size, y // win_size)] += 1
            
            # using median instead of mean because many windows can be empty (space is not fully occupied by tissue)
            avg_covered = np.median(list(cell_to_loc.values()))
            
            if MIN_COVERED < avg_covered < MAX_COVERED:
                break

            delta = np.sign(AVG_COVERED_GOAL - avg_covered) * (
                AVG_COVERED_GOAL / avg_covered if AVG_COVERED_GOAL > avg_covered else avg_covered / AVG_COVERED_GOAL
                ) * delta_multiplier
            win_size += delta
            iters += 1
        
        #closest doubly even number so that sliding step is also even number
        win_size = round(win_size)
        win_size = win_size + ((win_size & 0b11) ^ 0b11) + 1 if win_size & 0b11 else win_size
        
        if iters == MAX_ITERS:
            logging.warn(f"Optimal window size not obtained in {MAX_ITERS} iterations.")
        self.log_win_size_info_per_slice(self.slices[0], self.file_names[0], win_size, win_size // 2, x_range, y_range)
        
        return (str(win_size), str(win_size // 2))
     
    def plot_all_slices(self, img_name, clustering=False):
        """
        Plot all slices using the specified algorithms and annotations.

        Parameters:
        - img_name (str): The name of the output image file.
        - clustering (bool, optional): Whether to plot clustering or cell type annotation. Defaults to False.

        """
        number_of_samples = len(self.algo_list)
        number_of_rows = 2 if number_of_samples % 2 == 0 and number_of_samples > 2 else 1
        number_of_columns = (number_of_samples // 2) if number_of_samples % 2 == 0 and number_of_samples > 2 else number_of_samples

        figure, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, squeeze=False, layout='constrained', figsize=(10,6))
        h_d = {}
        unknown_label = []
        for (algo, ax) in zip(self.algo_list, axes.flatten()):
            palette = algo.cluster_palette if clustering else algo.annotation_palette
            annotation = f'tissue_{self.algo_list[0].method_key}' if clustering else self.algo_list[0].annotation
            plot_spatial(algo.adata, annotation=annotation, palette=palette, spot_size=algo.spot_size, ax=ax)
            ax.get_legend().remove()
            ax.set_title(f'{algo.filename}', fontsize=6, loc='center', wrap=True)
            hands, labs = ax.get_legend_handles_labels()
            for h, l in zip(hands, labs):
                h._sizes = [11]
                if l == 'unknown':
                    unknown_label = np.array([[h, l]])
                    continue
                if l not in h_d.values():
                    h_d[h] = l
        try:
            handles = np.array([[h, int(l)] for (h, l) in h_d.items()])
        except:
            handles = np.array([[h, l] for (h, l) in h_d.items()])

        handles = handles[handles[:, 1].argsort()]
        handles[:, 1] = handles[:, 1].astype('str')

        if len(unknown_label)>0:
            handles = np.concatenate((handles, unknown_label), axis=0) 
        
        legend_ncols = 1 if len(handles) <= 12 else 2
        figure.legend(handles[:, 0], handles[:, 1], bbox_to_anchor=(1.15, 0.5), loc='center', fontsize=4, frameon=False, borderaxespad=0., ncol=legend_ncols, labelspacing=1)
        figure.savefig(f'{self.params["out_path"]}/{img_name}', dpi=self.params['dpi'], bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()

    @timeit
    def plot_all_annotation(self):
        self.plot_all_slices('cell_type_per_slice.png')

    @timeit
    def plot_all_clustering(self):
        self.plot_all_slices('clustering_per_slice.png', True)

    @timeit 
    def plot_celltype_mixtures_total(self, cell_mixtures):
        """
        Plot the total cell type mixtures.

        Parameters:
        - cell_mixtures (list): A list of dictionaries containing cell type mixtures.

        """
        def merge_dicts(dict1, dict2):
            return { cluster: dict1.get(cluster, 0) + dict2.get(cluster, 0) for cluster in set(dict1) | set(dict2) }
        def merge_dicts_of_dicts(dict1, dict2):
            return { celltype: merge_dicts(dict1.get(celltype, {}), dict2.get(celltype, {})) for celltype in set(dict1) | set(dict2) }

        total_dict = reduce(merge_dicts_of_dicts, cell_mixtures)
        total = pd.DataFrame(total_dict).fillna(0)

        total['total_counts'] = np.array([sum(total.loc[row, :]) for row in total.index]).astype(int)

        cell_type_counts = {ct:[int(sum(total[ct]))] for ct in total.columns}
        total = pd.concat([total, pd.DataFrame(cell_type_counts, index=['total_cells'])])

        total.iloc[:-1, :-1] = total.iloc[:-1, :-1].div(total['total_counts'][:-1], axis=0).mul(100)
        total['perc_of_all_cells'] = np.around(total['total_counts'] / total['total_counts'][-1] * 100, decimals=1)
        total = total.loc[sorted(total.index.values, key=lambda x: float(x) if x != "total_cells" else float('inf'))]

        set_figure_params(dpi=self.params['dpi'], facecolor='white')
        sns.set(font_scale=1.5)

        ncols = len(total.columns)
        fig, axes = plt.subplots(ncols=ncols, figsize=(30,20))
        fig.subplots_adjust(wspace=0)

        vmax_perc = np.max(np.ravel(total.iloc[:-1,:-2]))
        for i, ax in enumerate(axes[:-2]):
            sns.heatmap(pd.DataFrame(total.iloc[:, i]), vmin=0.0, vmax=vmax_perc, linewidths=0, linecolor=None, annot=True, cbar=False, ax=ax, \
                            cmap="Greys", fmt='4.0f', xticklabels=True, yticklabels=True if i==0 else False, square=True)
        sns.heatmap(pd.DataFrame(total.iloc[:, -2]), annot=True, vmin=0, vmax=np.max(total.iloc[:-1, -2]), linewidths=0, linecolor=None, \
            cbar=False, cmap='Greens', ax=axes[-2], fmt='4.0f', xticklabels=True, yticklabels=False, square=True)
        sns.heatmap(pd.DataFrame(total.iloc[:, -1]), annot=True, vmin=0, vmax=np.max(total.iloc[:-1, -1]), linewidths=0, linecolor=None, cbar=False, \
            cmap='Greens', ax=axes[-1], fmt='4.0f', xticklabels=True, yticklabels=False, square=True)
        
        for ax in axes:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=70)
            ax.xaxis.tick_top() 
        
        plt.savefig(os.path.join(self.params['out_path'], f'total_cell_mixtures_table.png'), bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()


    @timeit
    def plot_cell_perc_in_community_per_slice(self):
        """
        Plots the percentage of cells in each community per slice.
        """
        cells_in_comm_per_slice = {algo.filename: algo.get_community_labels().value_counts(normalize=True).rename(algo.filename) for algo in self.algo_list}
        df = pd.concat(cells_in_comm_per_slice.values(), axis=1).fillna(0).mul(100).T
        df = df[sorted(df.columns.values, key=lambda x: float(x) if x != "unknown" else float('inf'))]
        set_figure_params(dpi=self.params['dpi'], facecolor='white')
        sns.set(font_scale=1.5)
        plt.figure(figsize=(30,20))

        ax = sns.heatmap(df, annot=True, fmt="4.0f", cmap="Greys", xticklabels=True, yticklabels=True, square=True, cbar=False)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')

        plt.savefig(os.path.join(self.params['out_path'], 'cell_perc_in_community_per_slice.png'), bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()


    @timeit
    def plot_cell_abundance_total(self):
        """
        Plots the total cell abundance for each algorithm.
        """
        fig, ax = plt.subplots(figsize=(20,20))
        fig.subplots_adjust(wspace=0)
        set_figure_params(dpi=self.params['dpi'], facecolor='white')

        greys=cycle(['darkgray','gray','dimgray','lightgray'])
        colors = [next(greys) for _ in range(len(self.algo_list))]
        cell_percentage_dfs = []
        plot_columns = []
        for algo in self.algo_list:
            cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs[algo.annotation].value_counts(normalize=True).mul(100).rename(algo.filename)))
            plot_columns.append(algo.filename)

        cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0)
        cummulative_df.plot(y=plot_columns, kind="bar", rot=70, ax=ax, xlabel="", color=colors)

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        ax.grid(False)
        ax.set_facecolor('white')
        plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))

        plt.savefig(os.path.join(self.params['out_path'], f'cell_abundance_all_slices.png'), bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()


    @timeit
    def plot_cell_abundance_per_slice(self):
        """
        Plots the cell abundance for each algorithm per slice.
        """
        number_of_samples = len(self.algo_list)
        if number_of_samples <=2:
            number_of_rows = 1
            number_of_columns = number_of_samples
        else:
            number_of_rows = 2 if number_of_samples % 2 == 0 else 1
            number_of_columns = number_of_samples // 2 if number_of_samples % 2 == 0 else number_of_samples
        fig, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, figsize=(20,20), squeeze=False)
        axes = axes.ravel()
        fig.subplots_adjust(wspace=0)
        set_figure_params(dpi=self.params['dpi'], facecolor='white')

        cell_percentage_dfs = []
        plot_columns = []
        for algo in self.algo_list:
            cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs[algo.annotation].value_counts(normalize=True).mul(100).rename(algo.filename)))
            plot_columns.append(algo.filename)

        cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0)

        for i in range(number_of_rows * number_of_columns):
            axes[i].yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
            axes[i].set_facecolor('white')
            axes[i].set_title(plot_columns[i])
            cummulative_df.plot(y=plot_columns[i], kind="bar", rot=70, ax=axes[i], xlabel="", color="grey", legend=False)
            axes[i].grid(False)

        for ax in axes:
            ax.grid(False)

        plt.savefig(os.path.join(self.params['out_path'], f'cell_abundance_per_slice.png'), bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()

    @timeit 
    def plot_cluster_abundance_total(self):
        """
        Plots the total cluster abundance for each algorithm.
        """
        fig, ax = plt.subplots(figsize=(20,20))
        fig.subplots_adjust(wspace=0)
        set_figure_params(dpi=self.params['dpi'], facecolor='white')

        greys=cycle(['darkgray','gray','dimgray','lightgray'])
        colors = [next(greys) for _ in range(len(self.algo_list))]
        cell_percentage_dfs = []
        plot_columns = []
        for algo in self.algo_list:
            cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs[f'tissue_{algo.method_key}'].value_counts(normalize=True).mul(100).rename(algo.filename)))
            plot_columns.append(algo.filename)

        cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0)
        cummulative_df = cummulative_df.loc[sorted(cummulative_df.index.values, key=lambda x: float(x) if x != "unknown" else float('inf'))]
        cummulative_df.plot(y=plot_columns, kind="bar", rot=0, ax=ax, xlabel="", color=colors)

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        ax.grid(False)
        ax.set_facecolor('white')
        plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))

        plt.savefig(os.path.join(self.params['out_path'], f'cluster_abundance_all_slices.png'), bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()

    @timeit
    def plot_cluster_abundance_per_slice(self):
        """
        Plots the cluster abundance for each algorithm per slice.
        """
        number_of_samples = len(self.algo_list)
        if number_of_samples <= 2:
            number_of_rows = 1
            number_of_columns = number_of_samples
        else:
            number_of_rows = 2 if number_of_samples % 2 == 0 else 1
            number_of_columns = number_of_samples // 2 if number_of_samples % 2 == 0 else number_of_samples
        fig, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, figsize=(20,20), squeeze=False)
        axes = axes.ravel()
        fig.subplots_adjust(wspace=0)
        set_figure_params(dpi=self.params['dpi'], facecolor='white')

        cell_percentage_dfs = []
        plot_columns = []
        for algo in self.algo_list:
            cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs[f'tissue_{algo.method_key}'].value_counts(normalize=True).mul(100).rename(algo.filename)))
            plot_columns.append(algo.filename)

        cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0)
        cummulative_df = cummulative_df.loc[sorted(cummulative_df.index.values, key=lambda x: float(x) if x != "unknown" else float('inf'))]

        for i in range(number_of_rows * number_of_columns):
            axes[i].yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
            axes[i].set_facecolor('white')
            axes[i].set_title(plot_columns[i])
            cummulative_df.plot(y=plot_columns[i], kind="bar", rot=0, ax=axes[i], xlabel="", color="grey", legend=False)
            axes[i].grid(False)

        for ax in axes:
            ax.grid(False)
            
        plt.savefig(os.path.join(self.params['out_path'], f'cluster_abundance_per_slice.png'), bbox_inches='tight')
        if not self.params['hide_plots']:
            plt.show()
        plt.close()
    
    def plot(self, function_ind : str, slice_id = 0, community_id = None):
        if function_ind == "all_annotations": self.plot_all_annotation()
        if function_ind == "all_clustering": self.plot_all_clustering()
        if function_ind == "cell_type_mixtures_total": self.plot_celltype_mixtures_total([algo.get_cell_mixtures().to_dict() for algo in self.algo_list])
        if function_ind == "cell_perc_in_community_per_slice": self.plot_cell_perc_in_community_per_slice()
        if function_ind == "cell_abundance_total": self.plot_cell_abundance_total()
        if function_ind == "cell_abundance_per_slice": self.plot_cell_abundance_per_slice()
        if function_ind == "cluster_abundance_total": self.plot_cluster_abundance_total()
        if function_ind == "cluster_abundance_per_slice": self.plot_cluster_abundance_per_slice()

        if function_ind == "annotation": self.algo_list[slice_id].plot_annotation()
        if function_ind == "clustering": self.algo_list[slice_id].plot_clustering()
        if function_ind == "colorplot": self.algo_list[slice_id].colorplot_stats(self.params['color_plot_system'], community_id)
        if function_ind == "colorplot_cell_type": self.algo_list[slice_id].colorplot_stats_per_cell_types()
        if function_ind == "cell_types_table": self.algo_list[slice_id].plot_celltype_table()
        if function_ind == "boxplot": self.algo_list[slice_id].boxplot_stats(community_id)
        if function_ind == "cell_types_images": self.algo_list[slice_id].plot_celltype_images()
        if function_ind == "histogram_cell_sums": self.algo_list[slice_id].plot_histogram_cell_sum_window()
        if function_ind == "cluster_mixtures": self.algo_list[slice_id].plot_cluster_mixtures(community_id)
        if function_ind == "cell_mixture_table": self.algo_list[slice_id].plot_stats()


