import anndata as ad
import numpy as np
import matplotlib.ticker as mtick
import os
import pandas as pd
import seaborn as sns
import scanpy as sc
import random
import logging

from functools import reduce
from itertools import cycle
from matplotlib import pyplot as plt
from collections import defaultdict

from anndata import AnnData
from typing import List
from ccd import *


class CommunityDetection():
    def __init__(
            self,
            slices: List[AnnData],
            annotation: str,
            file_names: List[str] = None,
            **kwargs) -> None:
        self.params = { **COMMUNITY_DETECTION_DEFAULTS, **kwargs }
        self.params['annotation'] = annotation
        self.slices = slices
        self.file_names = file_names if file_names != None else [f"Slice_{id}" for id in range(len(slices))]
        if self.params['win_sizes'] == 'NA' or self.params['sliding_steps'] == 'NA':
            self.params['win_sizes'], self.params['sliding_steps'] = self.calc_optimal_win_size_and_slide_step()

    @timeit
    def run(self):
        if not os.path.exists(self.params['out_path']):
            os.makedirs(self.params['out_path'])

        algo_list = []
        win_sizes = "_".join([i for i in self.params['win_sizes'].split(',')])
        self.params['project_name_orig'] = self.params['project_name']
        self.params['out_path_orig'] = self.params['out_path']
        self.params['project_name'] += f"_r{self.params['resolution']}_ws{win_sizes}_en{self.params['entropy_thres']}_sct{self.params['scatter_thres']}_dwr{self.params['downsample_rate']}_mcc{self.params['min_cells_coeff']}"
        self.params['out_path'] = os.path.join(self.params['out_path'], self.params['project_name'])

        if not os.path.exists(self.params['out_path']):
            os.makedirs(self.params['out_path'])

        for slice_id, (slice, file) in enumerate(zip(self.slices, self.file_names)):
            slice.uns['slice_id'] = slice_id
            if 'X_spatial' in slice.obsm:
                slice.obsm['spatial'] = slice.obsm['X_spatial'].copy()
            elif 'spatial_stereoseq' in slice.obsm:
                slice.obsm['spatial'] = np.array(slice.obsm['spatial_stereoseq'].copy())

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

            # add algo object for each slice to a list
            algo_list.append(algo)
        
        if self.params['plotting'] > 0 and len(algo_list) > 1:
            self.plot_all_annotation(self.params['out_path'], algo_list)

        # MERGE TISSUE ANNDATA
        # each tissue has slice_id as 3rd coordinate in tissue.obsm['spatial']
        merged_tissue = ad.concat([a.get_tissue() for a in algo_list], axis=0, join='outer')
        # if tissues have different sets of cell those columns are filled with NANs
        # this is corrected by writing 0s
        merged_tissue.X[np.isnan(merged_tissue.X)] = 0.0

        # CLUSTERING (WINDOW_LABELS)
        sc.pp.neighbors(merged_tissue, use_rep='X')
        sc.tl.leiden(merged_tissue, resolution=self.params['resolution'])

        for slice_id, algo in enumerate(algo_list):
            # extract clustering data from merged_tissue
            algo.set_clustering_labels(
                merged_tissue.obs.loc[merged_tissue.obsm['spatial'][:, 2] == slice_id, 'leiden'])

            # COMMUNITY CALLING (MAJORITY VOTING)
            algo.community_calling()

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
                    algo.colorplot_stats(color_system=self.params['color_plot_system'])
                if self.params['plotting'] > 3:
                    algo.colorplot_stats_per_cell_types()
                # save final tissue with stats
                algo.save_tissue(suffix='_stats')
        
        if self.params['plotting'] > 0 and len(algo_list) > 1:
            self.plot_all_clustering(self.params['out_path'], algo_list)
        if self.params['plotting'] > 2:
            self.plot_celltype_mixtures_total([algo.get_cell_mixtures().to_dict() for algo in algo_list], self.params['out_path'])
            self.plot_cell_abundance_total(algo_list, self.params['out_path'])
            self.plot_cluster_abundance_total(algo_list, self.params['out_path'])
        if self.params['plotting'] > 3:
            self.plot_cell_abundance_per_slice(algo_list, self.params['out_path'])
            self.plot_cluster_abundance_per_slice(algo_list, self.params['out_path'])
        if self.params['plotting'] > 4:
            self.plot_cell_perc_in_community_per_slice(algo_list, self.params['out_path'])

        generate_report(self.params)

    @timeit
    def calc_optimal_win_size_and_slide_step(self):
        """
        Method for calculating the optimal window size and sliding step.
        Window size is calculated such that it covers between MIN_COVERED and MAX_COVERED cells.
        Sliding step is set to the half of the window size.

        """
        MAX_ITERS = 10
        MIN_COVERED = 30
        MAX_COVERED = 50
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
            
            num_selected = int(len(cell_to_loc) * 0.1)
            avg_covered = np.median(random.choices(list(cell_to_loc.values()), k=num_selected))
            
            if MIN_COVERED < avg_covered < MAX_COVERED:
                break

            delta = np.sign(AVG_COVERED_GOAL - avg_covered) * (
                AVG_COVERED_GOAL / avg_covered if AVG_COVERED_GOAL > avg_covered else avg_covered / AVG_COVERED_GOAL
                ) * delta_multiplier
            win_size += delta
            iters += 1
        
        #closest doubly even number so that sliding step is also even number
        win_size = int(win_size)
        win_size = win_size + ((win_size & 0b11) ^ 0b11) + 1 if win_size & 0b11 else win_size
        
        return (str(win_size), str(win_size // 2))
    
    def plot_all_slices(self, out_path, algo_list, annotation, img_name, clustering=False):
        """
        Plot all slices using the specified algorithms and annotations.

        Parameters:
        - out_path (str): The output path where the image will be saved.
        - algo_list (list): A list of algorithm objects to plot.
        - annotation (str): The annotation to use for coloring.
        - img_name (str): The name of the output image file.
        - clustering (bool, optional): Whether clustering is enabled. Defaults to False.

        """
        number_of_samples = len(algo_list)
        number_of_rows = 2 if number_of_samples % 2 == 0 and number_of_samples > 2 else 1
        number_of_columns = (number_of_samples // 2) if number_of_samples % 2 == 0 and number_of_samples > 2 else number_of_samples

        figure, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, squeeze=False, layout='constrained')
        h_d = {}
        unknown_label = []
        for (algo, ax) in zip(algo_list, axes.flatten()):
            palette = algo.cluster_palette if clustering else algo.annotation_palette
            plot_spatial(algo.adata, color=[annotation], palette=palette, spot_size=algo.spot_size, ax=ax, show=False, frameon=False)
            ax.get_legend().remove()
            ax.set_title(f'{algo.filename}', fontsize=6, loc='center', wrap=True)
            hands, labs = ax.get_legend_handles_labels()
            for h, l in zip(hands, labs):
                h._sizes = [11]
                if l=='unknown':
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
        figure.legend(handles[:, 0], handles[:, 1], bbox_to_anchor=(1.15, 0.5), loc='center', fontsize=4, frameon=False, borderaxespad=0., ncol=legend_ncols, labelspacing=1, scatterpoints=10)
        figure.savefig(f'{out_path}/{img_name}', dpi=150, bbox_inches='tight')
        plt.close()

    @timeit
    def plot_all_annotation(self, out_path, algo_list):
        self.plot_all_slices(out_path, algo_list, algo_list[0].annotation, 'cell_type_per_slice.png')

    @timeit
    def plot_all_clustering(self, out_path, algo_list):
        self.plot_all_slices(out_path, algo_list, f'tissue_{algo_list[0].method_key}', 'clustering_per_slice.png', True)

    @timeit 
    def plot_celltype_mixtures_total(self, cell_mixtures, path):
        """
        Plot the total cell type mixtures.

        Parameters:
        - cell_mixtures (list): A list of dictionaries containing cell type mixtures.
        - path (str): The path where the plot will be saved.

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

        sc.settings.set_figure_params(dpi=300, facecolor='white')
        sns.set(font_scale=1.5)

        ncols = len(total.columns)
        fig, axes = plt.subplots(ncols=ncols, figsize=(20,15))
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
        
        plt.tight_layout()
        plt.savefig(os.path.join(path, f'total_cell_mixtures_table.png'), bbox_inches='tight')
        plt.close()


    @timeit
    def plot_cell_perc_in_community_per_slice(self, algos, path):
        """
        Plots the percentage of cells in each community per slice.

        Parameters:
        - algos (list): A list of algorithms.
        - path (str): The path to save the plot.
        """
        cells_in_comm_per_slice = {algo.filename: algo.get_community_labels().value_counts(normalize=True).rename(algo.filename) for algo in algos}
        df = pd.concat(cells_in_comm_per_slice.values(), axis=1).fillna(0).mul(100).T
        df = df[sorted(df.columns.values, key=lambda x: float(x) if x != "unknown" else float('inf'))]
        sc.settings.set_figure_params(dpi=200, facecolor='white')
        sns.set(font_scale=1.5)
        plt.figure(figsize=(15, 15))

        ax = sns.heatmap(df, annot=True, fmt="4.0f", cmap="Greys", xticklabels=True, yticklabels=True, square=True, cbar=False)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        plt.tight_layout()
        plt.savefig(os.path.join(path, 'cell_perc_in_community_per_slice.png'), bbox_inches='tight')
        plt.close()


    @timeit
    def plot_cell_abundance_total(self, algos, path):
        """
        Plots the total cell abundance for each algorithm.

        Parameters:
        - algos (list): A list of algorithms.
        - path (str): The path to save the plot.
        """
        fig, ax = plt.subplots(figsize=(20,10))
        fig.subplots_adjust(wspace=0)
        sc.settings.set_figure_params(dpi=300, facecolor='white')

        greys=cycle(['darkgray','gray','dimgray','lightgray'])
        colors = [next(greys) for _ in range(len(algos))]
        cell_percentage_dfs = []
        plot_columns = []
        for algo in algos:
            cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs[algo.annotation].value_counts(normalize=True).mul(100).rename(algo.filename)))
            plot_columns.append(algo.filename)

        cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0)
        cummulative_df.plot(y=plot_columns, kind="bar", rot=70, ax=ax, xlabel="", color=colors)

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        ax.grid(False)
        ax.set_facecolor('white')
        plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
        plt.tight_layout()
        plt.savefig(os.path.join(path, f'cell_abundance_all_slices.png'))
        plt.close()


    @timeit
    def plot_cell_abundance_per_slice(self, algos, path):
        """
        Plots the cell abundance for each algorithm per slice.

        Parameters:
        - algos (list): A list of algorithms.
        - path (str): The path to save the plot.
        """
        number_of_samples = len(algos)
        if number_of_samples <=2:
            number_of_rows = 1
            number_of_columns = number_of_samples
        else:
            number_of_rows = 2 if number_of_samples % 2 == 0 else 1
            number_of_columns = number_of_samples // 2 if number_of_samples % 2 == 0 else number_of_samples
        fig, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, figsize=(20,10), squeeze=False)
        axes = axes.ravel()
        fig.subplots_adjust(wspace=0)
        sc.settings.set_figure_params(dpi=300, facecolor='white')

        cell_percentage_dfs = []
        plot_columns = []
        for algo in algos:
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
        plt.tight_layout()
        plt.savefig(os.path.join(path, f'cell_abundance_per_slice.png'))
        plt.close()

    @timeit 
    def plot_cluster_abundance_total(self, algos, path):
        """
        Plots the total cluster abundance for each algorithm.

        Parameters:
        - algos (list): A list of algorithms.
        - path (str): The path to save the plot.
        
        """
        fig, ax = plt.subplots(figsize=(20,10))
        fig.subplots_adjust(wspace=0)
        sc.settings.set_figure_params(dpi=300, facecolor='white')

        greys=cycle(['darkgray','gray','dimgray','lightgray'])
        colors = [next(greys) for _ in range(len(algos))]
        cell_percentage_dfs = []
        plot_columns = []
        for algo in algos:
            cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs[f'tissue_{algo.method_key}'].value_counts(normalize=True).mul(100).rename(algo.filename)))
            plot_columns.append(algo.filename)

        cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0)
        cummulative_df = cummulative_df.loc[sorted(cummulative_df.index.values, key=lambda x: float(x) if x != "unknown" else float('inf'))]
        cummulative_df.plot(y=plot_columns, kind="bar", rot=0, ax=ax, xlabel="", color=colors)

        ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        ax.grid(False)
        ax.set_facecolor('white')
        plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
        plt.tight_layout()
        plt.savefig(os.path.join(path, f'cluster_abundance_all_slices.png'))
        plt.close()

    @timeit
    def plot_cluster_abundance_per_slice(self, algos, path):
        """
        Plots the cluster abundance for each algorithm per slice.

        Parameters:
        - algos (list): A list of algorithms.
        - path (str): The path to save the plot.
        
        """
        number_of_samples = len(algos)
        if number_of_samples <= 2:
            number_of_rows = 1
            number_of_columns = number_of_samples
        else:
            number_of_rows = 2 if number_of_samples % 2 == 0 else 1
            number_of_columns = number_of_samples // 2 if number_of_samples % 2 == 0 else number_of_samples
        fig, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, figsize=(20, 10), squeeze=False)
        axes = axes.ravel()
        fig.subplots_adjust(wspace=0)
        sc.settings.set_figure_params(dpi=100, facecolor='white')

        cell_percentage_dfs = []
        plot_columns = []
        for algo in algos:
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
        plt.tight_layout()
        plt.savefig(os.path.join(path, f'cluster_abundance_per_slice.png'))
        plt.close()

