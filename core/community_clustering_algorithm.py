import logging
import os
from itertools import cycle
from abc import ABC, abstractmethod

import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from .utils import timeit


class CommunityClusteringAlgo(ABC):
    def __init__(self, adata, slice_id, input_file_path, **params):
        sc.settings.verbosity = 3 if params['verbose'] else params['verbose']
        sc.settings.set_figure_params(dpi=300, facecolor='white')
        self.adata = adata
        self.slice_id = slice_id
        self.adata.uns['algo_params'] = params
        self.adata.uns['sample_name'] = os.path.basename(input_file_path.rsplit(".", 1)[0])
        for key, value in params.items():
            setattr(self, key, value)

        self.tissue = None

        cell_count_limit = (self.min_count_per_type*len(self.adata)) // 100
        cell_over_limit = []
        for cell_tp in self.adata.obs[self.annotation].cat.categories:
            cell_num = sum(self.adata.obs[self.annotation]==cell_tp)
            if cell_num > cell_count_limit:
                cell_over_limit.append(cell_tp)
            else:
                logging.info(f'{cell_tp} cell type excluded, due to insufficient cells of that type: {cell_num} cells < {int(cell_count_limit)} ({self.min_count_per_type} % of {len(self.adata)})')
        
        self.adata = self.adata[self.adata.obs[self.annotation].isin(cell_over_limit),:]
        self.unique_cell_type = list(self.adata.obs[self.annotation].cat.categories)

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def calc_feature_matrix(self):
        pass

    @abstractmethod
    def community_calling(self):
        pass

    def get_tissue(self):
        return self.tissue
    
    def set_clustering_labels(self, labels):
        self.tissue.obs.loc[:, 'leiden'] = labels

    def plot_annotation(self):
        figure, ax = plt.subplots(nrows=1, ncols=1)
        sc.pl.spatial(self.adata, color=[self.annotation], palette=None, spot_size=self.spot_size, ax=ax, show=False, frameon=False)
        figure.savefig(os.path.join(self.dir_path, f'cell_type_annotation.png'), dpi=300, bbox_inches='tight')
        plt.close()

    def plot_histogram_cell_sum_window(self):
        figure, ax = plt.subplots(nrows=1, ncols=1)
        plt.hist(self.tissue.obs['window_cell_sum'].values)
        figure.savefig(os.path.join(self.dir_path, f'window_cell_num_hist_ws_{"_".join([str(i) for i in self.win_sizes_list])}.png'), dpi=300, bbox_inches='tight')
        plt.close()
   
    def cell_type_filtering(self):
        # extract binary image of cell positions for each cell type in the slice
        var_use = self.tissue.var.loc[(self.tissue.var['entropy']<self.entropy_thres) & (self.tissue.var['scatteredness']<self.scatter_thres)].index
        self.tissue.raw = self.tissue
        self.tissue = self.tissue[:, var_use]

    def plot_celltype_images(self):
        for cell_t in self.unique_cell_type:
            plt.imsave(fname=os.path.join(self.dir_path, f'tissue_window_{cell_t}_{self.params_suffix}.png'), arr=self.tissue.uns['cell_t_images'][cell_t], vmin=0, vmax=1, cmap='gray', dpi=250)
    
    @timeit
    def plot_clustering(self):
        # # plot initial clustering for each window
        # sc.pl.spatial(self.tissue, color='leiden', spot_size=1)
        # # plot clustering after majority voting for each subwindow
        # sc.pl.spatial(self.tissue, color='leiden_max_vote', spot_size=1)
        figure, ax = plt.subplots(nrows=1, ncols=1)
        sc.pl.spatial(self.adata, color=[f'tissue_{self.method_key}'], palette=None, spot_size=self.spot_size, ax=ax, show=False, frameon=False)
        figure.savefig(os.path.join(self.dir_path, f'clusters_cellspots_{self.params_suffix}.png'), dpi=300, bbox_inches='tight')
        plt.close()

    def calculate_spatial_cell_type_metrics(self):
        pass
    
    def calculate_cell_mixture_stats(self):

        # extract information on leiden clustering labels and cell types to create cell communities statistics
        clustering_labels = f'tissue_{self.method_key}'
        cell_types_communities = self.adata.obs[[clustering_labels, self.annotation]]

        stats_table = {}
        # calculate cell type mixtures for every cluster
        for label, cluster_data in cell_types_communities.groupby(clustering_labels):
            cell_type_dict = {ct:0 for ct in self.unique_cell_type}
            for cell in cluster_data[self.annotation]:
                cell_type_dict[cell]+=1
            # remove excluded cell types
            cell_type_dict = {k:cell_type_dict[k] for k in self.tissue.var.index.sort_values()}

            stats_table[label] = {k:cell_type_dict[k] for k in cell_type_dict}

            stats_table[label]['total_counts'] = int(sum(cell_type_dict.values()))

        stats = pd.DataFrame(stats_table).T
        stats.columns.name = "cell types"

        # add final row with total counts per cell types
        cell_type_counts = {ct:[int(sum(stats[ct].values))] for ct in self.tissue.var.index}
        stats = pd.concat([stats, pd.DataFrame(cell_type_counts, index=['total_cells'])])

        # divide each row with total sum of cells per cluster
        for i in range(len(stats.index.values[:-1])):
            if stats.iloc[i,-1] > 0:
                stats.iloc[i, :-1] = (100 * stats.iloc[i, :-1] / stats.iloc[i, -1]).astype(int)

        # add column with percentage of all cells belonging to a cluster
        stats['perc_of_all_cells'] = np.around(stats['total_counts'] / stats['total_counts'].sum() * 100, decimals=1)

        # save cell mixture statistics to tissue
        self.tissue.uns['cell mixtures'] = stats.iloc[:, :]

    def plot_stats(self):
        stats = self.tissue.uns['cell mixtures']
        sc.settings.set_figure_params(dpi=400, facecolor='white')
        sns.set(font_scale=0.5)

        ncols = len(stats.columns) # we want to separately print the total_counts column
        fig, axes = plt.subplots(ncols=ncols)

        # no space between columns
        fig.subplots_adjust(wspace=0)

        # put colormaps of your choice in a list:
        cmap_cycle = cycle(['Greens', 'Reds', 'Blues', 'Oranges', 'Purples'])

        for i, ax in enumerate(axes):
            g = sns.heatmap(pd.DataFrame(stats.iloc[:, i]), vmin=0.0, vmax=50, linewidths=0, linecolor=None, annot=True, cbar=False, ax=ax, \
                            cmap=cmap_cycle.__next__(),fmt='4.0f', xticklabels=True, yticklabels=True if i==0 else False, square=True)
            g.set_xticklabels(g.get_xticklabels(), rotation=70)
            g.xaxis.tick_top() # x axis on top
        # final column should have the sum of all cells per cluster
        sns.heatmap(pd.DataFrame(stats.iloc[:, -1]), annot=True, linewidths=0, linecolor=None, cbar=False, cmap=None, ax=ax, \
                    fmt='4.0f', xticklabels=True, yticklabels=False, square=True)
        plt.savefig(os.path.join(self.dir_path, f'cell_mixture_table_{self.params_suffix}.png'), dpi=400)
        plt.close()

        # plot each cluster and its cells mixture
        sc.settings.set_figure_params(dpi=100, facecolor='white')

        new_stats = stats.copy()
        new_stats = new_stats.drop(labels=['total_counts', 'perc_of_all_cells'], axis=1)
        new_stats = new_stats.drop(labels='total_cells', axis=0)
        for cluster in new_stats.iterrows():
            ct_perc = cluster[1].sort_values(ascending=False)
            # only display clusters with more than min_cells_in_cluster cells
            if stats.loc[cluster[0]]['total_counts'] > self.min_cluster_size:
                # only cell types which have more than min_perc_to_show abundance will be shown
                ct_ind = [x for x in ct_perc.index[ct_perc>self.min_perc_to_show]]
                
                fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,6))
                fig.subplots_adjust(wspace=0.35)

                sc.pl.spatial(self.adata, groups=ct_ind, color=self.annotation, spot_size=self.spot_size, ax=ax[0], show=False, frameon=False)
                ax[0].legend([f'{ind.get_text()} ({ct_perc[ind.get_text()]})' for ind in ax[0].get_legend().texts[:-1]], bbox_to_anchor=(1.0, 0.5), loc='center left', frameon=False, fontsize=12)
                sc.pl.spatial(self.adata, groups=[cluster[0]], color=f'tissue_{self.method_key}', spot_size=self.spot_size, ax=ax[1], show=False, frameon=False)
                ax[1].legend([f'{ind.get_text()} ({stats.loc[ind.get_text(), "perc_of_all_cells"]})' for ind in ax[1].get_legend().texts[:-1]], bbox_to_anchor=(1.0, 0.5), loc='center left', frameon=False, fontsize=12)
                fig.savefig(os.path.join(self.dir_path, f'cmixtures_{self.params_suffix}_c{cluster[0]}.png'), bbox_inches='tight')

                plt.close()

    def boxplot_stats(self):
        stats = self.tissue.uns['cell mixtures']
        
        # box plot per cluster of cell type percentages distribution
        sc.settings.set_figure_params(dpi=100, facecolor='white')

        # drop total_counts, perc_of_all_cells columns and total_cells row
        new_stats = stats.copy()
        new_stats = new_stats.drop(labels=['total_counts', 'perc_of_all_cells'], axis=1)
        new_stats = new_stats.drop(labels='total_cells', axis=0)
        for cluster in new_stats.iterrows():


            win_cell_distrib = self.tissue[self.tissue.obs['leiden'] == cluster[0]] # this is going to be more noisy because
            # the feature vector contains the cell mixture of the window which top left corner is in this subwindow
            # thus, the label of the subwindow is correlated with the feature vector of bin_ratio times bigger window
            # covering the the down right area from the subwindow.
            win_cell_distrib_df = pd.DataFrame(win_cell_distrib.X / (self.total_cell_norm/100), columns=win_cell_distrib.var.index)

            # Reshape data into long format
            cell_type_distrib = pd.melt(win_cell_distrib_df, var_name='Cell Type', value_name='Percentage')
            cell_type_distrib = cell_type_distrib.sort_values(by='Cell Type')
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,6))
            # plot boxplot of cell type percentages per mixture
            ax = sns.boxplot(x='Cell Type', y='Percentage', data=cell_type_distrib)
            # overlap with a plot of specific percentage values. 
            # Jitter allows dots to move left and right for better visibility of all points
            ax = sns.stripplot(x='Cell Type', y='Percentage', data=cell_type_distrib, jitter=True, color='black', size=2)
            # remove top and right frame of the plot
            sns.despine(top=True, right=True)
            ax.set_title(f'mixutre {cluster[0]} ({self.adata.uns["sample_name"]})')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.xaxis.tick_bottom() # x axis on the bottom
            fig.savefig(os.path.join(self.dir_path, f'boxplot_{self.params_suffix}_c{cluster[0]}.png'), bbox_inches='tight')

            plt.close()
            
    def save_metrics(self):
        # save metrics results in csv format
        self.tissue.var[['entropy', 'scatteredness']].to_csv(os.path.join(self.dir_path, f'spatial_metrics_{self.params_suffix}.csv'))

    def save_tissue(self, suffix=''):
        # save anndata file
        self.tissue.write_h5ad(os.path.join(self.dir_path, f'tissue_{self.filename}{suffix}.h5ad'), compression="gzip")

        logging.info(f'Saved clustering result tissue_{self.filename}{suffix}.h5ad.')

    def save_adata(self, suffix=''):
        # save anndata file
        self.adata.write_h5ad(os.path.join(self.dir_path, f'{self.filename}{suffix}.h5ad'), compression="gzip")

        logging.info(f'Saved clustering result as a part of original anndata file {self.filename}{suffix}.h5ad.')

    def save_mixture_stats(self):
        # save cell mixture statistics
        self.tissue.uns['cell mixtures'].to_csv(os.path.join(self.dir_path, f'cell_mixture_stats_{self.params_suffix}.csv'))