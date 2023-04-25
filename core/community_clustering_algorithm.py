import logging
import os

import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from matplotlib import pyplot as plt
# import skimage.measure
# import scipy.ndimage.measurements
from itertools import cycle

from abc import ABC, abstractmethod

class CommunityClusteringAlgo(ABC):
    def __init__(self, adata, slice_id, input_file_path, **params):
        self.adata = adata
        self.slice_id = slice_id
        self.adata.uns['algo_params'] = params
        self.adata.uns['sample_name'] = os.path.basename(input_file_path.rsplit(".", 1)[0])
        for key, value in params.items():
            setattr(self, key, value)

        self.unique_cell_type = list(self.adata.obs[self.annotation].cat.categories)
        self.tissue = None
        # [NOTE] this should be included later
        # # define if each cell type needs a minimum amount of cells to be considered in cell mixtures and what is the minimum value
        # min_count_per_type_limit = False
        # min_count_per_type = 100

        # # define if only cell types with more than min_count_per_type cells are used for cell communities extraction
        # if min_count_per_type_limit:

        #     cell_over_limit = []
        #     for cell_tp in adata.obs[annotation_label].cat.categories:
        #         cell_num = sum(adata.obs[annotation_label]==cell_tp)
        #         if cell_num > min_count_per_type:
        #             cell_over_limit.append(cell_tp)

        #     adata = adata[adata.obs[annotation_label].isin(cell_over_limit),:]

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
        # to prevent warning about appending data to a view of obs
        self.tissue.obs = self.tissue.obs.copy()
        self.tissue.obs['leiden'] = pd.Series(labels, index=self.tissue.obs.index)
   
    def cell_type_filtering(self):
        # extract binary image of cell positions for each cell type in the slice
        var_use = self.tissue.var.loc[(self.tissue.var['entropy']<self.entropy_thres) & (self.tissue.var['scatteredness']<self.scatter_thres)].index
        self.tissue.raw = self.tissue
        self.tissue = self.tissue[:, var_use]
    
    def plot_celltype_images(self):
        for cell_t in self.unique_cell_type:
            plt.imsave(fname=os.path.join(self.dir_path, f'tissue_window_{cell_t}_{self.params_suffix}.png'), arr=self.tissue.uns['cell_t_images'][cell_t], vmin=0, vmax=1, cmap='gray', dpi=250)

    def plot_clustering(self):
        # # plot initial clustering for each window
        # sc.pl.spatial(self.tissue, color='leiden', spot_size=1)
        # # plot clustering after majority voting for each subwindow
        # sc.pl.spatial(self.tissue, color='leiden_max_vote', spot_size=1)

        figure, ax = plt.subplots(nrows=1, ncols=1)
        sc.pl.spatial(self.adata, color=[f'tissue_{self.method_key}'], spot_size=self.spot_size, ax=ax, show=False)
        ax.axis('off')
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
            cell_type_dict = {k:cell_type_dict[k] for k in self.tissue.var.index}
            stats_table[label] = {k:cell_type_dict[k] for k in cell_type_dict}

            stats_table[label]['total_counts'] = int(sum(cell_type_dict.values()))

        stats = pd.DataFrame(stats_table).T
        # stats.index.name='clusters'
        stats.columns.name="cell types"

        # add final row with total counts per cell types
        cell_type_counts = {ct:[int(sum(stats[ct].values))] for ct in self.tissue.var.index}
        stats = pd.concat([stats, pd.DataFrame(cell_type_counts, index=['total_cells'])])

        # divide each row with total sum of cells per cluster
        for i in range(len(stats.index.values[:-1])):
            if stats.iloc[i,-1] > 0:
                stats.iloc[i, :-1] = (100 * stats.iloc[i, :-1] / stats.iloc[i, -1]).astype(int)
        # save cell mixture statistics to csv file and to tissue
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
            g = sns.heatmap(pd.DataFrame(stats.iloc[:, i]), vmin=0.0, vmax=50, linewidths=0, linecolor=None, annot=True, cbar=False, ax=ax, cmap=cmap_cycle.__next__(),fmt='4.0f', xticklabels=True, yticklabels=True if i==0 else False, square=True)
            g.set_xticklabels(g.get_xticklabels(), rotation=70)
            g.xaxis.tick_top() # x axis on top
        # final column should have the sum of all cells per cluster
        sns.heatmap(pd.DataFrame(stats.iloc[:, -1]), annot=True, linewidths=0, linecolor=None, cbar=False, cmap=None, ax=ax, fmt='4.0f', xticklabels=True, yticklabels=False, square=True)
        plt.savefig(os.path.join(self.dir_path, f'cell_mixture_table_{self.params_suffix}.png'), dpi=400)
        plt.close()

        min_cell_types = 3
        min_perc = 15
        min_perc_to_show = 5
        min_cells_in_cluster = 500

        sc.settings.set_figure_params(dpi=100, facecolor='white')

        new_stats = stats.copy()
        new_stats = new_stats.drop(labels='total_counts', axis=1)
        new_stats = new_stats.drop(labels='total_cells', axis=0)
        for cluster in new_stats.iterrows():
            ct_perc = cluster[1].sort_values(ascending=False)
            if ct_perc[min_cell_types-1] > min_perc and stats.loc[cluster[0]]['total_counts']>min_cells_in_cluster:
                ct_ind = [x for x in ct_perc.index[ct_perc>min_perc_to_show]]
                
                fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,6))
                fig.subplots_adjust(wspace=0.35)

                sc.pl.spatial(self.adata, groups=ct_ind, color=self.annotation, spot_size=self.spot_size, ax=ax[0], show=False, frameon=False)
                ax[0].legend([f'{ind.get_text()} ({ct_perc[ind.get_text()]})' for ind in ax[0].get_legend().texts[:-1]], bbox_to_anchor=(1.0, 0.5), loc='center left', frameon=False, fontsize=12)
                sc.pl.spatial(self.adata, groups=[cluster[0]], color=f'tissue_{self.method_key}', spot_size=self.spot_size, ax=ax[1], show=False, frameon=False,)
                fig.savefig(os.path.join(self.dir_path, f'cmixtures_{self.params_suffix}_c{cluster[0]}.png'), bbox_inches='tight')

                plt.close()

    def save_metrics(self):
        # save metrics results in csv format
        # print(self.tissue.var[['entropy', 'scatteredness']])
        self.tissue.var[['entropy', 'scatteredness']].to_csv(os.path.join(self.dir_path, f'spatial_metrics_{self.params_suffix}.csv'))

    def save_tissue(self, suffix=''):
        # save anndata file
        self.tissue.write_h5ad(os.path.join(self.dir_path, f'tissue_{self.filename}{suffix}.h5ad'), compression="gzip")

        logging.info(f'Saved clustering result tissue_{self.filename}.h5ad.')

    def save_mixture_stats(self):
        # save cell mixture statistics
        self.tissue.uns['cell mixtures'].to_csv(os.path.join(self.dir_path, f'cell_mixture_stats_{self.params_suffix}.csv'))