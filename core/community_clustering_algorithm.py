import logging
import os
from itertools import cycle
from abc import ABC, abstractmethod

from skimage import color
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.colors as mcolors
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

        self.unique_cell_type = list(self.adata.obs[self.annotation].cat.categories)
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
    ## CALCULATE_CELL_MIXTURE_STATS
    # brief - Creates a pandas DataFrame of cell type percentages per class
    # percentages are calculated globaly for all cells with single class label.
    # This is saved in adata.uns['cell mixtures'] for further use by plot fn.
    # Columns of total cell count per class and percentage of tissue per cluster
    # are added. Row of total cell type count is added.
    # DataFrame with additional columns and row is saved in adata.uns['cell mixture stats']
    def calculate_cell_mixture_stats(self):

        # extract information on leiden clustering labels and cell types to create cell communities statistics
        clustering_labels = f'tissue_{self.method_key}'
        cell_types_communities = self.adata.obs[[clustering_labels, self.annotation]]
        # # remove cells with unknown cell community label
        # cell_types_communities = cell_types_communities[self.adata.obs[clustering_labels] != 'unknown']

        stats_table = {}
        # calculate cell type mixtures for every cluster
        for label, cluster_data in cell_types_communities.groupby(clustering_labels):
            cell_type_dict = {ct:0 for ct in self.unique_cell_type}
            for cell in cluster_data[self.annotation]:
                cell_type_dict[cell]+=1
            # remove excluded cell types
            cell_type_dict = {k:cell_type_dict[k] for k in self.tissue.var.index.sort_values()}

            stats_table[label] = {k:cell_type_dict[k] for k in cell_type_dict}

        stats = pd.DataFrame(stats_table).T
        stats.columns.name = "cell types"

        # save absolute cell mixtures to tissue
        self.tissue.uns['cell mixtures'] = stats.iloc[:,:].copy()

        # add column with total cell count per cluster
        stats['total_counts'] = np.array([sum(stats.loc[row, :]) for row in stats.index]).astype(int)

        # add row with total counts per cell types
        cell_type_counts = {ct:[int(sum(stats[ct]))] for ct in self.tissue.var.index}
        stats = pd.concat([stats, pd.DataFrame(cell_type_counts, index=['total_cells'])])

        # divide each row with total sum of cells per cluster and mul by 100 to get percentages
        stats.iloc[:-1, :-1] = stats.iloc[:-1, :-1].div(stats['total_counts'][:-1], axis=0).mul(100).astype(int)

        # add column with percentage of all cells belonging to a cluster
        stats['perc_of_all_cells'] = np.around(stats['total_counts'] / stats['total_counts'].sum() * 100, decimals=1)

        # save cell mixture statistics to tissue
        self.tissue.uns['cell mixtures stats'] = stats.iloc[:, :]

    def plot_stats(self):
        stats = self.tissue.uns['cell mixtures stats']
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

    def plot_cluster_mixtures(self):
        # plot each cluster and its cells mixture
        sc.settings.set_figure_params(dpi=100, facecolor='white')
        stats = self.tissue.uns['cell mixtures stats']

        new_stats = stats.copy()
        new_stats = new_stats.drop(labels=['total_counts', 'perc_of_all_cells'], axis=1)
        new_stats = new_stats.drop(labels='total_cells', axis=0)
        for cluster in new_stats.iterrows():
            # only display clusters with more than min_cells_in_cluster cells
            if stats.loc[cluster[0]]['total_counts'] > self.min_cluster_size:
                # sort cell types by their abundnce in the cluster
                ct_perc = cluster[1].sort_values(ascending=False)
                # only cell types which have more than min_perc_to_show abundance will be shown
                ct_ind = [x for x in ct_perc.index[ct_perc>self.min_perc_to_show]]
                
                fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,6))
                fig.subplots_adjust(wspace=0.35)

                sc.pl.spatial(self.adata, groups=ct_ind, color=self.annotation, spot_size=self.spot_size, ax=ax[0], show=False, frameon=False)
                ax[0].set_title(f'cell types')
                ax[0].legend([f'{ind.get_text()} ({ct_perc[ind.get_text()]}%)' for ind in ax[0].get_legend().texts[:-1]], bbox_to_anchor=(1.0, 0.5), loc='center left', frameon=False, fontsize=12)
                sc.pl.spatial(self.adata, groups=[cluster[0]], color=f'tissue_{self.method_key}', spot_size=self.spot_size, ax=ax[1], show=False, frameon=False)
                ax[1].set_title(f'cell community')
                ax[1].legend([f'{ind.get_text()} ({stats.loc[ind.get_text(), "perc_of_all_cells"]}%)' for ind in ax[1].get_legend().texts[:-1]], bbox_to_anchor=(1.0, 0.5), loc='center left', frameon=False, fontsize=12)
                fig.savefig(os.path.join(self.dir_path, f'cmixtures_{self.params_suffix}_c{cluster[0]}.png'), bbox_inches='tight')

                plt.close()

    def boxplot_stats(self, stripplot=False):
        # box plot per cluster of cell type percentages distribution
        sc.settings.set_figure_params(dpi=100, facecolor='white')

        cluster_list = np.unique(self.tissue.obs['leiden'])
        
        for cluster in cluster_list:
            # for each window size a box plot is provided per cluster
            cl_win_cell_distrib = self.tissue[self.tissue.obs['leiden'] == cluster]
            for window_size in self.win_sizes_list:
                # extract only windows of specific size
                win_cell_distrib = cl_win_cell_distrib[cl_win_cell_distrib.obsm['spatial'][:,3] == window_size]
                # a DataFrame with cell percentages instead of normalized cel number is created
                win_cell_distrib_df = pd.DataFrame(win_cell_distrib.X / (self.total_cell_norm/100), columns=win_cell_distrib.var.index)

                # Reshape data into long format
                cell_type_distrib = pd.melt(win_cell_distrib_df, var_name='Cell Type', value_name='Percentage')
                cell_type_distrib = cell_type_distrib.sort_values(by='Cell Type')
                
                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,6))
                # plot boxplot of cell type percentages per mixture
                ax = sns.boxplot(x='Cell Type', y='Percentage', data=cell_type_distrib)
                if stripplot:
                    # overlap with a plot of specific percentage values. 
                    # Jitter allows dots to move left and right for better visibility of all points
                    ax = sns.stripplot(x='Cell Type', y='Percentage', data=cell_type_distrib, jitter=True, color='black', size=2)
                # remove top and right frame of the plot
                sns.despine(top=True, right=True)
                ax.set_title(f'Cell community {cluster} ({self.adata.uns["sample_name"]})')
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                ax.xaxis.tick_bottom() # x axis on the bottom
                fig.savefig(os.path.join(self.dir_path, f'boxplot_c{cluster}_ws{window_size}.png'), bbox_inches='tight')

                plt.close()

    def colorplot_stats(self, color_system='hsv'):
        supported_color_systems = ['hsv', 'rgb']
        if color_system in supported_color_systems:
            stats = self.tissue.uns['cell mixtures']
            # divide each row with total sum of cells per cluster and mul by 100 to get percentages
            stats.iloc[:, :] = stats.iloc[:, :].div(np.array([sum(stats.loc[row, :]) for row in stats.index]), axis=0).mul(100).astype(int)

            cx_min = int(np.min(self.adata.obsm['spatial'][:,0]))
            cy_min = int(np.min(self.adata.obsm['spatial'][:,1]))
            cx_max = int(np.max(self.adata.obsm['spatial'][:,0]))
            cy_max = int(np.max(self.adata.obsm['spatial'][:,1]))

            for cluster in stats.iterrows():
                ct_perc = cluster[1].sort_values(ascending=False)
                top_three_ct = ct_perc.index.values[0:3]

                cl_win_cell_distrib = self.tissue[self.tissue.obs['leiden'] == cluster[0]]
                cl_win_cell_distrib = cl_win_cell_distrib[:, top_three_ct]

                # for each pair of window size and sliding step a separate color plot should be made
                for window_size, sliding_step in zip(self.win_sizes_list, self.sliding_steps_list):
                    # extract data for a specfic window size
                    win_cell_distrib = cl_win_cell_distrib[cl_win_cell_distrib.obsm['spatial'][:,3] == window_size]
                    # data is a DataFrame with rows for each window, columns of 3 top most cell types for
                    # current cell mixture cluster, with data on cell type in percentages [0.00-1.00]
                    # The last is achieved by dividing the features with self.total_cell_norm
                    data_df = pd.DataFrame(win_cell_distrib.X/self.total_cell_norm, columns=win_cell_distrib.var.index, index=win_cell_distrib.obs.index)
                    # init image
                    mixture_image = np.zeros(shape=(cy_max-cy_min+1, cx_max-cx_min+1, 3), dtype=np.float32)

                    for window in data_df.iterrows():
                        wx = int(window[0].split("_")[0])
                        wy = int(window[0].split("_")[1])
                        mixture_image[int(wy*sliding_step-cy_min) : int(wy*sliding_step + window_size-cy_min), int(wx*sliding_step-cx_min) : int(wx*sliding_step + window_size-cx_min), :] = 1 - window[1].values.astype(np.float32)
                    
                    # convert image of selected color representation to rgb
                    if color_system == 'hsv':
                        # if hsv display the 1 - percentage since the colors will be too dark
                        rgb_image = color.hsv2rgb(mixture_image)
                    elif color_system == 'rgb':
                        rgb_image = mixture_image
                    # plot the colored window image of the cell scatterplot
                    fig, ax = plt.subplots(nrows=1, ncols=1)
                    # cell scatterplot for visual spatial reference
                    plt.scatter(x = self.adata.obsm['spatial'][:,0]-cx_min, y=self.adata.obsm['spatial'][:,1]-cy_min, c='#CCCCCC', marker='.', s=0.5, zorder=1)
                    # mask of window positions
                    window_mask = rgb_image[:,:,1] != 0
                    # mask adjusted to alpha channel and added to rgb image
                    window_alpha = (window_mask==True).astype(int)[..., np.newaxis]
                    rgba_image = np.concatenate([rgb_image, window_alpha], axis=2)
                    # plot windows, where empty areas will have alpha=0, making them transparent
                    plt.imshow(rgba_image, zorder=2)
                    plt.axis('off')
                    ax.grid(visible=False)
                    ax.set_title(f'{color_system} of community {cluster[0]} win size {window_size}, step {sliding_step} - top 3 cell types\n({self.adata.uns["sample_name"]})')
                    
                    if color_system == 'hsv':
                        plane_names = ['H\'', 'S\'', 'V\'']
                    elif color_system == 'rgb':
                        plane_names = ['R\'', 'G\'', 'B\'']
                    
                    ax.text(1.05, 0.5, f'{plane_names[0]} - {top_three_ct[0]} ({ct_perc[top_three_ct[0]]}%)\n{plane_names[1]} - {top_three_ct[1]} ({ct_perc[top_three_ct[1]]}%)\n{plane_names[2]} - {top_three_ct[2]} ({ct_perc[top_three_ct[2]]}%)', \
                                transform=ax.transAxes, fontsize=12, va='center', ha='left')
            
                    fig.savefig(os.path.join(self.dir_path, f'colorplot_{color_system}_c{cluster[0]}_ws{window_size}_ss{sliding_step}.png'), bbox_inches='tight', dpi=200)

                    plt.close()
        else:
            logging.warn(f'Unsupported color system: {color_system}.')

    def plot_celltype_table(self):
        sc.settings.set_figure_params(dpi=300, facecolor='white')
        sns.set(font_scale=0.5)

        stats = self.tissue.uns['cell mixtures']

        # calculate percentage of all cells of one cell type belonging to each cluster
        ct_perc_per_celltype = stats.iloc[:,:].div(np.array([sum(stats.loc[:, col]) for col in stats.columns]), axis=1).mul(100).astype(int)
        ct_perc_per_cluster = stats.iloc[:,:].div(np.array([sum(stats.loc[row, :]) for row in stats.index]), axis=0).mul(100).astype(int)

        # divide cell numbers with total number of cells per cluster to obtain ct perc per cluster
        for cluster in stats.iterrows():
            ct_perc_per_cluster_sorted = ct_perc_per_cluster.loc[cluster[0], :].sort_values(ascending=False)
            if (ct_perc_per_cluster_sorted[self.min_num_celltype-1] < self.min_perc_celltype) or (sum(stats.loc[cluster[0], :]) < self.min_cluster_size):
                # remove all clusters that have low number of cells, or high abundance of single cell type
                stats = stats.drop(labels=cluster[0], axis=0)
                ct_perc_per_celltype = ct_perc_per_celltype.drop(labels=cluster[0], axis=0)
                ct_perc_per_cluster = ct_perc_per_cluster.drop(labels=cluster[0], axis=0)

        # remove cell types that are not a significant part of any heterogeneous cluster
        for celltype in stats.columns:
            if (max(ct_perc_per_celltype.loc[:, celltype]) < self.min_perc_to_show):
                stats = stats.drop(labels=celltype, axis=1)
                ct_perc_per_celltype = ct_perc_per_celltype.drop(labels=celltype, axis=1)
                ct_perc_per_cluster = ct_perc_per_cluster.drop(labels=celltype, axis=1)

        ncols = len(stats)
        # table will have a clumn for each cluster and first column for cell types
        fig, axes = plt.subplots(nrows=1, ncols=ncols+1)
        # no space between columns
        fig.subplots_adjust(wspace=0, hspace=0)

        # create a list of colors to use as the row background
        scplspatial_colors = matplotlib.rcParams["axes.prop_cycle"].by_key()['color']

        # create a dictionary mapping each cluster to its corresponding color
        cluster_color = dict(zip(stats.columns, [scplspatial_colors[int(x)] for x in stats.index.values]))

        # cell type colors from adata.uns['annotation_colors']
        cell_type_colors = [self.adata.uns[f'{self.annotation}_colors'][np.unique(self.adata.obs[self.annotation])==ct][0] for ct in stats.columns]
        cmap = mcolors.ListedColormap([mcolors.hex2color(hexc) for hexc in cell_type_colors])

        for i, ax in enumerate(axes):
            if i == 0:
                g = sns.heatmap(np.array(range(len(stats.columns)))[:,np.newaxis], linewidths=0.5, linecolor='gray', \
                                annot=np.array([column for column in stats.columns])[:, np.newaxis], ax=ax, cbar=False, \
                                      cmap=cmap, fmt="", xticklabels=False, yticklabels=False, square=None)        
            else:
                table_annotation = np.array([f'{ct_perc_per_cluster.iloc[i-1, int(x)]}%\n({ct_perc_per_celltype.iloc[i-1, int(x)]}%)' for x in range(len(stats.columns))])[:, np.newaxis]
                g = sns.heatmap(np.array(stats.iloc[i-1, :])[:, np.newaxis], linewidths=0.5, linecolor='gray', annot=table_annotation, cbar=False, cmap=[cluster_color[stats.columns[i-1]]], ax=ax, fmt='', xticklabels=True, yticklabels=False, square=None)
                g.set_xticklabels([f'cluster {stats.index[i-1]}'], rotation=0)
                g.xaxis.tick_top() # x axis on top
        axes[i//2].set_title('Cell type abundance per cluster (and per cel type set)')
        fig.savefig(os.path.join(self.dir_path, f'celltype_table_{self.params_suffix}.png'), bbox_inches='tight')

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