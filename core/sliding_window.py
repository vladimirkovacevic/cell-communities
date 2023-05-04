import logging
import os

import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData

from .utils import timeit
from core import CommunityClusteringAlgo


class SlidingWindow(CommunityClusteringAlgo):
    def __init__(self, adata, slice_id, input_file_path, **params):
        super().__init__(adata, slice_id, input_file_path,  **params)
        self.params_suffix = f"_sldwin_sl{self.slice_id}_r{self.resolution}_ws{self.win_size}_ss{self.sliding_step}_en{self.entropy_thres}_sct{self.scatter_thres}_dwr{self.downsample_rate}_mcc{self.min_cells_coeff}"
        self.filename = self.adata.uns['sample_name'] + self.params_suffix
        self.dir_path = os.path.join(self.adata.uns['algo_params']['out_path'], self.filename)
        # create results folder
        if (not os.path.exists(self.dir_path)):
            os.mkdir(self.dir_path)

        self.method_key = 'sliding_window'
    
    @timeit
    def run(self):
        if self.tfile==None:
            self.calc_feature_matrix()
        else:
            if self.tfile.endswith('.h5ad'):
                self.tissue = sc.read(self.tfile)
            else:
                raise AttributeError(f"File '{self.tfile}' extension is not .h5ad")


    def calc_feature_matrix(self):
        # window size needs to be a multiple of sliding step
        self.sliding_step = (self.win_size/int((self.win_size/self.sliding_step))) if self.sliding_step!=None else self.win_size
        bin_slide_ratio = int(self.win_size/self.sliding_step)

        # create centroids for each sliding step of windows
        # .obs data is assigned as pd.Series since sometimes the new column added to .obs Dataframe can have 'nan' values
        # if the index of data doesn't match the index of .obs
        self.adata.obs['Centroid_X'] = pd.Series(((self.adata.obsm['spatial'][:,0])/self.sliding_step).astype(int), index=self.adata.obs_names)
        self.adata.obs['Centroid_Y'] = pd.Series(((self.adata.obsm['spatial'][:,1])/self.sliding_step).astype(int), index=self.adata.obs_names)
        # need to understand borders and padding
        # subwindows belonging to borders will not have a complete cell count
        x_max = self.adata.obs['Centroid_X'].max()
        y_max = self.adata.obs['Centroid_Y'].max()


        self.adata.obs['window_spatial'] = self.adata.obs['Centroid_X'].astype(str) +'_'+self.adata.obs['Centroid_Y'].astype(str) + '_' + str(self.slice_id)
        
        tmp = self.adata.obs[['window_spatial', self.annotation]]
        ret = {}
        # calculate features for each subwindow
        for sw_ind, sw_data in tmp.groupby('window_spatial'):
            templete_dic = {ct:0 for ct in self.unique_cell_type}
            for cell in sw_data[self.annotation]:
                templete_dic[cell]+=1
            ret[sw_ind] = templete_dic
        # merge features by windows
        feature_matrix = {}
        for subwindow in ret.keys():
            # index of window is in the top left corner of the whole window
            feature_matrix[subwindow] = {}
            x_curr = int(subwindow.split("_")[0])
            y_curr = int(subwindow.split("_")[1])
            z_curr = int(subwindow.split("_")[2])

            for slide_x in range(0, np.min([bin_slide_ratio, x_max-x_curr+1])):
                for slide_y in range(0, np.min([bin_slide_ratio, y_max-y_curr+1])):  # starts from 1 since values with coordinates (0,0) are already written by initializing with ret[subwindow]
                    if (f'{x_curr + slide_x}_{y_curr + slide_y}_{z_curr}') in ret.keys():
                        feature_matrix[subwindow] = {k: 
                                                    feature_matrix[subwindow].get(k, 0) + ret[f'{x_curr + slide_x}_{y_curr + slide_y}_{z_curr}'].get(k, 0)
                                                    for k in set(feature_matrix[subwindow]).union(ret[f'{x_curr + slide_x}_{y_curr + slide_y}_{z_curr}'])}

                
        feature_matrix = pd.DataFrame(feature_matrix).T
        # feature_matrix is placd in AnnData object with specified spatial cooridnated of the sliding windows
        self.tissue = AnnData(feature_matrix.astype(np.float32), dtype=np.float32)
        # spatial coordinates are expanded with 3rd dimension with slice_id 
        # this should enable calculation of multislice cell communities
        self.tissue.obsm['spatial'] = np.array([x.split('_') for x in feature_matrix.index]).astype(int)
        self.tissue.obs = self.tissue.obs.copy()
        self.tissue.obs['window_cell_sum'] = np.sum(self.tissue.X, axis=1)
        # scale the feature vector by the total numer of cells in it
        self.tissue.X = ((self.tissue.X.T * self.total_cell_norm) / self.tissue.obs['window_cell_sum'].values).T
        # remove feature vectors which have less than a specified amount of cells
        mean_cell_sum = np.mean(self.tissue.obs['window_cell_sum'].values)
        stddev_cell_sum = np.std(self.tissue.obs['window_cell_sum'].values)
        min_cells_per_window = mean_cell_sum - self.min_cells_coeff * stddev_cell_sum
        self.tissue_pruned = self.tissue[self.tissue.obs['window_cell_sum'].values >= min_cells_per_window, :].copy()

    ## COMMUNITY CALLING
    # Define subwindow cluster label based on labels of all overlapping windows
    # Functions goes through all subwindow positions and gathers clustering
    # labels of all windows that contain it. Final label of the subwindow is
    # the one with majority vote. Window cluster labels are in self.tissue_pruned.obs['leiden']
    # and the subwindow labels are placed in self.tissue.obs['leiden_max_vote']
    # self.tissue_pruned is used only for clustering and is discarded
    def community_calling(self):
        bin_slide_ratio = int(self.win_size/self.sliding_step)
        x_min = self.adata.obs['Centroid_X'].min()
        y_min = self.adata.obs['Centroid_Y'].min()
        # max voting on cluster labels
        # init the new obs column
        self.tissue.obs = self.tissue.obs.copy()
        self.tissue.obs['leiden_max_vote'] = np.nan
        for x_curr, y_curr, z_curr in self.tissue.obsm['spatial']:
            # index of subwindow is in the top left corner of the whole window
            subwindow_labels = {}
            for slide_x in range(0, np.min([bin_slide_ratio, x_curr - x_min + 1])):
                for slide_y in range(0, np.min([bin_slide_ratio, y_curr - y_min + 1])):
                    # check if location exist (spatial area is not complete)
                    if (f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}') in self.tissue_pruned.obs.index:
                        new_value = self.tissue_pruned.obs.loc[f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}', 'leiden']
                        subwindow_labels[new_value] = subwindow_labels[new_value] + 1 if new_value in subwindow_labels.keys() else 1
            
            # MAX VOTE
            # max vote should be saved in a new obs column so that it does not have diagonal effect on
            # other labels during refinement
            # max_voting result is created for each subwindow, while the 'leiden' clustering was defined for each window
            self.tissue.obs.loc[f'{x_curr}_{y_curr}_{z_curr}', 'leiden_max_vote'] = max(subwindow_labels, key=subwindow_labels.get) if subwindow_labels!={} else np.nan

        # copy clustering results from subwindows to cells of those subwindows in adata object
        self.adata.obs.loc[:, f'tissue_{self.method_key}'] = list(self.tissue.obs.loc[self.adata.obs['window_spatial'], 'leiden_max_vote'])

        logging.info(f'Sliding window cell mixture calculation done. Added results to adata.obs["tissue_{self.method_key}"]')

        # delete helper anndata object
        del self.tissue_pruned
