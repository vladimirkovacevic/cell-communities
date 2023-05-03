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
        self.params_suffix = f"_sldwin_sl{self.slice_id}_r{self.resolution}_ws{self.win_sizes}_ss{self.sliding_steps}_en{self.entropy_thres}_sct{self.scatter_thres}_dwr{self.downsample_rate}_tn{self.total_cell_norm}_mcc{self.min_cells_coeff}"
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
        whole_feature_matrix = pd.DataFrame()

        for win_size, sliding_step in zip(self.win_sizes_list, self.sliding_steps_list):
            # window size needs to be a multiple of sliding step
            sliding_step = (win_size/int((win_size/sliding_step))) if sliding_step!=None else win_size
            bin_slide_ratio = int(win_size/sliding_step)

            # create centroids for each sliding step of windows
            # .obs data is assigned as pd.Series since sometimes the new column added to .obs Dataframe can have 'nan' values
            # if the index of data doesn't match the index of .obs
            self.adata.obs['Centroid_X'] = pd.Series(((self.adata.obsm['spatial'][:,0])/sliding_step).astype(int), index=self.adata.obs_names)
            self.adata.obs['Centroid_Y'] = pd.Series(((self.adata.obsm['spatial'][:,1])/sliding_step).astype(int), index=self.adata.obs_names)
            # need to understand borders and padding
            # subwindows belonging to borders will not have a complete cell count
            x_max = self.adata.obs['Centroid_X'].max()
            y_max = self.adata.obs['Centroid_Y'].max()

            self.adata.obs['x_y'] = self.adata.obs['Centroid_X'].astype(str) +'_'+self.adata.obs['Centroid_Y'].astype(str)
            
            tmp = self.adata.obs[['x_y', self.annotation]]
            ret = {}
            # calculate features for each subwindow
            for sw_ind, sw_data in tmp.groupby('x_y'):
                templete_dic = {ct:0 for ct in self.unique_cell_type}
                for cell in sw_data[self.annotation]:
                    templete_dic[cell]+=1
                ret[sw_ind] = templete_dic
            # merge features by windows
            feature_matrix = {}
            for subwindow in ret.keys():
                # index of window is in the top left corner of the whole window
                subwindow_id = str(subwindow) + '_' + str(win_size) 
                feature_matrix[subwindow_id] = {}
                x_curr = int(subwindow.split("_")[0])
                y_curr = int(subwindow.split("_")[1])
                # # count the number of subwindows participating in full window feature vec
                # num_subw = 0

                for slide_x in range(0, np.min([bin_slide_ratio, x_max-x_curr+1])):
                    for slide_y in range(0, np.min([bin_slide_ratio, y_max-y_curr+1])):  # starts from 1 since values with coordinates (0,0) are already written by initializing with ret[subwindow]
                        if (f'{x_curr + slide_x}_{y_curr + slide_y}') in ret.keys():
                            # num_subw += 1
                            feature_matrix[subwindow_id] = {k: 
                                                        feature_matrix[subwindow_id].get(k, 0) + ret[f'{x_curr + slide_x}_{y_curr + slide_y}'].get(k, 0)
                                                        for k in set(feature_matrix[subwindow_id]).union(ret[f'{x_curr + slide_x}_{y_curr + slide_y}'])}

                    
            feature_matrix = pd.DataFrame(feature_matrix).T
            pd.concat([whole_feature_matrix, feature_matrix], axis = 0)
        

        # feature_matrix is placd in AnnData object with specified spatial cooridnated of the sliding windows
        self.tissue = AnnData(whole_feature_matrix.astype(np.float32), dtype=np.float32)
        # spatial coordinates are expanded with 3rd dimension with slice_id 
        # this should enable calculation of multislice cell communities
        self.tissue.obsm['spatial'] = np.array([[x.split('_')[0], x.split('_')[1], self.slice_id] for x in whole_feature_matrix.index]).astype(int)
        self.tissue.obsm['window_size'] = np.array([x.split('_')[2] for x in whole_feature_matrix.index]).astype(int)
        
        self.tissue.obs = self.tissue.obs.copy()
        self.tissue.obs['window_cell_sum'] = np.sum(self.tissue.X, axis=1)
        # remove feature vectors which have less than a specified amount of time
        mean_cell_sum = np.mean(self.tissue.obs['window_cell_sum'].values)
        stddev_cell_sum = np.std(self.tissue.obs['window_cell_sum'].values)
        min_cells_per_window = mean_cell_sum - self.min_cells_coeff * stddev_cell_sum
        self.tissue = self.tissue[self.tissue.obs['window_cell_sum'].values >= min_cells_per_window, :]
        # scale the feature vector by the total numer of cells in it
        self.tissue.X = ((self.tissue.X.T * self.total_cell_norm) / self.tissue.obs['window_cell_sum'].values).T

    def community_calling(self):
        win_size = self.win_sizes_list[0]
        sliding_step = self.sliding_steps_list[0]

        bin_slide_ratio = int(win_size/sliding_step)
        
        x_min = self.adata.obs['Centroid_X'].min()
        y_min = self.adata.obs['Centroid_Y'].min()
        # max voting on cluster labels
        # init the new obs column
        self.tissue.obs['leiden_max_vote'] = list('x' for x in range(len(self.tissue.obs.index)))
        for x_curr, y_curr, _ in self.tissue.obsm['spatial']:
            # index of subwindow is in the top left corner of the whole window
            subwindow_labels = {}
            for slide_x in range(0, np.min([bin_slide_ratio, x_curr - x_min + 1])):
                for slide_y in range(0, np.min([bin_slide_ratio, y_curr - y_min + 1])):
                    # check if location exist (spatial area is not complete)
                    if (f'{x_curr - slide_x}_{y_curr - slide_y}') in self.tissue.obs.index:
                        new_value = self.tissue.obs.loc[f'{x_curr - slide_x}_{y_curr - slide_y}', 'leiden']
                        subwindow_labels[new_value] = subwindow_labels[new_value] + 1 if new_value in subwindow_labels.keys() else 1
            
            # max vote
            # max vote should be saved in a new obs column so that it does not have diagonal effect on
            # other labels during refinement
            # max_voting result is created for each subwindow, while the 'leiden' clustering was defined for each window
            self.tissue.obs.loc[f'{x_curr}_{y_curr}', 'leiden_max_vote'] = max(subwindow_labels, key=subwindow_labels.get)

        # since tissue does not have all the indexes, for the ones without enough cells are removed we need to initialize
        # self.adata.obs for method results
        self.adata.obs[f'tissue_{self.method_key}'] = np.nan
        idx_mask = self.adata.obs['x_y'].isin(self.tissue.obs.index)
        self.adata.obs.loc[idx_mask, f'tissue_{self.method_key}'] = list(self.tissue.obs.loc[self.adata.obs.loc[idx_mask, 'x_y'], 'leiden_max_vote'])

        logging.info(f'Sliding window cell mixture calculation done. Added results to adata.obs["tissue_{self.method_key}"]')


class SlidingWindowMultipleSizes(SlidingWindow):
    def __init__(self, adata, slice_id, input_file_path, **params):
        super().__init__(adata, slice_id, input_file_path,  **params)
    
    def community_calling(self):
        if len(self.win_sizes_list)==1 and len(self.sliding_steps_list)==1:
            super().community_calling()
        else: # TODO
            pass
