import logging
import os
from collections import Counter
from collections import defaultdict
from tqdm.auto import tqdm

import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
import anndata as ad

from .utils import timeit
from core import CommunityClusteringAlgo


class SlidingWindow(CommunityClusteringAlgo):
    def __init__(self, adata, slice_id, input_file_path, **params):
        super().__init__(adata, slice_id, input_file_path,  **params)
        self.params_suffix = f"_sldwin_sl{self.slice_id}_r{self.resolution}_ws{self.win_sizes}_ss{self.sliding_steps}_en{self.entropy_thres}_sct{self.scatter_thres}_dwr{self.downsample_rate}_mcc{self.min_cells_coeff}"
        self.filename = self.adata.uns['sample_name'] + self.params_suffix
        self.dir_path = os.path.join(self.adata.uns['algo_params']['out_path'], self.filename)
        # create results folder
        if (not os.path.exists(self.dir_path)):
            os.mkdir(self.dir_path)

        self.method_key = 'sliding_window'
    
    @timeit
    def run(self):
        if self.tfile==None:
            self.calc_feature_matrix(self.win_sizes_list[0], self.sliding_steps_list[0])
        else:
            if self.tfile.endswith('.h5ad'):
                self.tissue = sc.read(self.tfile)
            else:
                raise AttributeError(f"File '{self.tfile}' extension is not .h5ad")


    def calc_feature_matrix(self, win_size, sliding_step):
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


        self.adata.obs['window_spatial'] = self.adata.obs['Centroid_X'].astype(str) +'_'+self.adata.obs['Centroid_Y'].astype(str) + '_' + str(self.slice_id) + '_' + str(win_size)
        self.adata.obs[f'window_spatial_{win_size}'] = self.adata.obs['Centroid_X'].astype(str) +'_'+self.adata.obs['Centroid_Y'].astype(str) + '_' + str(self.slice_id) + '_' + str(win_size)
        
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
            w_curr = int(subwindow.split("_")[3])

            for slide_x in range(0, np.min([bin_slide_ratio, x_max-x_curr+1])):
                for slide_y in range(0, np.min([bin_slide_ratio, y_max-y_curr+1])):  # starts from 1 since values with coordinates (0,0) are already written by initializing with ret[subwindow]
                    if (f'{x_curr + slide_x}_{y_curr + slide_y}_{z_curr}_{w_curr}') in ret.keys():
                        feature_matrix[subwindow] = {k: 
                                                    feature_matrix[subwindow].get(k, 0) + ret[f'{x_curr + slide_x}_{y_curr + slide_y}_{z_curr}_{w_curr}'].get(k, 0)
                                                    for k in set(feature_matrix[subwindow]).union(ret[f'{x_curr + slide_x}_{y_curr + slide_y}_{z_curr}_{w_curr}'])}

        feature_matrix = pd.DataFrame(feature_matrix).T
        # feature_matrix is placd in AnnData object with specified spatial cooridnated of the sliding windows
        self.tissue = AnnData(feature_matrix.astype(np.float32), dtype=np.float32)
        # spatial coordinates are expanded with 3rd dimension with slice_id 
        # this should enable calculation of multislice cell communities
        self.tissue.obsm['spatial'] = np.array([x.split('_') for x in feature_matrix.index]).astype(int)
        # self.tissue.obsm['window_size'] = np.array([win_size for _ in feature_matrix.index])
        self.tissue.obs['window_size'] = np.array([win_size for _ in feature_matrix.index])
        self.tissue.obs = self.tissue.obs.copy()
        self.tissue.obs['window_cell_sum'] = np.sum(self.tissue.X, axis=1)
        # scale the feature vector by the total numer of cells in it
        self.tissue.X = ((self.tissue.X.T * self.total_cell_norm) / self.tissue.obs['window_cell_sum'].values).T
        # remove feature vectors which have less than a specified amount of cells
        mean_cell_sum = np.mean(self.tissue.obs['window_cell_sum'].values)
        stddev_cell_sum = np.std(self.tissue.obs['window_cell_sum'].values)
        min_cells_per_window = mean_cell_sum - self.min_cells_coeff * stddev_cell_sum
        self.tissue = self.tissue[self.tissue.obs['window_cell_sum'].values >= min_cells_per_window, :]

    ## COMMUNITY CALLING
    # Define subwindow cluster label based on labels of all overlapping windows
    # Functions goes through all subwindow positions and gathers clustering
    # labels of all windows that contain it. Final label of the subwindow is
    # the one with majority vote. Window cluster labels are in self.tissue_pruned.obs['leiden']
    # and the subwindow labels are placed in self.tissue.obs['leiden_max_vote']
    # self.tissue_pruned is used only for clustering and is discarded
    def community_calling(self):
        win_size = self.win_sizes_list[0]
        sliding_step = self.sliding_steps_list[0]

        bin_slide_ratio = int(win_size/sliding_step)
        x_min = self.adata.obs['Centroid_X'].min()
        y_min = self.adata.obs['Centroid_Y'].min()

        # max voting on cluster labels
        subwindow_locations = np.unique(self.adata.obs['window_spatial'])
        # variable for final subwindow labels
        leiden_max_vote = pd.Series(index=subwindow_locations, name='leiden_max_vote', dtype=np.float64)
        for location in subwindow_locations:
            # extract x,y,z coordinates from location string
            x_curr, y_curr, z_curr, _ = np.array(location.split("_")).astype(int)
            # index of subwindow is in the top left corner of the whole window
            subwindow_labels = {}
            for slide_x in range(0, np.min([bin_slide_ratio, x_curr - x_min + 1])):
                for slide_y in range(0, np.min([bin_slide_ratio, y_curr - y_min + 1])):
                    # check if location exist (spatial area is not complete)
                    if (f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}_{win_size}') in self.tissue.obs.index:
                        new_value = self.tissue.obs.loc[f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}_{win_size}', 'leiden']
                        subwindow_labels[new_value] = subwindow_labels[new_value] + 1 if new_value in subwindow_labels.keys() else 1
            
            # MAX VOTE
            # max vote is saved in a new variable (not changed in tissue.obs) so that it does not have diagonal effect on other labels during refinement
            # max_voting result is created for each subwindow, while the 'leiden' clustering was defined for each window
            leiden_max_vote.loc[location] = max(subwindow_labels, key=subwindow_labels.get) if subwindow_labels!={} else np.nan

        # copy clustering results from subwindows to cells of those subwindows in adata object
        self.adata.obs.loc[:, f'tissue_{self.method_key}'] = list(leiden_max_vote.loc[self.adata.obs['window_spatial']])

        logging.info(f'Sliding window cell mixture calculation done. Added results to adata.obs["tissue_{self.method_key}"]')


class SlidingWindowMultipleSizes(SlidingWindow):
    def __init__(self, adata, slice_id, input_file_path, **params):
        super().__init__(adata, slice_id, input_file_path,  **params)

    def run(self):
        tissue_list = []

        if self.tfile==None:
            n = len(self.win_sizes_list)
            for i in range(n):
                super().calc_feature_matrix(self.win_sizes_list[i], self.sliding_steps_list[i])
                tissue_list.append(self.tissue)
            
            self.tissue = ad.concat(tissue_list, axis=0, join='outer')
        else:
            if self.tfile.endswith('.h5ad'):
                self.tissue = sc.read(self.tfile)
            else:
                raise AttributeError(f"File '{self.tfile}' extension is not .h5ad")
    
    def community_calling(self):
        self.win_sizes_list = [int(w) for w in self.win_sizes.split(',')]
        self.sliding_steps_list = [int(s) for s in self.sliding_steps.split(',')]
        if len(self.win_sizes_list) == 1 and len(self.sliding_steps_list) == 1:
            super().community_calling()
        else:
            self.community_calling_multiple_window_sizes_per_cell()
    
    # def community_calling_multiple_window_sizes_per_window(self):
    #     x_min = self.adata.obs['Centroid_X'].min()
    #     y_min = self.adata.obs['Centroid_Y'].min()

    #     for win_size, sliding_step in zip(self.win_sizes_list, self.sliding_steps_list):
    #         bin_slide_ratio = int(win_size/sliding_step)

    #         # max voting on cluster labels, take only current window size
    #         subwindow_locations = np.unique(self.adata.obs[f'window_spatial_{str(win_size)}'])
    #         # variable for final subwindow labels
    #         leiden_max_vote = pd.Series(index=subwindow_locations, name='leiden_max_vote', dtype=np.float64)
    #         for location in subwindow_locations:
    #             # extract x,y,z coordinates from location string
    #             x_curr, y_curr, z_curr, w_size = np.array(location.split("_")).astype(int)
    #             # index of subwindow is in the top left corner of the whole window
    #             subwindow_labels = {}
    #             for slide_x in range(0, np.min([bin_slide_ratio, x_curr - x_min + 1])):
    #                 for slide_y in range(0, np.min([bin_slide_ratio, y_curr - y_min + 1])):
    #                     # check if location exist (spatial area is not complete)
    #                     if (f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}_{w_size}') in self.tissue.obs.index:
    #                         new_value = self.tissue.obs.loc[f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}_{w_size}', 'leiden']
    #                         subwindow_labels[new_value] = subwindow_labels[new_value] + 1 if new_value in subwindow_labels.keys() else 1
                
    #             # MAX VOTE PER WINDOW SIZE
    #             # max vote is saved in a new variable (not changed in tissue.obs) so that it does not have diagonal effect on other labels during refinement
    #             # max_voting result is created for each subwindow, while the 'leiden' clustering was defined for each window
    #             leiden_max_vote.loc[location] = max(subwindow_labels, key=subwindow_labels.get) if subwindow_labels!={} else np.nan

    #         # copy clustering results from subwindows to cells of those subwindows in adata object
    #         self.adata.obs.loc[:, f'tissue_{self.method_key}_{str(win_size)}'] = list(leiden_max_vote.loc[self.adata.obs[f'window_spatial_{str(win_size)}']])

    #         logging.info(f'Sliding window cell mixture calculation done for window size = {str(win_size)}. Added results to adata.obs["tissue_{self.method_key}_{str(win_size)}"]')
        
        # perform final max - label with most appearances in all window sizes will be chosen
        # leiden_final_max_vote = pd.Series(index=self.adata.obs.index, name='leiden_final_max_vote', dtype=np.int32)
        # cols = [col for col in self.adata.obs.columns if f"tissue_{self.method_key}" in col]
        # df_voting = self.adata.obs.loc[:, [col for col in self.adata.obs.columns if f"tissue_{self.method_key}" in col]]
        # df_voting[f'tissue_{self.method_key}'] = df_voting.apply(lambda row : Counter(list(row)).most_common(1)[0][0], axis=1)
        # self.adata.obs.loc[:, f'tissue_{self.method_key}'] = df_voting[f'tissue_{self.method_key}']
        # self.adata.obs = self.adata.obs.drop(columns=cols, axis=1)

    #     logging.info(f'Sliding window cell mixture calculation done. Added results to adata.obs["tissue_{self.method_key}"]')

    @timeit
    def community_calling_multiple_window_sizes_per_cell(self):
        x_min = self.adata.obs['Centroid_X'].min()
        y_min = self.adata.obs['Centroid_Y'].min()
        self.adata.obs[f'tissue_{self.method_key}'] = pd.Series(index=self.adata.obs.index, dtype=str)
        cache = {}

        for index, cell in tqdm(self.adata.obs.iterrows(), desc="Per cell computation... ", total=self.adata.obs.shape[0]):
            cell_labels_all = defaultdict(int)

            for win_size, sliding_step in zip(self.win_sizes_list, self.sliding_steps_list):
                sliding_step = (win_size/int((win_size/sliding_step))) if sliding_step!=None else win_size
                bin_slide_ratio = int(win_size/sliding_step)
                cell_labels = defaultdict(int)
                x_curr, y_curr, z_curr, w_size = [int(num) for num in cell[f'window_spatial_{win_size}'].split("_")]

                if (x_curr, y_curr, z_curr, w_size) in cache:
                    cell_labels = cache[(x_curr, y_curr, z_curr, w_size)]
                    cell_labels_all = {
                        key: cell_labels_all.get(key, 0) + cell_labels.get(key, 0)
                        for key in set(cell_labels_all) | set(cell_labels)
                    }
                    continue
                # index of cell is in the top left corner of the whole window
                for slide_x in range(0, np.min([bin_slide_ratio, x_curr - x_min + 1])):
                    for slide_y in range(0, np.min([bin_slide_ratio, y_curr - y_min + 1])):
                        # check if location exist (spatial area is not complete)
                        if (f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}_{win_size}') in self.tissue.obs.index:
                            win_label = self.tissue.obs.loc[f'{x_curr - slide_x}_{y_curr - slide_y}_{z_curr}_{win_size}', 'leiden']
                            cell_labels[win_label] += 1

                cache[(x_curr, y_curr, z_curr, w_size)] = cell_labels
                cell_labels_all = {
                    key: cell_labels_all.get(key, 0) + cell_labels.get(key, 0)
                    for key in set(cell_labels_all) | set(cell_labels)
                }
            max_vote_label = max(cell_labels_all, key=cell_labels_all.get) if cell_labels_all != {} else np.nan
            self.adata.obs.loc[index, f'tissue_{self.method_key}'] = max_vote_label
        
        self.adata.obs[f'tissue_{self.method_key}'] = pd.Categorical(self.adata.obs[f'tissue_{self.method_key}'].fillna('unknown'))
        
