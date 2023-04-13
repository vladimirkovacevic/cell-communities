import logging
import os

import scanpy as sc
import numpy as np
import pandas as pd

from anndata import AnnData
from matplotlib import pyplot as plt
import skimage.measure
import scipy.ndimage.measurements

from .utils import timeit
from core import CommunityClusteringAlgo

class SlidingWindow(CommunityClusteringAlgo):
    def __init__(self, adata, **params):
        super().__init__(adata, **params)
        self.params_suffix = f"_sldwin_r{self.resolution}_ws{self.win_size}_ss{self.sliding_step}_entt{self.entropy_thres}_scatt{self.scatter_thres}_dwnsr{self.downsample_rate}"
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
            self.calculate_spatial_cell_type_metrics()

            var_use = self.tissue.var.loc[(self.tissue.var['entropy']<self.entropy_thres) & (self.tissue.var['scatteredness']<self.scatter_thres)].index
            self.tissue.raw = self.tissue
            self.tissue = self.tissue[:, var_use]
        else:
            if self.tfile.endswith('.h5ad'):
                self.tissue = sc.read(self.tfile)
            else:
                raise AttributeError(f"File '{self.tfile}' extension is not .h5ad")
        if f'tissue_{self.method_key}' not in self.adata.obs.keys():
            self.cluster()
        
        self.plot_clustering()
        self.calculate_cell_mixture_stats()
        self.plot_stats()
        self.save_results()



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
            feature_matrix[subwindow] = {}
            x_curr = int(subwindow.split("_")[0])
            y_curr = int(subwindow.split("_")[1])
            # # count the number of subwindows participating in full window feature vec
            # num_subw = 0

            for slide_x in range(0, np.min([bin_slide_ratio, x_max-x_curr+1])):
                for slide_y in range(0, np.min([bin_slide_ratio, y_max-y_curr+1])):  # starts from 1 since values with coordinates (0,0) are already written by initializing with ret[subwindow]
                    if (f'{x_curr + slide_x}_{y_curr + slide_y}') in ret.keys():
                        # num_subw += 1
                        feature_matrix[subwindow] = {k: 
                                                    feature_matrix[subwindow].get(k, 0) + ret[f'{x_curr + slide_x}_{y_curr + slide_y}'].get(k, 0)
                                                    for k in set(feature_matrix[subwindow]).union(ret[f'{x_curr + slide_x}_{y_curr + slide_y}'])}
            # # scale the feature values by the number of summed subwindows that form it (it could be useful as feature vector normalization)
            # feature_matrix[subwindow] = {k:feature_matrix[subwindow][k]/num_subw for k in feature_matrix[subwindow].keys()}
            
            # scale the feature vector by the total numer of cells in it
            norm_factor = self.total_cell_norm/sum(feature_matrix[subwindow].values())
            for k in feature_matrix[subwindow]:
                feature_matrix[subwindow][k] = np.float32(feature_matrix[subwindow][k] * norm_factor)
                
        feature_matrix = pd.DataFrame(feature_matrix).T
        # feature_matrix is placd in AnnData object with specified spatial cooridnated of the sliding windows
        self.tissue = AnnData(feature_matrix.astype(np.float32), dtype=np.float32)
        self.tissue.obsm['spatial'] = np.array([[x.split('_')[0], x.split('_')[1]] for x in feature_matrix.index]).astype(int)

    def calculate_spatial_cell_type_metrics(self):
        # calculate cell type specific global metrics
        # [NOTE] this could go into a function of the class
        # entropy & scatteredness
        self.adata.obs['x_coor'] = (self.adata.obsm['spatial'][:,0])
        self.adata.obs['y_coor'] = (self.adata.obsm['spatial'][:,1])
        cx_min = np.min(self.adata.obs['x_coor'])
        cx_max = np.max(self.adata.obs['x_coor'])
        cy_min = np.min(self.adata.obs['y_coor'])
        cy_max = np.max(self.adata.obs['y_coor'])

        self.tissue.var['entropy'] = pd.Series(index=self.unique_cell_type, name='entropy', dtype=np.float64)
        self.tissue.var['scatteredness'] = pd.Series(index=self.unique_cell_type, name='scatteredness', dtype=np.float64)
        for cell_t in self.unique_cell_type:
            tissue_window = np.zeros(shape=(int(np.ceil((cx_max-cx_min+1)/self.downsample_rate)), int(np.ceil((cy_max-cy_min+1)/self.downsample_rate))), dtype=np.int8)
            tissue_window[((self.adata.obs['x_coor'][self.adata.obs[self.annotation] == cell_t] - cx_min)/self.downsample_rate).astype(int), ((self.adata.obs['y_coor'][self.adata.obs[self.annotation] == cell_t] - cy_min)/self.downsample_rate).astype(int)] = 1
            self.tissue.uns[f'tw_{cell_t}'] = tissue_window
            
            self.tissue.var['entropy'].loc[cell_t] = skimage.measure.shannon_entropy(tissue_window)
            _, num_objects = scipy.ndimage.measurements.label(tissue_window, structure=None, output=None) # this assumes 4 neighbors connectivity
            # # idea for scatteredness was to compute the number of connected components and divide it with number of existing non-zero elements
            # # but this measure does not contain the information on percentage of non-zero elements in the matrix.
            # # thus we multiply it with non-zero percentage (num non-zero / total num) creating just this formula
            # # num_object/image.size
            # # max value is based on neighbors size (if 4 then 1/4, if 8, 1/8), min value is 0 if there are no non-zero elements
            # [NOTE] add neighbourhood size for scatteredness calculation to params
            # [NOTE] try to find a heuristic to control the downsampling rate based on the proportion of cell number to area pixel number
            self.tissue.var['scatteredness'].loc[cell_t] = num_objects/tissue_window.size *4

    def cluster(self):

        sc.pp.neighbors(self.tissue, use_rep='X')
        sc.tl.leiden(self.tissue, resolution=self.resolution)

        bin_slide_ratio = int(self.win_size/self.sliding_step)
        x_min = self.adata.obs['Centroid_X'].min()
        y_min = self.adata.obs['Centroid_Y'].min()
        # max voting on cluster labels
        # init the new obs column
        self.tissue.obs['leiden_max_vote'] = list('x' for x in range(len(self.tissue.obs.index)))
        for x_curr, y_curr in self.tissue.obsm['spatial']:
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
            # other labels during refinment
            self.tissue.obs.loc[f'{x_curr}_{y_curr}', 'leiden_max_vote'] = max(subwindow_labels, key=subwindow_labels.get)

        self.adata.obs[f'tissue_{self.method_key}'] = list(self.tissue.obs.loc[self.adata.obs['x_y'], 'leiden_max_vote'])

        logging.info(r"Sliding window cell mixture calculation done. Added results to adata.obs['sliding_window']")

        
    def save_results(self):
        
        for cell_t in self.unique_cell_type:
            plt.imsave(fname=os.path.join(self.dir_path, f'tissue_window_{cell_t}_{self.params_suffix}.png'), arr=self.tissue.uns[f'tw_{cell_t}'], vmin=0, vmax=1, cmap='gray', dpi=250)
        
        # save metrics results in csv format
        # print(self.tissue.var[['entropy', 'scatteredness']])
        self.tissue.var[['entropy', 'scatteredness']].to_csv(os.path.join(self.dir_path, f'spatial_metrics_{self.params_suffix}.csv'))

        # save anndata file
        self.tissue.write_h5ad(os.path.join(self.dir_path, f'tissue_{self.filename}.h5ad'), compression="gzip")

        logging.info(f'Saved clustering result tissue_{self.filename}.h5ad.')

        if 'cell mixtures' in self.tissue.uns.keys():
            self.tissue.uns['cell mixtures'].to_csv(os.path.join(self.dir_path, f'cell_mixture_stats_{self.params_suffix}.csv'))
