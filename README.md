# Cell Communities Clusterization on ST Data

Sequencing technology advances applied in spatial transcriptomics opened a new window into the tissue microenvironment. By capturing gene expression profiles of individual cells along with their spatial coordinates, researchers have gained the ability to accurately ascertain cell types and their respective functions. While some cell types are known to coexist within various regions of the tissue, the co-occurrence of many others remains elusive. The proposed Cell Community Detection (CCD) algorithm addresses this challenge by offering a novel computational approach for identifying tissue domains with a significant mixture of particular cell types. The CCD algorithm divides the tissue using sliding windows, quantifies the proportion of each cell type within each window, and groups together windows with similar cell type mixtures. By employing majority voting, the algorithm assigns a community label to each cell based on the covered windows. Notably, CCD accommodates multiple window sizes and enables the simultaneous analysis of multiple slides from the same tissue. Its Python implementation with a flexible user interface enables the processing of datasets with tens of thousands of cells in sub-minute execution time.

**Keywords**: spatial transcriptomics, cell co-occurance, cell community, bioinformatics,

## Prerequisites
- Python 3.7 or higher
- leidenalg
- scikit-image
- matplotlib
- numpy
- pandas
- scanpy
- scikit-learn
- scipy
- seaborn
## Installation

You can install the required Python packages by running:

```

pip install -r requirements.txt
```


## Usage

```css

python main.py [-h] -f FILE [-t TFILE] -a ANNOTATION -m METHODS [-o OUT_PATH] [-r RESOLUTION] [-s SPOT_SIZE] [-v VERBOSE] [--total_cell_norm TOTAL_CELL_NORM] [--downsample_rate DOWNSAMPLE_RATE] [--entropy_thres ENTROPY_THRES] [--scatter_thres SCATTER_THRES] [-w WIN_SIZE] [--sliding_step SLIDING_STEP]
```


### Required arguments: 
- `-f, --files`: csv list of file paths to file that contain data to be analyzed/clustered. 
- `-a, --annotation`: Annotation label for cell types. 
### Optional arguments: 
- `-t, --tfile`: File path to Anndata object with calculated cell mixtures for data windows, output of `calc_feature_matrix`. 
- `-o, --out_path`: Absolute path to store outputs.
- `-c, --cluster_algo`: Clustering algorithm for cell communities. Available: leiden, spectral, agglomerative. Default is `leiden`.
- `-r, --resolution`: Resolution of leiden clustering algorithm. Ignored for spectral and agglomerative. Default is `0.2`. 
- `-s, --spot_size`: Size of the spot on plot. Default is `30`. 
- `-v, --verbose`: Show logging messages. `0` shows warnings, `>0` shows info, `<0` shows no output generated. Default is `0`. 
- `p, --plotting`: Save plots flag.`0`- No plotting/saving, `1` - save clustering plot, `2` - save all plots (cell type images statistics, and cell mixture plots). Default is `2`.
- `--skip_stats`: Skip statistics calculation on cell community clustering result. A table of cell mixtures and comparative spatial plots of cell types and mixtures will not be created. Default is False.
- `--total_cell_norm`: Total number of cells per window mixture after normalization. Default is `10000`.
- `--downsample_rate`: Rate by which the binary image of cells is downsampled before calculating the entropy and scatteredness metrics. Default is `50`.
- `--num_threads`: Number of threads that will be used to speed up community calling. Default is `5`.
- `--entropy_thres`: Threshold value for spatial cell type entropy for filtering out overdispersed cell types. Default is `1.0`.
- `--scatter_thres`: Threshold value for spatial cell type scatteredness for filtering out overdispersed cell types. Default is `1.0`.
- `-w, --win_sizes`: Comma-separated list of window sizes for analyzing the cell community. Default is `NA`. If window sizes are not provided, the algorithm uses iterative approach to obtain an optimal size of window base on average cell number.
- `--sliding_steps`: Comma-separated list of sliding steps for sliding window. Default is `NA`. If sliding step is not provided it will be defined as half of window size.
- `--min_cluster_size`: Minimum number of cells for a cluster to be plotted in the plot_stats() function. Default value is `200`.
- `--min_perc_to_show`: Minimum percentage of a cell type within a cluster for the cell type to be plotted in the plot_stats() function. Default value is `4`.
- `--min_cells_coeff`: Multiple of standard deviations from mean values that determines the cutoff for a certain value. Default value is `1.5`.
- `--save_adata`: Save Anndata file with resulting .obs column of cell community labels. Default value is False.
- `--min_count_per_type`: Minimum number of cells per cell type needed to use the cell type for cell communities extraction (in percentages). Default is `0.1`.
- `--project_name`: Project name that is used to name a directory containing all the slices used. Default is `community`.
- `--min_num_celltype`: Minimum number of cell types that have more than --min_perc_celltype in a cluster, for a cluster to be shown in plot_celltype_table(). Default is `1`.
- `--min_perc_celltype`: Minimum percentage of cells of a cell type which at least min_num_celltype cell types need to have to show a cluster in plot_celltype_table(). Default is `10`.
- `--color_plot_system`: Color system for display of cluster specific windows. Available: rgb, hsv. Default is `rgb`.
- `--n_clusters`: Number of clusters for spectral and agglomerative clustering. Ignored for leiden. Default is `10`.
- `--hide_plots`, Stop plots from displaying in notebooks or standard ouput. Used for batch processing. Default is `False`.

## Example

#### Download data: 
```css
wget https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E1S3_cell_bin_whole_brain.h5ad
```

#### Run the algorithm:
```css
python main.py -f E16.5_E1S3_cell_bin_whole_brain.h5ad -o results/whole_brain -a "sim anno" --scatter_thres 0.12 --resolution 0.25 --downsample_rate 80 --plotting 3
```


