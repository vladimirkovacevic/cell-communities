from functools import wraps
import time
import seaborn as sns
import scanpy as sc
import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
import os


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__} took {total_time:.4f}s')
        return result
    return timeit_wrapper

def plot_all_slices(out_path, algo_list, annotation, img_name, clustering=False):
    number_of_samples = len(algo_list)
    number_of_rows = 2 if number_of_samples%2==0 and number_of_samples>2 else 1
    number_of_columns = (number_of_samples // 2) if number_of_samples % 2==0 and number_of_samples>2 else number_of_samples

    figure, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, squeeze=False, layout='constrained')
    axes_list = axes.flatten()
    h_d = {}
    handles, labels = [], []
    for (algo, ax) in zip(algo_list, axes_list):
        palette = algo.cluster_palette if clustering else algo.annotation_palette
        sc.pl.spatial(algo.adata, color=[annotation], palette=palette, spot_size=algo.spot_size, ax=ax, show=False, frameon=False)
        ax.get_legend().remove()
        ax.set_title(f'{algo.filename}', fontsize=6, loc='center', wrap=True)
        hands, labs = ax.get_legend_handles_labels()
        for h, l in zip(hands, labs):
            h._sizes = [11]
            if l not in h_d.values():
                h_d[h] = l
    
    for h, l in h_d.items():
        handles.append(h) 
        labels.append(l)
    legend_ncols = 1 if len(labels)<=12 else 2
    figure.legend(handles, labels, bbox_to_anchor=(1.15, 0.5), loc='center', fontsize=4, frameon=False, borderaxespad=0., ncols=legend_ncols, labelspacing=1, scatterpoints=10)
    figure.savefig(f'{out_path}/{img_name}', dpi=150, bbox_inches='tight')
    plt.close()


@timeit
def plot_all_annotation(out_path, algo_list):
    plot_all_slices(out_path, algo_list, algo_list[0].annotation, 'cell_type_per_slice.png')
    

@timeit
def plot_all_clustering(out_path, algo_list):
    plot_all_slices(out_path, algo_list, f'tissue_{algo_list[0].method_key}', 'clustering_per_slice.png', True)


@timeit
def plot_cell_perc_in_community_per_slice(df, path):
    sc.settings.set_figure_params(dpi=200, facecolor='white')
    sns.set(font_scale=2)
    plt.figure(figsize=(15, 15))

    ax = sns.heatmap(df, annot=True, fmt=".2%", cmap="Greys", xticklabels=True, yticklabels=True, square=True, cbar=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.savefig(os.path.join(path, 'cell_perc_in_community_per_slice.png'))
    plt.close()

@timeit
def calculate_something(num):
    """
    Simple function that returns sum of all numbers up to the square of num. To demonstrate the use of timeit decorator
    """
    total = sum([x for x in range(0, num**2)])
    return total

if __name__ == '__main__':
    calculate_something(10)
    calculate_something(100)
    calculate_something(10000)