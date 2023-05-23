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


@timeit
def plot_all_annotation(out_path, algo_list):
    number_of_samples = len(algo_list)
    if number_of_samples<=2:
        number_of_rows = 1
        number_of_columns = number_of_samples
    else:
        number_of_rows = 2 if number_of_samples % 2==0 else 1
        number_of_columns = number_of_samples // 2 if number_of_samples % 2==0 else number_of_samples
    figure, axes = plt.subplots(nrows=number_of_rows, ncols=number_of_columns, squeeze=False, layout='constrained')

    i, j = 0, 0
    for algo in algo_list:
        sc.pl.spatial(algo.adata, color=[algo.annotation], palette=algo.annotation_palette, spot_size=algo.spot_size, ax=axes[i, j], show=False, frameon=False, title="")
        axes[i, j].get_legend().remove()
        axes[i, j].set_title(f'{algo.filename}', fontsize=4, loc='center', wrap=True)

        i = i if j+1<number_of_columns else i+1
        j = j+1 if j+1<number_of_columns else 0

    figure.savefig(f'{out_path}/cell_type_per_slice.png', dpi=300, bbox_inches='tight')
    plt.close()


@timeit
def plot_all_cluster_label():
    pass


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