import time
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
import os
import math
import matplotlib.ticker as mtick

from functools import wraps
from functools import reduce
from matplotlib import pyplot as plt



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
def celltype_mixtures_total_plot(cell_mixtures, path):
    def merge_dicts(dict1, dict2):
        return { cluster: dict1.get(cluster, 0) + dict2.get(cluster, 0) for cluster in set(dict1) | set(dict2) }
    def merge_dicts_of_dicts(dict1, dict2):
        return { celltype: merge_dicts(dict1.get(celltype, {}), dict2.get(celltype, {})) for celltype in set(dict1) | set(dict2) }

    total_dict = reduce(merge_dicts_of_dicts, cell_mixtures)
    total = pd.DataFrame(total_dict)

    total['total_counts'] = np.array([sum(total.loc[row, :]) for row in total.index]).astype(int)

    cell_type_counts = {ct:[int(sum(total[ct]))] for ct in total.columns}
    total = pd.concat([total, pd.DataFrame(cell_type_counts, index=['total_cells'])])

    total.iloc[:-1, :-1] = total.iloc[:-1, :-1].div(total['total_counts'][:-1], axis=0).mul(100)
    total['perc_of_all_cells'] = np.around(total['total_counts'] / total['total_counts'][-1] * 100, decimals=1)
    total = total.loc[sorted(total.index.values, key=lambda x: float(x) if x != "total_cells" else float('inf'))]

    sc.settings.set_figure_params(dpi=300, facecolor='white')
    sns.set(font_scale=1.5)

    ncols = len(total.columns)
    fig, axes = plt.subplots(ncols=ncols, figsize=(15,15))
    fig.subplots_adjust(wspace=0)

    vmax_perc = np.max(np.ravel(total.iloc[:-1,:-2]))
    for i, ax in enumerate(axes[:-2]):
        sns.heatmap(pd.DataFrame(total.iloc[:, i]), vmin=0.0, vmax=vmax_perc, linewidths=0, linecolor=None, annot=True, cbar=False, ax=ax, \
                        cmap="Greys", fmt='4.0f', xticklabels=True, yticklabels=True if i==0 else False, square=True)
    sns.heatmap(pd.DataFrame(total.iloc[:, -2]), annot=True, vmin=0, vmax=np.max(total.iloc[:-1, -2]), linewidths=0, linecolor=None, \
        cbar=False, cmap='Greens', ax=axes[-2], fmt='4.0f', xticklabels=True, yticklabels=False, square=True)
    sns.heatmap(pd.DataFrame(total.iloc[:, -1]), annot=True, vmin=0, vmax=np.max(total.iloc[:-1, -1]), linewidths=0, linecolor=None, cbar=False, \
        cmap='Greens', ax=axes[-1], fmt='4.0f', xticklabels=True, yticklabels=False, square=True)
    
    for ax in axes:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=70)
        ax.xaxis.tick_top() 

    plt.savefig(os.path.join(path, f'total_cell_mixtures_table.png'), dpi=400)
    plt.close()


@timeit
def plot_cell_perc_in_community_per_slice(df, path):
    sc.settings.set_figure_params(dpi=200, facecolor='white')
    sns.set(font_scale=1.5)
    plt.figure(figsize=(15, 15))

    ax = sns.heatmap(df, annot=True, fmt="4.0f", cmap="Greys", xticklabels=True, yticklabels=True, square=True, cbar=False)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.savefig(os.path.join(path, 'cell_perc_in_community_per_slice.png'))
    plt.close()

@timeit
def plot_cell_abundance_total(algos, path):
    fig, ax = plt.subplots(figsize=(20,10))
    fig.subplots_adjust(wspace=0)

    cell_percentage_dfs = []
    plot_columns = []
    for algo in algos:
        cell_percentage_dfs.append(pd.DataFrame(algo.get_adata().obs['author_cell_type'].value_counts(normalize=True).mul(100).rename(algo.filename)))
        plot_columns.append(algo.filename)

    cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0).reset_index()
    cummulative_df.plot(x="index", y=plot_columns, kind="bar", rot=70, ax=ax, xlabel="")

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax.grid(False)
    ax.set_facecolor('white')
    plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig(os.path.join(path, f'cell_abundance_all_slices.png'), dpi=400)
    plt.close()

@timeit
def cell_abundance_per_slice(algos, path):
    # columns = rows = math.ceil(math.sqrt(len(samples)))
    # plt.figure(figsize=(10, 10))
    # for i in range(1, columns*rows +1):
    #     ax = plt.subplot(2,2,i)
    #     sns.countplot(x=samples[i], dodge=False, ax=ax, hue=samples[i])
    #     ax.set(xticklabels=[])
    # plt.show()
    pass

@timeit 
def cluster_abundance_total(algos, path):
    fig, ax = plt.subplots(figsize=(20,10))
    fig.subplots_adjust(wspace=0)

    cell_percentage_dfs = []
    plot_columns = []
    for algo in algos:
        cell_percentage_dfs.append(pd.DataFrame(algo.get_tissue().obs['leiden'].value_counts(normalize=True).mul(100).rename(algo.filename)))
        plot_columns.append(algo.filename)

    cummulative_df = pd.concat(cell_percentage_dfs, axis=1).fillna(0).reset_index()
    cummulative_df.plot(x="index", y=plot_columns, kind="bar", rot=70, ax=ax, xlabel="")

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax.grid(False)
    ax.set_facecolor('white')
    plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig(os.path.join(path, f'cluster_abundance_all_slices.png'), dpi=400)
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