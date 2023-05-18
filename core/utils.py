from functools import wraps
import time
import seaborn as sns
import scanpy as sc
from matplotlib import pyplot as plt
from itertools import cycle
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
def plot_cell_perc_in_community_per_slice(df, path):
    sc.settings.set_figure_params(dpi=300, facecolor='white')
    sns.set(font_scale=1)

    ncols = len(df.columns) # we want to separately print the total_counts column
    fig, axes = plt.subplots(ncols=ncols, figsize=(15,15))

    fig.subplots_adjust(wspace=0)

    # cmap_cycle = cycle(['Greens', 'Reds', 'Blues', 'Oranges', 'Purples'])

    for i, ax in enumerate(axes):
        g = sns.heatmap(pd.DataFrame(df.iloc[:, i]), vmin=df.to_numpy().min(), vmax=df.to_numpy().max(), linewidths=0, linecolor=None, annot=True, cbar=False, ax=ax, \
                        cmap="Greys", fmt='.2%', xticklabels=True, yticklabels=True if i==0 else False, square=True)
        g.set_xticklabels(g.get_xticklabels(), rotation=70)
        g.xaxis.tick_top()
    plt.savefig(os.path.join(path, 'cell_perc_in_community_per_slice.png'), dpi=400)
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