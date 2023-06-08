import time
import scanpy as sc
import numpy as np

from matplotlib.axes import Axes
from typing import List
from functools import wraps
from collections import defaultdict

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
def plot_spatial(
    adata,
    color: List[str],
    ax: Axes,
    spot_size: float,
    title: str = None,
    groups = None,
    palette: List[str] = None,
    show: bool = False,
    frameon: bool = False
):
    """
    Scatter plot in spatial coordinates. A wrapper around scanpy's pl.spatial to correct inversion of y-axis.

    """
    ax.invert_yaxis()
    sc.pl.spatial(adata, color=color, spot_size=spot_size, ax=ax, show=show, frameon=frameon, title=title, palette=palette, groups=groups)

@timeit
def calc_optimal_win_size_and_slide_step(slice):
    x_min, x_max = np.min(slice.obsm['spatial'][0]), np.max(slice.obsm['spatial'][0])
    y_min, y_max = np.min(slice.obsm['spatial'][1]), np.max(slice.obsm['spatial'][1])
    x_range, y_range = abs(abs(x_max) - abs(x_min)), abs(abs(y_max) - abs(y_min))

    win_size = int(x_range // 100 if x_range < y_range else y_range // 100)
    win_size = win_size + 1 if win_size & 1 else win_size

    avg_covered = -1
    iters = 0
    MAX_ITERS = 10 
    MIN_COVERED = 30
    MAX_COVERED = 50
    while (avg_covered < MIN_COVERED or avg_covered > MAX_COVERED) and iters < MAX_ITERS :
        cell_to_loc = defaultdict(int)
        for x, y in slice.obsm['spatial']:
            cell_to_loc[(x // win_size, y // win_size)] += 1
        
        cell_to_loc.items()

        iters += 1



