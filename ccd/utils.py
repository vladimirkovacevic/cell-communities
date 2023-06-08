import time
import scanpy as sc

from matplotlib.axes import Axes
from typing import List
from functools import wraps

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