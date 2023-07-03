import time
import seaborn as sns

from matplotlib import rcParams
from matplotlib.axes import Axes
from functools import wraps


import scanpy as sc

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
    annotation,
    ax: Axes,
    spot_size: float,
    palette = None,
    title: str = "",
    groups=None
):
    """
    Scatter plot in spatial coordinates. A wrapper around scanpy's pl.spatial to correct inversion of y-axis. 
    Standard plots have coordinates in 'lower left' while images are considered
    to start from 'upper left'. Scanpy assumes that we will always
    plot some sort of staining image. If we do not provide image, scanpy
    will flip yaxis for us to get back to the standard plot coordinates.
    That causes inversion of our plotting figures so we wrapped scanpy's
    pl.spatial. utils.py is created for such wrappers, timeit decorator and
    everything else that should be globally accessible and does not belong
    to any specific class.

    """
    spot_size = spot_size * 0.5
    data = adata[adata.obs[annotation].isin(groups)] if groups else adata
    ax = sns.scatterplot(data=data.obs, hue=annotation, x=data.obsm['spatial'][:, 0], y=data.obsm['spatial'][:, 1], ax=ax, s=spot_size, linewidth=0, palette=palette, marker='.')
    ax.set(yticklabels=[], xticklabels=[], title=title)
    ax.tick_params(bottom=False, left=False)
    ax.set_aspect("equal")
    sns.despine(bottom=True, left=True, ax=ax)

def set_figure_params(
        dpi: int,
        facecolor: str,
):
    rcParams['figure.facecolor'] = facecolor
    rcParams['axes.facecolor'] = facecolor
    rcParams["figure.dpi"] = dpi
