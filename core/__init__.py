import logging

from .community_clustering_algorithm import CommunityClusteringAlgo
from .utils import timeit, plot_cell_perc_in_community_per_slice, celltype_mixtures_total_plot
from .spatial_metrics import calculate_spatial_metrics

try:
    from .sliding_window import SlidingWindow
except ImportError:
    logging.warn("Module SlidingWindow is not present.")


try:
    from .sliding_window import SlidingWindowMultipleSizes
except ImportError:
    logging.warn("Module SlidingWindowMultipleSizes is not present.")
