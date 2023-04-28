import logging

from .community_clustering_algorithm import CommunityClusteringAlgo
from .utils import timeit
from .spatial_metrics import calculate_spatial_metrics

try:
    from .sliding_window import SlidingWindow
except ImportError:
    logging.warn("Module SlidingWindow is not present.")

try:
    from .sliding_window_multiple_sizes import SlidingWindowMultipleSizes
except ImportError:
    logging.warn("Module SlidingWindowMultipleSizes is not present.")
