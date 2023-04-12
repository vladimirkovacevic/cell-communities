import sys
import logging

from .community_clustering_algorithm import CommunityClusteringAlgo
from .utils import timeit

try:
    from .sliding_window import SlidingWindow
except ImportError:
    logging.warn("Module SlidingWindow is not present.")
