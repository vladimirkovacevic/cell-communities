import logging

from .community_clustering_algorithm import CommunityClusteringAlgo
from .metrics import *
from .html_report import generate_report
from .metrics import calculate_spatial_metrics
from .community_detection import CommunityDetection

try:
    from .sliding_window import SlidingWindow
except ImportError:
    logging.warn("Module SlidingWindow is not present.")


try:
    from .sliding_window import SlidingWindowMultipleSizes
except ImportError:
    logging.warn("Module SlidingWindowMultipleSizes is not present.")
