__version__ = "0.1.0"

from .data import SampleData, to_binary
from .clustering import (
    custom_silhouette,
    create_blocks,
    spatial_silhouette_analysis,
    spatial_constrained_hac,
)
from .mmd import two_sample_mmd, multiple_comparison
from . import visualization
from .visualization import pairwise_results_to_matrix, plot_spatial_clusters

__all__ = [
    "SampleData",
    "to_binary",
    "custom_silhouette",
    "create_blocks",
    "spatial_silhouette_analysis",
    "spatial_constrained_hac",
    "two_sample_mmd",
    "multiple_comparison",
    "visualization",
    "plot_spatial_clusters",
    "pairwise_results_to_matrix",
]
