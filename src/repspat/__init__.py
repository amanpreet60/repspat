# src/spatialstats/__init__.py

# Data-related
from .data import SampleData, to_binary

# Clustering functions
from .clustering import custom_silhouette, create_blocks, spatial_silhouette_analysis, spatial_constrained_hac

# MMD functions
from .mmd import two_sample_mmd, multiple_comparison

from .visualization import plot_spatial_clusters, pairwise_results_to_matrix

# Optional: define __all__ for clarity
__all__ = [
    "SampleData", "to_binary",
    "custom_silhouette", "create_blocks", "spatial_silhouette_analysis", "spatial_constrained_hac",
    "two_sample_mmd", "multiple_comparison",
    "plot_spatial_clusters", "pairwise_results_to_matrix"
]
