# src/spatialstats/__init__.py

# Data-related
from .data import SampleData

# Clustering functions
from .clustering import custom_silhouette, create_blocks, spatial_silhouette_analysis

# MMD functions
from .mmd import two_sample_mmd, multiple_comparison, pairwise_results_to_matrix

# Optional: define __all__ for clarity
__all__ = [
    "SampleData", "to_binary",
    "custom_silhouette", "create_blocks",
    "two_sample_mmd", "multiple_comparison", "pairwise_results_to_matrix"
]
