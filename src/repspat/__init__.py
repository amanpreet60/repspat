__version__ = "0.1.0"

from .data import SampleData, to_binary
from .clustering import (
    cluster_feature_presence,
    custom_silhouette,
    create_blocks,
    spatial_silhouette_analysis,
    spatial_constrained_hac,
)
from .mmd import two_sample_mmd, multiple_comparison

__all__ = [
    "SampleData",
    "to_binary",
    "cluster_feature_presence",
    "custom_silhouette",
    "create_blocks",
    "spatial_silhouette_analysis",
    "spatial_constrained_hac",
    "two_sample_mmd",
    "multiple_comparison",
]
