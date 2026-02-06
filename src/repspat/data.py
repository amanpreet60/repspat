import pandas as pd
import scanpy as sc
import numpy as np
from scipy.spatial.distance import pdist, squareform
import warnings

def to_binary(column: pd.Series, marker_name: str, thresholds: dict) -> pd.Series:
    """Convert continuous marker values to binary using provided threshold dict."""
    if marker_name not in thresholds:
        warnings.warn(f"No threshold defined for marker '{marker_name}'. Column left unchanged.")
        return column
    
    threshold = thresholds[marker_name]
    return (column >= threshold).astype(int)


class SampleData:
    def __init__(
        self,
        sample_column,
        sample_name,
        adata_path=None,
        adata_obj=None,
        metric="euclidean",
        thresholds=None
    ):
        """
        Parameters
        ----------
        thresholds : dict or None
            User-provided threshold map {marker_name: cutoff}.
            If None, no binarization is applied.
        """
        # Load AnnData
        if adata_obj is not None:
            adata = adata_obj
        elif adata_path is not None:
            adata = sc.read_h5ad(adata_path)
        else:
            raise ValueError("Provide either adata_path or adata_obj")
        
        # Subset sample
        sample = adata[adata.obs[sample_column] == sample_name].copy()
        
        # Feature matrix
        X = sample.X.A if hasattr(sample.X, "A") else sample.X
        self.feature_mat = pd.DataFrame(X, index=sample.obs_names, columns=sample.var_names)

        # ---- Apply thresholds ONLY if user provided them ----
        if thresholds is not None:
            for col in self.feature_mat.columns:
                if col in thresholds:
                    self.feature_mat[col] = to_binary(self.feature_mat[col], col, thresholds)

        # Spatial coordinates
        self.coords_mat = sample.obsm["spatial"]

        # Distance matrix
        self.dist_matrix = pd.DataFrame(
            squareform(pdist(self.feature_mat, metric=metric)),
            index=self.feature_mat.index,
            columns=self.feature_mat.index
        )

        # Cell type
        self.cell_type = sample.obs[["mm"]].copy()
        self.cell_type.rename(columns={"mm": "cell_type"}, inplace=True)

        # Store AnnData subset
        self.sample_adata = sample

    def summary(self):
        return {
            "feature_mat": self.feature_mat.shape,
            "coords_mat": self.coords_mat.shape,
            "dist_matrix": self.dist_matrix.shape,
            "cell_type": self.cell_type.shape,
            "sample_adata": self.sample_adata.shape,
        }
