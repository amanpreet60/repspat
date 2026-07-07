import pandas as pd
from scipy.spatial.distance import pdist, squareform

def _as_array(matrix):
    """Return a dense array from dense, sparse, or matrix-like AnnData storage."""
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    if hasattr(matrix, "A"):
        return matrix.A
    return matrix

def _feature_array(adata):
    layer = adata.uns.get("repspat", {}).get("layer")
    if layer is not None:
        return _as_array(adata.layers[layer])
    return _as_array(adata.X)

class SampleData:
    def __init__(
        self,
        adata,
        sample_column=None,
        sample_name=None,
        *,
        metric,
        layer=None,
        cell_type_column="mm"
    ):
        """
        Parameters
        ----------
        adata
            AnnData object, or path to an .h5ad file.
        sample_column : str or None
            Column in adata.obs that identifies samples. If None, use all rows.
        sample_name : str or None
            Sample value to select from sample_column. If None, use all rows.
        metric : str
            Distance metric to use for the pairwise feature distance matrix.
        layer : str or None
            Name of an AnnData layer to use as the feature matrix instead of
            adata.X. Useful when thresholded/binary data is already stored in
            the input file.
        """
        if isinstance(adata, str):
            import scanpy as sc

            adata = sc.read_h5ad(adata)

        if (sample_column is None) != (sample_name is None):
            raise ValueError("sample_column and sample_name must be provided together.")

        if sample_column is not None:
            sample = adata[adata.obs[sample_column] == sample_name].copy()
        else:
            sample = adata

        # Feature matrix
        if layer is not None:
            if layer not in sample.layers:
                raise KeyError(f"Layer '{layer}' not found in AnnData layers.")
            X = _as_array(sample.layers[layer])
        else:
            X = _as_array(sample.X)

        self.feature_mat = pd.DataFrame(X, index=sample.obs_names, columns=sample.var_names)

        # Spatial coordinates
        coords = sample.obsm["spatial"]
        if hasattr(coords, "to_numpy"):
            coords = coords.to_numpy()
            sample.obsm["spatial"] = coords
        self.coords_mat = coords

        # Distance matrix
        dist_values = squareform(pdist(self.feature_mat, metric=metric))
        self.dist_matrix = pd.DataFrame(
            dist_values,
            index=self.feature_mat.index,
            columns=self.feature_mat.index
        )
        sample.obsp["repspat_distances"] = dist_values
        sample.uns["repspat"] = {
            "layer": layer,
            "metric": metric,
            "sample_column": sample_column,
            "sample_name": sample_name,
            "distance_key": "repspat_distances",
            "region_key": "repspat_region",
            "polygon_key": "repspat_polygon_id",
        }

        # Cell type
        self.cell_type = sample.obs[[cell_type_column]].copy()
        self.cell_type.rename(columns={cell_type_column: "cell_type"}, inplace=True)

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
