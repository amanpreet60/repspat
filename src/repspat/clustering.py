import numpy as np
import pandas as pd
import warnings
from sklearn.cluster import KMeans, AgglomerativeClustering
from .data import _as_array, _feature_array


def _spatial_neighbors(*args, **kwargs):
    import squidpy as sq

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Calling `spatial_neighbors` is deprecated.*",
            category=FutureWarning,
        )
        return sq.gr.spatial_neighbors(*args, **kwargs)

def spatial_silhouette_analysis(adata, n_neighbors_list=[6,8], n_clusters_range=range(4,9)):
    results = []

    # Distance matrix of features
    dist_matrix = _as_array(adata.obsp["repspat_distances"])
    feature_values = _feature_array(adata)

    for knn in n_neighbors_list:
        # Construct spatial neighbors graph
        _spatial_neighbors(
            adata,
            n_neighs=knn,
            coord_type="generic",
            delaunay=False
        )

        connectivity_sparse = adata.obsp["spatial_connectivities"].tocsr()
        adjacency = connectivity_sparse.toarray().astype(bool)

        for n_clusters in n_clusters_range:
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric="euclidean",
                linkage="ward",
                connectivity=connectivity_sparse
            )

            cluster_labels = clustering.fit_predict(feature_values)
            labels = np.asarray(cluster_labels)
            silhouettes = np.zeros(labels.shape[0], dtype=float)
            label_indices = {
                label: np.flatnonzero(labels == label)
                for label in np.unique(labels)
            }
            neighbor_labels = {}
            for label, indices in label_indices.items():
                connected = adjacency[indices].any(axis=0)
                neighbor_labels[label] = [
                    other
                    for other in np.unique(labels[connected])
                    if other != label
                ]

            for label, members in label_indices.items():
                if len(members) > 1:
                    intra = dist_matrix[np.ix_(members, members)]
                    a = (intra.sum(axis=1) - np.diag(intra)) / (len(members) - 1)
                else:
                    a = np.zeros(len(members), dtype=float)

                neighbors = neighbor_labels.get(label, [])
                if not neighbors:
                    continue

                b_candidates = [
                    dist_matrix[np.ix_(members, label_indices[neighbor])].mean(axis=1)
                    for neighbor in neighbors
                ]
                b = np.minimum.reduce(b_candidates)
                denom = np.maximum(a, b)
                silhouettes[members] = np.divide(
                    b - a,
                    denom,
                    out=np.zeros_like(denom, dtype=float),
                    where=denom != 0,
                )

            sil_scores = silhouettes
            avg_sil = np.mean(sil_scores)

            results.append({
                "n_neighbors": knn,
                "n_clusters": n_clusters,
                "avg_silhouette": avg_sil
            })

    results_df = pd.DataFrame(results)
    adata.uns["repspat_silhouette"] = results_df
    return results_df


def create_blocks(adata, knn: int, region_key: str = "repspat_region",
                  polygon_key: str = "repspat_polygon_id") -> pd.DataFrame:
    """Create KMeans blocks within each region."""
    if region_key not in adata.obs:
        raise KeyError(f"'{region_key}' not found in adata.obs.")

    feature_df = pd.DataFrame(
        _feature_array(adata),
        index=adata.obs_names,
        columns=adata.var_names,
    )
    polygon_ids = pd.Series(index=adata.obs_names, dtype="Int64")

    for region in adata.obs[region_key].unique():
        obs_names = adata.obs_names[adata.obs[region_key] == region]
        df = feature_df.loc[obs_names]
        num_blks = len(df) // knn
        num_blks = 1 if num_blks == 0 or num_blks >= len(df.drop_duplicates()) else num_blks
        if num_blks == 1:
            polygon_ids.loc[obs_names] = 1
        else:
            polygon_ids.loc[obs_names] = (
                KMeans(n_clusters=num_blks, n_init=10, random_state=0).fit(df).labels_ + 1
            )

    adata.obs[polygon_key] = polygon_ids
    adata.uns.setdefault("repspat", {})["polygon_key"] = polygon_key
    return adata.obs[[region_key, polygon_key]].copy()


def spatial_constrained_hac(adata, n_clusters: int = 7, n_neighs: int = 8,
                            coord_type: str = "generic", delaunay: bool = False,
                            region_key: str = "repspat_region",
):
    _spatial_neighbors(
        adata,
        n_neighs=n_neighs,
        coord_type=coord_type,
        delaunay=delaunay,
    )

    connectivity = adata.obsp["spatial_connectivities"].tocsr()

    model = AgglomerativeClustering(
        n_clusters=n_clusters,
        linkage="ward",
        connectivity=connectivity,
        compute_distances=True,
    )

    X = _as_array(adata.X)
    labels = model.fit_predict(X) + 1
    adata.obs[region_key] = pd.Series(labels, index=adata.obs_names).astype("category")
    adata.uns.setdefault("repspat", {})["region_key"] = region_key
    adata.uns["repspat_hac"] = {
        "n_clusters": n_clusters,
        "n_neighs": n_neighs,
        "coord_type": coord_type,
        "delaunay": delaunay,
    }

    return labels, adata, model
