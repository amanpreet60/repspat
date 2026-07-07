import numpy as np
import pandas as pd
import warnings
from scipy.sparse import csr_matrix
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
    dist_features_df = pd.DataFrame(
        dist_matrix,
        index=adata.obs_names,
        columns=adata.obs_names
    )

    for knn in n_neighbors_list:
        # Construct spatial neighbors graph
        _spatial_neighbors(
            adata,
            n_neighs=knn,
            coord_type="generic",
            delaunay=False
        )

        adjacency = adata.obsp["spatial_connectivities"].toarray()
        adjacency_df = pd.DataFrame(
            adjacency,
            index=adata.obs_names,
            columns=adata.obs_names
        )
        connectivity_sparse = csr_matrix(adjacency)

        for n_clusters in n_clusters_range:
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric="euclidean",
                linkage="ward",
                connectivity=connectivity_sparse
            )

            cluster_labels = clustering.fit_predict(_as_array(adata.X))
            clusters = pd.Series(cluster_labels, index=adata.obs_names)

            sil_scores = custom_silhouette(clusters, dist_features_df, adjacency_df)
            avg_sil = np.mean(sil_scores)

            results.append({
                "n_neighbors": knn,
                "n_clusters": n_clusters,
                "avg_silhouette": avg_sil
            })

    results_df = pd.DataFrame(results)
    adata.uns["repspat_silhouette"] = results_df
    return results_df


def custom_silhouette(clusters, dist_matrix, adjacency):
    # Convert distance and adjacency matrices to DataFrames if they are NumPy arrays
    dist = (pd.DataFrame(dist_matrix, index=clusters.index, columns=clusters.index)
            if isinstance(dist_matrix, np.ndarray) else dist_matrix)
    adj = (pd.DataFrame(adjacency, index=clusters.index, columns=clusters.index)
           if isinstance(adjacency, np.ndarray) else adjacency)

    # Precompute neighboring clusters for each cluster based on adjacency
    neighbors = {}
    for cl in clusters.unique():
        obs = clusters[clusters == cl].index  # indices of current cluster
        connected = adj.loc[obs].any()        # boolean series of connected cells
        neighbors[cl] = [c for c in clusters[connected].unique() if c != cl]  # neighbor clusters

    silhouettes = []
    for i in clusters.index:
        cl = clusters[i]                       # cluster of current observation
        members = clusters[clusters == cl].index.drop(i, errors="ignore")  # other members in same cluster
        a = dist.loc[i, members].mean() if len(members) else 0.0            # mean intra-cluster distance
        neighs = neighbors.get(cl, [])                                        # neighboring clusters
        if not neighs:
            silhouettes.append(0.0)
            continue

        b = min(dist.loc[i, clusters == n].mean() for n in neighs)           # min mean distance to neighbors
        silhouettes.append((b - a) / max(a, b) if max(a, b) else 0.0)        # silhouette score for i

    return np.array(silhouettes)  # return all silhouette scores as NumPy array

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
