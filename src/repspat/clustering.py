import numpy as np
import pandas as pd
import warnings
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans, AgglomerativeClustering


def _spatial_neighbors(*args, **kwargs):
    import squidpy as sq

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Calling `spatial_neighbors` is deprecated.*",
            category=FutureWarning,
        )
        return sq.gr.spatial_neighbors(*args, **kwargs)

def spatial_silhouette_analysis(sample_data, n_neighbors_list=[6,8], n_clusters_range=range(4,9)):
    results = []

    # Distance matrix of features
    dist_features_df = pd.DataFrame(
        sample_data.dist_matrix,
        index=sample_data.feature_mat.index,
        columns=sample_data.feature_mat.index
    )

    for knn in n_neighbors_list:
        # Construct spatial neighbors graph
        _spatial_neighbors(
            sample_data.sample_adata,
            n_neighs=knn,
            coord_type="generic",
            delaunay=False
        )

        adjacency = sample_data.sample_adata.obsp["spatial_connectivities"].toarray()
        adjacency_df = pd.DataFrame(
            adjacency,
            index=sample_data.sample_adata.obs_names,
            columns=sample_data.sample_adata.obs_names
        )
        connectivity_sparse = csr_matrix(adjacency)

        for n_clusters in n_clusters_range:
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric="euclidean",
                linkage="ward",
                connectivity=connectivity_sparse
            )

            cluster_labels = clustering.fit_predict(sample_data.sample_adata.X)
            clusters = pd.Series(cluster_labels, index=sample_data.feature_mat.index)

            sil_scores = custom_silhouette(clusters, dist_features_df, adjacency_df)
            avg_sil = np.mean(sil_scores)

            results.append({
                "n_neighbors": knn,
                "n_clusters": n_clusters,
                "avg_silhouette": avg_sil
            })

    return pd.DataFrame(results)


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

def create_blocks(feature_mat: pd.DataFrame, num_features: int, knn: int) -> pd.DataFrame:
    """Create KMeans blocks within each region."""
    feature_mat = feature_mat.copy()
    feature_mat['idx'] = np.arange(len(feature_mat))  # preserve original order
    blk_data = []

    for region in feature_mat['region'].unique():
        region_data = feature_mat[feature_mat['region'] == region].copy()
        df = region_data.iloc[:, :num_features]
        num_blks = len(df) // knn
        num_blks = 1 if num_blks == 0 or num_blks >= len(df.drop_duplicates()) else num_blks
        # Assign polygon IDs
        region_data['polygon_id'] = 1 if num_blks == 1 else KMeans(n_clusters=num_blks, n_init=10, random_state=0).fit(df).labels_ + 1
        blk_data.append(region_data)

    # Combine and restore original order
    return pd.concat(blk_data).sort_values('idx').drop(columns=['idx'])

def cluster_feature_presence(
    feature_mat: pd.DataFrame,
    clusters=None,
    cluster_column: str = "region",
    top_n: int = 5,
    min_presence: float = 0.0,
    feature_columns=None,
) -> pd.DataFrame:
    """Rank mostly present thresholded features within each cluster.

    Presence is calculated as the mean of each feature column inside a cluster.
    For binary thresholded columns, this is the fraction of cells where the
    marker is present.
    """
    data = feature_mat.copy()

    if clusters is None:
        if cluster_column not in data.columns:
            raise ValueError(
                f"cluster_column '{cluster_column}' not found. "
                "Pass clusters=... or provide a feature table with this column."
            )
        cluster_labels = data[cluster_column]
    else:
        cluster_labels = pd.Series(clusters, index=data.index, name=cluster_column)

    if feature_columns is None:
        excluded = {cluster_column, "polygon_id"}
        feature_columns = [
            col for col in data.select_dtypes(include=[np.number]).columns
            if col not in excluded
        ]

    if not feature_columns:
        raise ValueError("No numeric feature columns found to summarize.")

    rows = []
    for cluster_id in pd.Series(cluster_labels).dropna().unique():
        mask = cluster_labels == cluster_id
        cluster_features = data.loc[mask, feature_columns]
        cluster_size = int(mask.sum())
        presence_rates = cluster_features.mean(axis=0).sort_values(ascending=False)

        if min_presence > 0:
            presence_rates = presence_rates[presence_rates >= min_presence]

        if top_n is not None:
            presence_rates = presence_rates.head(top_n)

        for feature, presence_rate in presence_rates.items():
            rows.append({
                "cluster": cluster_id,
                "feature": feature,
                "presence_rate": float(presence_rate),
                "present_count": int(cluster_features[feature].sum()),
                "cluster_size": cluster_size,
            })

    return pd.DataFrame(rows)

def spatial_constrained_hac(adata, feature_df: pd.DataFrame, n_clusters: int = 7, 
                            n_neighs: int = 8, coord_type: str = "generic", delaunay: bool = False
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

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    labels = model.fit_predict(X) + 1
    feature_df['region'] = pd.Series(labels, index=feature_df.index).astype("category")

    return labels, feature_df, model
