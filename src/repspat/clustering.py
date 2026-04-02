import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans, AgglomerativeClustering
import squidpy as sq

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
        sq.gr.spatial_neighbors(
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
        if neighs:
            b = min(dist.loc[i, clusters == n].mean() for n in neighs)       # min mean distance to neighbors
        else:
            b = np.inf
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

def spatial_constrained_hac(adata, feature_df: pd.DataFrame, n_clusters: int = 7, 
                            n_neighs: int = 8, coord_type: str = "generic", delaunay: bool = False
):
    sq.gr.spatial_neighbors(
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
