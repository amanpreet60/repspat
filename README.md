# repSpat

Spatial statistics and clustering for single-cell imaging data.

Built for analyzing how cell types are spatially distributed in tissue samples — things like whether two cell populations are more mixed or segregated than expected. Works with AnnData `.h5ad` files (the standard format for single-cell data).

## Install

```bash
pip install repspat
```

## Usage

```python
from repspat import SampleData, spatial_silhouette_analysis, spatial_constrained_hac
from repspat import plot_spatial_clusters, create_blocks, multiple_comparison, pairwise_results_to_matrix
```

### Load a sample

```python
data = SampleData(
    sample_column="sample_id",
    sample_name="Sample_04",
    adata_path="03_TNBC_2018_spe.h5ad",
    cell_type_column="mm"
)

data.summary()
# {
#   'feature_mat': (5381, 36),
#   'coords_mat':  (5381, 2),
#   'dist_matrix': (5381, 5381),
#   'cell_type':   (5381, 1),
#   'sample_adata':(5381, 36)
# }
```

### Find the right number of clusters

Runs spatially-aware silhouette analysis over a range of k and neighbor configs:

```python
results = spatial_silhouette_analysis(data, n_neighbors_list=[6, 8], n_clusters_range=range(4, 9))
#    n_neighbors  n_clusters  avg_silhouette
# 0            6           4        0.136721
# 5            8           4        0.161823
```

### Cluster cells

```python
labels, feature_df, model = spatial_constrained_hac(
    data.sample_adata,
    feature_df=data.feature_mat,
    n_clusters=7,
    n_neighs=8
)
```

### Plot

```python
plot_spatial_clusters(data.coords_mat.centroidX, data.coords_mat.centroidY, labels=labels)
```

### Run MMD tests between clusters

```python
blocked_data = create_blocks(feature_df, num_features=36, knn=8)
results_df = multiple_comparison(blocked_data, data.dist_matrix, kernel="IMQ")

# pairs that are not significantly different
print(results_df[results_df["adj_p"] >= 0.05])
```

### Similarity network

```python
matrix = pairwise_results_to_matrix(results_df)
```

See `example.ipynb` for a full walkthrough on a TNBC dataset.

## API

### `SampleData(sample_column, sample_name, adata_path=None, adata_obj=None, metric="euclidean", thresholds=None, cell_type_column="mm")`
Loads and subsets an AnnData object for one sample. Computes the feature matrix, spatial coordinates, and pairwise distance matrix. Pass `thresholds` as `{marker: cutoff}` to binarize continuous markers.

### `spatial_silhouette_analysis(sample_data, n_neighbors_list, n_clusters_range)`
Returns a DataFrame of silhouette scores across all `(n_neighbors, n_clusters)` combinations. Use this to pick clustering params.

### `spatial_constrained_hac(adata, feature_df, n_clusters, n_neighs, coord_type, delaunay)`
Ward HAC with a spatial connectivity constraint. Returns `(labels, feature_df, model)`.

### `create_blocks(feature_mat, num_features, knn)`
Splits regions into KMeans blocks for block permutation testing.

### `two_sample_mmd(sample1_idx, sample2_idx, dist_matrix, patient_data, kernel, kernel_param, nperm)`
MMD² between two groups with a permutation null. Returns observed statistic, null distribution, and p-value.

### `multiple_comparison(patient_data, dist_matrix, kernel, kernel_param, nperm, adj_p)`
Pairwise MMD across all cluster pairs. `adj_p` can be `"BH"`, `"bonferroni"`, or `"holm"`.

### `plot_spatial_clusters(x, y, labels)`
Scatter plot colored by cluster.

### `pairwise_results_to_matrix(df, plot)`
Builds an adjacency matrix and network graph from MMD results. Edges connect clusters that are not significantly different.

### `to_binary(column, marker_name, thresholds)`
Converts a continuous marker column to binary using a threshold dict.

## Requirements

- Python >= 3.9
- numpy, pandas, scipy, scikit-learn, statsmodels, scanpy, squidpy >= 1.2, networkx, matplotlib

## License

MIT
