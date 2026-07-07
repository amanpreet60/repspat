# repSpat

A Python package for detecting repeated spatial patterns in spatial omics data, tissue regions that are spatially separated but share similar underlying cell-type or biological profile distributions. Built on a nonparametric statistical inference framework that integrates constrained clustering with a block-permutation procedure using the maximum mean discrepancy (MMD) statistic, enabling formal hypothesis testing for repeated spatial domains. Works natively with AnnData .h5ad files, the standard format for single-cell and spatial omics workflows.

## Install

```bash
pip install repspat
```

## Usage

```python
import scanpy as sc
from repspat import SampleData, spatial_silhouette_analysis, spatial_constrained_hac
from repspat import create_blocks, multiple_comparison
from repspat import (
    pairwise_results_to_matrix,
    plot_cluster_feature_presence,
    plot_spatial_clusters,
)
```

### Load a sample

```python
adata = sc.read_h5ad("03_TNBC_2018_spe.h5ad")
data = SampleData(
    adata,
    sample_column="sample_id",
    sample_name="Sample_04",
    metric="euclidean",
    cell_type_column="mm"
)
sample = data.sample_adata

data.summary()
# {
#   'feature_mat': (5381, 36),
#   'coords_mat':  (5381, 2),
#   'dist_matrix': (5381, 5381),
#   'cell_type':   (5381, 1),
#   'sample_adata':(5381, 36)
# }
```

If your `.h5ad` already contains thresholded binary features in a layer, pass
that layer name and RepSpat will use it directly:

```python
data = SampleData(
    adata,
    sample_column="sample_id",
    sample_name="Sample_04",
    layer="binary",
    metric="jaccard",
    cell_type_column="mm"
)
sample = data.sample_adata
```

### Find the right number of clusters

Runs spatially-aware silhouette analysis over a range of k and neighbor configs:

```python
results = spatial_silhouette_analysis(sample, n_neighbors_list=[6, 8], n_clusters_range=range(4, 9))
#    n_neighbors  n_clusters  avg_silhouette
# 0            6           4        0.136721
# 5            8           4        0.161823
```

### Cluster cells

```python
labels, sample, model = spatial_constrained_hac(sample, n_clusters=7, n_neighs=8)
# labels are also stored in sample.obs["repspat_region"]
```

### Plot

```python
plot_spatial_clusters(sample)
```

### Plot each cluster's top binary features

Creates one two-panel figure per cluster: a spatial highlight and the cluster's
most prevalent binary features.

```python
cluster_figures = plot_cluster_feature_presence(sample, top_n=5)
```

### Run MMD tests between clusters

```python
create_blocks(sample, knn=8)
results_df = multiple_comparison(sample, kernel="IMQ")

# pairs that are not significantly different
print(results_df[results_df["adj_p"] >= 0.05])
```

### Similarity network

```python
matrix = pairwise_results_to_matrix(sample)
```

See `example.ipynb` for a full walkthrough on a TNBC dataset.

## API

### `SampleData(adata, sample_column=None, sample_name=None, *, metric, layer=None, cell_type_column="mm")`
Loads and optionally subsets an AnnData object for one sample. Computes the feature matrix, spatial coordinates, and pairwise distance matrix using the required `metric` value. Pass `layer` to use an existing AnnData layer, such as a pre-thresholded binary matrix, instead of `adata.X`. Distances are stored in `adata.obsp["repspat_distances"]`.

### `spatial_silhouette_analysis(adata, n_neighbors_list, n_clusters_range)`
Returns a DataFrame of silhouette scores across all `(n_neighbors, n_clusters)` combinations and stores it in `adata.uns["repspat_silhouette"]`. Use this to pick clustering params.

### `spatial_constrained_hac(adata, n_clusters, n_neighs, coord_type, delaunay)`
Ward HAC with a spatial connectivity constraint. Stores labels in `adata.obs["repspat_region"]` and returns `(labels, adata, model)`.

### `create_blocks(adata, knn)`
Splits regions into KMeans blocks for block permutation testing and stores block IDs in `adata.obs["repspat_polygon_id"]`.

### `two_sample_mmd(sample1_idx, sample2_idx, dist_matrix, patient_data, kernel, kernel_param, nperm)`
MMD² between two groups with a permutation null. Returns observed statistic, null distribution, and p-value.

### `multiple_comparison(adata, kernel, kernel_param, nperm, adj_p)`
Pairwise MMD across all cluster pairs. Reads regions, block IDs, and distances from AnnData. Stores results in `adata.uns["repspat_mmd_results"]`. `adj_p` can be `"BH"`, `"bonferroni"`, or `"holm"`.

### `plot_spatial_clusters(adata)`
Scatter plot colored by `adata.obs["repspat_region"]`.

### `plot_cluster_feature_presence(adata, top_n)`
Creates a spatial highlight and top binary-feature prevalence chart for every cluster.

### `pairwise_results_to_matrix(adata, plot)`
Builds an adjacency matrix and network graph from `adata.uns["repspat_mmd_results"]`. Edges connect clusters that are not significantly different.

## Requirements

- Python >= 3.9
- numpy, pandas, scipy, scikit-learn, statsmodels, scanpy, squidpy >= 1.2, matplotlib, networkx

## License

MIT
