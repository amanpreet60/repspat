# Changelog

## [0.1.0] - 2026-04-01

Initial release.

- `SampleData` — load and subset AnnData samples, compute distance matrices
- `spatial_silhouette_analysis` — pick cluster count using spatial silhouette scores
- `spatial_constrained_hac` — spatially-constrained HAC clustering
- `create_blocks` — KMeans block creation for permutation testing
- `two_sample_mmd` — MMD² test with block permutation null
- `multiple_comparison` — pairwise MMD across all cluster pairs with BH/Bonferroni/Holm correction
- `plot_spatial_clusters` — scatter plot colored by cluster
- `pairwise_results_to_matrix` — similarity network from pairwise MMD results
- `to_binary` — binarize continuous markers with user-defined thresholds
