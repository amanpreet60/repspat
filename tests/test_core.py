import numpy as np
import pandas as pd
import pytest


def make_sample_adata():
    anndata = pytest.importorskip("anndata")

    adata = anndata.AnnData(
        X=np.array([[0.0, 0.0], [1.0, 1.0]]),
        obs=pd.DataFrame(
            {"sample_id": ["sample", "sample"], "mm": ["a", "b"]},
            index=["cell1", "cell2"],
        ),
    )
    adata.var_names = ["marker_a", "marker_b"]
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0]])
    return adata


def test_sample_data_chooses_metric_from_thresholds():
    from repspat import SampleData

    euclidean = SampleData(
        "sample_id", "sample", adata_obj=make_sample_adata()
    )
    jaccard = SampleData(
        "sample_id",
        "sample",
        adata_obj=make_sample_adata(),
        thresholds={"marker_a": 0.5, "marker_b": 0.5},
    )
    explicit = SampleData(
        "sample_id",
        "sample",
        adata_obj=make_sample_adata(),
        thresholds={"marker_a": 0.5, "marker_b": 0.5},
        metric="cityblock",
    )

    assert euclidean.dist_matrix.iloc[0, 1] == pytest.approx(np.sqrt(2))
    assert jaccard.dist_matrix.iloc[0, 1] == pytest.approx(1.0)
    assert explicit.dist_matrix.iloc[0, 1] == pytest.approx(2.0)


def test_to_binary_warns_when_threshold_missing():
    from repspat import to_binary

    column = pd.Series([0.1, 0.7, 1.2], name="marker_a")

    with pytest.warns(UserWarning, match="No threshold defined"):
        result = to_binary(column, "marker_a", thresholds={})

    pd.testing.assert_series_equal(result, column)


def test_to_binary_applies_threshold():
    from repspat import to_binary

    column = pd.Series([0.1, 0.7, 1.2], name="marker_a")
    result = to_binary(column, "marker_a", thresholds={"marker_a": 0.7})

    pd.testing.assert_series_equal(result, pd.Series([0, 1, 1], name="marker_a"))


def test_compute_mmd_sq_df_gaussian_identical_groups_is_zero():
    from repspat.mmd import compute_mmd_sq_df

    dist = pd.DataFrame(
        [[0.0, 1.0], [1.0, 0.0]],
        index=["a", "b"],
        columns=["a", "b"],
    )

    result = compute_mmd_sq_df(["a", "b"], ["a", "b"], dist, kernel="Gaussian")

    assert result == pytest.approx(0.0)


def test_compute_mmd_sq_df_rejects_unknown_kernel():
    from repspat.mmd import compute_mmd_sq_df

    dist = pd.DataFrame([[0.0]], index=["a"], columns=["a"])

    with pytest.raises(ValueError, match="kernel must be one of"):
        compute_mmd_sq_df(["a"], ["a"], dist, kernel="linear")


def test_custom_silhouette_uses_zero_for_clusters_without_neighbors():
    from repspat import custom_silhouette

    clusters = pd.Series([0, 0, 1], index=["a", "b", "c"])
    dist = pd.DataFrame(
        [[0.0, 1.0, 5.0], [1.0, 0.0, 5.0], [5.0, 5.0, 0.0]],
        index=clusters.index,
        columns=clusters.index,
    )
    adjacency = pd.DataFrame(0, index=clusters.index, columns=clusters.index)

    result = custom_silhouette(clusters, dist, adjacency)

    np.testing.assert_array_equal(result, np.zeros(3))


def test_pairwise_results_to_matrix_links_non_significant_pairs():
    from repspat.visualization import pairwise_results_to_matrix

    results = pd.DataFrame(
        {
            "region_1": [1, 1, 2],
            "region_2": [2, 3, 3],
            "obs_mmd_sq": [0.25, 0.75, 0.5],
            "adj_p": [0.10, 0.01, 0.20],
        }
    )

    matrix = pairwise_results_to_matrix(results, plot=False)

    assert matrix.loc[1, 2] == pytest.approx(0.25)
    assert matrix.loc[2, 3] == pytest.approx(0.5)
    assert matrix.loc[1, 3] == pytest.approx(0.0)


def test_plot_cluster_feature_presence_creates_ranked_plot_per_cluster():
    import matplotlib.pyplot as plt

    from repspat import plot_cluster_feature_presence

    feature_df = pd.DataFrame(
        {
            "marker_a": [1, 1, 0, 0, 1],
            "marker_b": [0, 1, 0, 1, 1],
            "marker_c": [0, 0, 1, 1, 1],
            "continuous": [0.2, 0.4, 0.6, 0.8, 1.0],
        }
    )

    figures = plot_cluster_feature_presence(
        x=[0, 1, 2, 3, 4],
        y=[0, 1, 0, 1, 0],
        labels=[1, 1, 2, 2, 2],
        feature_df=feature_df,
        top_n=2,
    )

    assert set(figures) == {1, 2}

    cluster_1_figure, cluster_1_axes = figures[1]
    spatial_ax, feature_ax = cluster_1_axes
    assert len(spatial_ax.collections[1].get_offsets()) == 2
    assert spatial_ax.get_legend() is None
    assert [tick.get_text() for tick in feature_ax.get_yticklabels()][-1] == "marker_a"
    assert feature_ax.patches[-1].get_width() == pytest.approx(100.0)

    for figure, _ in figures.values():
        plt.close(figure)
