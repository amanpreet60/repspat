import numpy as np
import pandas as pd
import pytest


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


def test_cluster_feature_presence_ranks_thresholded_features():
    from repspat import cluster_feature_presence

    feature_df = pd.DataFrame(
        {
            "marker_a": [1, 1, 0, 0, 1],
            "marker_b": [0, 1, 0, 1, 1],
            "marker_c": [0, 0, 1, 1, 1],
            "region": [1, 1, 2, 2, 2],
        },
        index=["cell1", "cell2", "cell3", "cell4", "cell5"],
    )

    result = cluster_feature_presence(feature_df, top_n=2)

    cluster_1 = result[result["cluster"] == 1].reset_index(drop=True)
    cluster_2 = result[result["cluster"] == 2].reset_index(drop=True)

    assert cluster_1.loc[0, "feature"] == "marker_a"
    assert cluster_1.loc[0, "presence_rate"] == pytest.approx(1.0)
    assert cluster_1.loc[0, "present_count"] == 2
    assert cluster_1.loc[0, "cluster_size"] == 2

    assert cluster_2.loc[0, "feature"] == "marker_c"
    assert cluster_2.loc[0, "presence_rate"] == pytest.approx(1.0)
    assert len(cluster_2) == 2


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
