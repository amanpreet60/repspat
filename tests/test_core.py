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


def make_layered_adata():
    adata = make_sample_adata()
    adata.X = np.array([[0.2, 4.0], [2.0, 0.1]])
    adata.layers["binary"] = np.array([[0, 1], [1, 0]])
    return adata


def test_sample_data_chooses_metric_from_layer():
    from repspat import SampleData

    euclidean = SampleData(
        make_sample_adata(), "sample_id", "sample", metric="euclidean"
    )
    jaccard = SampleData(
        make_layered_adata(),
        "sample_id",
        "sample",
        metric="jaccard",
        layer="binary",
    )
    explicit = SampleData(
        make_layered_adata(),
        "sample_id",
        "sample",
        metric="cityblock",
        layer="binary",
    )

    assert euclidean.dist_matrix.iloc[0, 1] == pytest.approx(np.sqrt(2))
    assert jaccard.dist_matrix.iloc[0, 1] == pytest.approx(1.0)
    assert explicit.dist_matrix.iloc[0, 1] == pytest.approx(2.0)


def test_sample_data_can_use_feature_layer():
    from repspat import SampleData

    data = SampleData(
        make_layered_adata(),
        "sample_id",
        "sample",
        metric="jaccard",
        layer="binary",
    )

    expected = pd.DataFrame(
        [[0, 1], [1, 0]],
        index=["cell1", "cell2"],
        columns=["marker_a", "marker_b"],
    )
    pd.testing.assert_frame_equal(data.feature_mat, expected)
    np.testing.assert_array_equal(
        data.sample_adata.X,
        np.array([[0.2, 4.0], [2.0, 0.1]]),
    )
    np.testing.assert_array_equal(data.sample_adata.obsp["repspat_distances"], data.dist_matrix.to_numpy())
    assert data.dist_matrix.iloc[0, 1] == pytest.approx(1.0)


def test_sample_data_keeps_original_x_when_using_feature_layer():
    from repspat import SampleData

    data = SampleData(
        make_layered_adata(),
        "sample_id",
        "sample",
        metric="jaccard",
        layer="binary",
    )

    assert data.sample_adata.uns["repspat"]["layer"] == "binary"
    np.testing.assert_array_equal(
        data.feature_mat.to_numpy(),
        np.array([[0, 1], [1, 0]]),
    )
    np.testing.assert_array_equal(
        data.sample_adata.X,
        np.array([[0.2, 4.0], [2.0, 0.1]]),
    )


def make_three_cell_binary_sample():
    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(
        X=np.array(
            [
                [0.2, 4.0],
                [2.0, 0.1],
                [3.0, 5.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "sample_id": ["sample", "sample", "sample"],
                "mm": ["a", "b", "c"],
            },
            index=["cell1", "cell2", "cell3"],
        ),
    )
    adata.var_names = ["marker_a", "marker_b"]
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    adata.layers["binary"] = np.array(
        [
            [1, 0],
            [1, 1],
            [0, 1],
        ]
    )
    return adata


def fake_spatial_neighbors(adata, **kwargs):
    from scipy.sparse import csr_matrix

    n_obs = adata.n_obs
    adata.obsp["spatial_connectivities"] = csr_matrix(
        np.ones((n_obs, n_obs)) - np.eye(n_obs)
    )


def test_spatial_constrained_hac_uses_selected_feature_layer_with_ward(monkeypatch):
    from repspat import SampleData
    import repspat.clustering as clustering

    data = SampleData(
        make_three_cell_binary_sample(),
        "sample_id",
        "sample",
        metric="jaccard",
        layer="binary",
    )
    captured = {}

    class FakeAgglomerativeClustering:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def fit_predict(self, X):
            captured["X"] = np.asarray(X)
            return np.array([0, 1, 1])

    monkeypatch.setattr(clustering, "_spatial_neighbors", fake_spatial_neighbors)
    monkeypatch.setattr(clustering, "AgglomerativeClustering", FakeAgglomerativeClustering)

    labels, clustered_adata, model = clustering.spatial_constrained_hac(
        data.sample_adata,
        n_clusters=2,
        n_neighs=1,
        linkage="ward",
    )

    np.testing.assert_array_equal(
        captured["X"],
        np.array(
            [
                [1, 0],
                [1, 1],
                [0, 1],
            ]
        ),
    )
    np.testing.assert_array_equal(labels, np.array([1, 2, 2]))
    assert clustered_adata.uns["repspat_hac"]["layer"] == "binary"
    assert clustered_adata.uns["repspat_hac"]["metric"] == "jaccard"
    assert clustered_adata.uns["repspat_hac"]["linkage"] == "ward"
    assert model.kwargs["n_clusters"] == 2
    assert model.kwargs["metric"] == "euclidean"
    assert model.kwargs["linkage"] == "ward"


def test_validate_hac_linkage_rejects_unknown_linkage():
    import repspat.clustering as clustering

    with pytest.raises(ValueError, match="linkage must be one of"):
        clustering._validate_hac_linkage("auto")


def test_sample_data_converts_spatial_dataframe_to_numpy():
    from repspat import SampleData

    adata = make_sample_adata()
    adata.obsm["spatial"] = pd.DataFrame(
        {"centroidX": [0.0, 1.0], "centroidY": [0.0, 1.0]},
        index=adata.obs_names,
    )

    data = SampleData(adata, "sample_id", "sample", metric="euclidean")

    assert isinstance(data.sample_adata.obsm["spatial"], np.ndarray)
    np.testing.assert_array_equal(
        data.sample_adata.obsm["spatial"],
        np.array([[0.0, 0.0], [1.0, 1.0]]),
    )


def test_sample_data_rejects_missing_feature_layer():
    from repspat import SampleData

    with pytest.raises(KeyError, match="Layer 'missing' not found"):
        SampleData(
            make_layered_adata(),
            "sample_id",
            "sample",
            metric="jaccard",
            layer="missing",
        )


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

    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(
        X=np.array(
            [
                [1, 0, 0, 0.2],
                [1, 1, 0, 0.4],
                [0, 0, 1, 0.6],
                [0, 1, 1, 0.8],
                [1, 1, 1, 1.0],
            ]
        ),
        obs=pd.DataFrame({"repspat_region": [1, 1, 2, 2, 2]}),
    )
    adata.var_names = ["marker_a", "marker_b", "marker_c", "continuous"]
    adata.obsm["spatial"] = np.array(
        [[0, 0], [1, 1], [2, 0], [3, 1], [4, 0]]
    )

    figures = plot_cluster_feature_presence(
        adata,
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


def test_plot_cluster_feature_presence_falls_back_to_enrichment_for_continuous_data():
    import matplotlib.pyplot as plt

    from repspat import plot_cluster_feature_presence

    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(
        X=np.array(
            [
                [1000.0, 2.0, 4.0],
                [1000.0, 2.0, 5.0],
                [1100.0, 0.0, 5.0],
                [900.0, 0.2, 4.0],
                [1000.0, 0.4, 4.0],
            ]
        ),
        obs=pd.DataFrame({"repspat_region": [1, 1, 2, 2, 2]}),
    )
    adata.var_names = ["large_scale_marker", "cluster_marker", "other_marker"]
    adata.obsm["spatial"] = np.array(
        [[0, 0], [1, 1], [2, 0], [3, 1], [4, 0]]
    )

    figures = plot_cluster_feature_presence(
        adata,
        top_n=2,
    )

    _, cluster_1_axes = figures[1]
    feature_ax = cluster_1_axes[1]
    assert feature_ax.get_xlabel() == "Cluster-vs-rest standardized enrichment"
    assert [tick.get_text() for tick in feature_ax.get_yticklabels()][-1] == "cluster_marker"
    assert feature_ax.patches[-1].get_width() > 0

    for figure, _ in figures.values():
        plt.close(figure)


def test_plot_spatial_clusters_reads_anndata_slots():
    import matplotlib.pyplot as plt

    from repspat import plot_spatial_clusters

    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame({"repspat_region": [1, 1, 2]}),
    )
    adata.obsm["spatial"] = np.array([[0, 0], [1, 1], [2, 0]])

    fig, ax = plot_spatial_clusters(adata)

    assert len(ax.collections) == 2
    plt.close(fig)


def test_plot_spatial_clusters_accepts_dataframe_spatial_coords():
    import matplotlib.pyplot as plt

    from repspat import plot_spatial_clusters

    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(
        X=np.ones((3, 2)),
        obs=pd.DataFrame({"repspat_region": [1, 1, 2]}),
    )
    adata.obsm["spatial"] = pd.DataFrame(
        {"centroidX": [0, 1, 2], "centroidY": [0, 1, 0]},
        index=adata.obs_names,
    )

    fig, ax = plot_spatial_clusters(adata)

    assert len(ax.collections) == 2
    plt.close(fig)


def test_create_blocks_writes_polygon_ids_to_adata():
    from repspat import create_blocks

    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(
        X=np.array([[0, 0], [0, 1], [5, 5], [5, 6]]),
        obs=pd.DataFrame({"repspat_region": [1, 1, 2, 2]}),
    )
    adata.var_names = ["marker_a", "marker_b"]

    blocks = create_blocks(adata, knn=8)

    assert "repspat_polygon_id" in adata.obs
    assert blocks["repspat_polygon_id"].tolist() == [1, 1, 1, 1]


def test_pairwise_results_to_matrix_reads_anndata_results():
    from repspat.visualization import pairwise_results_to_matrix

    anndata = pytest.importorskip("anndata")
    adata = anndata.AnnData(X=np.ones((2, 1)))
    adata.uns["repspat_mmd_results"] = pd.DataFrame(
        {
            "region_1": [1],
            "region_2": [2],
            "obs_mmd_sq": [0.25],
            "adj_p": [0.10],
        }
    )

    matrix = pairwise_results_to_matrix(adata, plot=False)

    assert matrix.loc[1, 2] == pytest.approx(0.25)
