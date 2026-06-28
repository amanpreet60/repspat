def test_package_imports_core_api():
    import repspat

    assert repspat.__version__ == "0.1.0"
    assert "SampleData" in repspat.__all__
    assert "multiple_comparison" in repspat.__all__
    assert "plot_spatial_clusters" not in repspat.__all__


def test_import_does_not_load_visualization():
    import sys

    sys.modules.pop("repspat", None)
    sys.modules.pop("repspat.visualization", None)

    import repspat  # noqa: F401

    assert "repspat.visualization" not in sys.modules
