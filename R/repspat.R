# repspat_wrappers.R

suppressPackageStartupMessages({
  library(reticulate)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(S4Vectors)
})

# -----------------------------
# Python setup
# -----------------------------

init_repspat <- function(python = NULL) {
  if (!is.null(python)) {
    use_python(python, required = TRUE)
  }
  
  repspat <- import("repspat")
  return(repspat)
}

default_repspat_thresholds <- function() {
  c(
    betaCatenin = 0.253063192,
    CD11b = 1.100536222,
    CD11c = 0.213556566,
    CD138 = 1.038502303,
    CD16 = 0.304460727,
    CD20 = 0.296460556,
    CD209 = 0.406669859,
    CD3 = 0.250762424,
    CD31 = 1.006381212,
    CD4 = 0.223762,
    CD45 = 0.327450495,
    CD45RO = 0.236334141,
    CD56 = 0.770183826,
    CD63 = 0.268126101,
    CD68 = 1.128510101,
    CD8 = 0.260237333,
    dsDNA = 0.20301697,
    EGFR = 1.069264889,
    FoxP3 = 0.34649899,
    H3K27me3 = 0.201893232,
    H3K9ac = 0.406071229,
    HLA_Class_1 = 0.23773137,
    HLADR = 0.4408116,
    IDO = 0.312824202,
    Keratin17 = 0.291843697,
    Keratin6 = 0.282480556,
    Ki67 = 0.210413737,
    Lag3 = 1.239097091,
    MPO = 0.324934545,
    p53 = 0.265704293,
    panKeratin = 0.244129571,
    PD1 = 0.972859909,
    PDL1 = 1.352024101,
    pS6 = 0.282699091,
    SMA = 0.250358024,
    Vimentin = 0.28239596
  )
}

as_python_thresholds <- function(thresholds) {
  if (is.null(thresholds)) {
    return(NULL)
  }

  if (is.atomic(thresholds)) {
    thresholds <- as.list(thresholds)
  }

  if (!is.list(thresholds) || is.null(names(thresholds)) || any(names(thresholds) == "")) {
    stop("thresholds must be NULL or a named numeric vector/list.")
  }

  r_to_py(thresholds)
}

as_r_object <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    return(py_to_r(x))
  }

  x
}

feature_matrix_for_anndata <- function(repspat_data) {
  feature_mat <- as_r_object(repspat_data$feature_mat)
  var_names <- as_r_object(repspat_data$sample_adata$var_names$to_list())
  missing_features <- setdiff(var_names, colnames(feature_mat))

  if (length(missing_features) > 0) {
    stop(
      "feature_mat is missing AnnData features: ",
      paste(missing_features, collapse = ", ")
    )
  }

  as.matrix(feature_mat[, var_names, drop = FALSE])
}

sync_anndata_x_to_feature_mat <- function(repspat_data) {
  np <- import("numpy")
  repspat_data$sample_adata$X <- np$array(feature_matrix_for_anndata(repspat_data))
  invisible(repspat_data)
}


# -----------------------------
# Direct AnnData / h5ad workflow
# -----------------------------

load_repspat_sample <- function(
    adata_path,
    sample_column,
    sample_name,
    cell_type_column = "mm",
    metric = "euclidean",
    thresholds = default_repspat_thresholds(),
    python = NULL
) {
  repspat <- init_repspat(python)
  
  data <- repspat$SampleData(
    sample_column = sample_column,
    sample_name = sample_name,
    adata_path = adata_path,
    metric = metric,
    thresholds = as_python_thresholds(thresholds),
    cell_type_column = cell_type_column
  )
  
  return(data)
}


run_repspat_clustering <- function(
    repspat_data,
    n_clusters = 7,
    n_neighs = 8,
    coord_type = "generic",
    delaunay = FALSE,
    cluster_on_feature_mat = TRUE,
    python = NULL
) {
  repspat <- init_repspat(python)

  if (isTRUE(cluster_on_feature_mat)) {
    sync_anndata_x_to_feature_mat(repspat_data)
  }
  
  result <- repspat$spatial_constrained_hac(
    repspat_data$sample_adata,
    repspat_data$feature_mat,
    as.integer(n_clusters),
    as.integer(n_neighs),
    coord_type,
    delaunay
  )
  
  labels <- py_to_r(result[[1]])
  feature_df <- py_to_r(result[[2]])
  model <- result[[3]]
  
  list(
    labels = labels,
    feature_df = feature_df,
    model = model
  )
}


run_repspat_mmd <- function(
    feature_df,
    dist_matrix,
    num_features,
    knn = 8,
    kernel = "IMQ",
    kernel_param = 1.0,
    nperm = 200,
    adj_p = "BH",
    python = NULL
) {
  repspat <- init_repspat(python)
  
  blocked_data <- repspat$create_blocks(
    r_to_py(feature_df),
    as.integer(num_features),
    as.integer(knn)
  )
  
  results <- repspat$multiple_comparison(
    blocked_data,
    dist_matrix,
    kernel,
    kernel_param,
    as.integer(nperm),
    adj_p
  )
  
  list(
    blocked_data = py_to_r(blocked_data),
    mmd_results = py_to_r(results)
  )
}


run_repspat_workflow_from_h5ad <- function(
    adata_path,
    sample_column,
    sample_name,
    cell_type_column = "mm",
    n_clusters = 7,
    n_neighs = 8,
    knn = 8,
    kernel = "IMQ",
    kernel_param = 1.0,
    nperm = 200,
    adj_p = "BH",
    thresholds = default_repspat_thresholds(),
    cluster_on_feature_mat = TRUE,
    python = NULL
) {
  data <- load_repspat_sample(
    adata_path = adata_path,
    sample_column = sample_column,
    sample_name = sample_name,
    cell_type_column = cell_type_column,
    thresholds = thresholds,
    python = python
  )
  
  clustering <- run_repspat_clustering(
    repspat_data = data,
    n_clusters = n_clusters,
    n_neighs = n_neighs,
    cluster_on_feature_mat = cluster_on_feature_mat,
    python = python
  )
  
  num_features <- ncol(clustering$feature_df) - 1
  
  mmd <- run_repspat_mmd(
    feature_df = clustering$feature_df,
    dist_matrix = data$dist_matrix,
    num_features = num_features,
    knn = knn,
    kernel = kernel,
    kernel_param = kernel_param,
    nperm = nperm,
    adj_p = adj_p,
    python = python
  )
  
  list(
    data = data,
    labels = clustering$labels,
    feature_df = clustering$feature_df,
    blocked_data = mmd$blocked_data,
    mmd_results = mmd$mmd_results
  )
}


# -----------------------------
# SpatialExperiment integration
# -----------------------------

spe_to_anndata <- function(
    spe,
    assay_name = NULL,
    spatial_columns = NULL,
    sample_column = "sample_id",
    sample_name = "sample_1"
) {
  ad <- import("anndata")
  pd <- import("pandas")
  np <- import("numpy")
  
  if (is.null(assay_name)) {
    assay_name <- assayNames(spe)[1]
  }
  
  mat <- assay(spe, assay_name)
  
  # SingleCellExperiment / SpatialExperiment stores features x cells.
  # AnnData expects cells x features.
  x <- t(as.matrix(mat))
  
  obs <- as.data.frame(colData(spe))
  var <- as.data.frame(rowData(spe))
  
  if (!sample_column %in% colnames(obs)) {
    obs[[sample_column]] <- sample_name
  }
  
  rownames(obs) <- colnames(spe)
  rownames(var) <- rownames(spe)
  
  obs_py <- pd$DataFrame(obs)
  var_py <- pd$DataFrame(var)
  
  obs_py$index <- rownames(obs)
  var_py$index <- rownames(var)
  
  adata <- ad$AnnData(
    X = np$array(x),
    obs = obs_py,
    var = var_py
  )
  
  adata$obs_names <- rownames(obs)
  adata$var_names <- rownames(var)
  
  if (is.null(spatial_columns)) {
    coords <- spatialCoords(spe)
  } else {
    coords <- as.matrix(colData(spe)[, spatial_columns])
  }
  
  adata$obsm["spatial"] <- np$array(as.matrix(coords))
  
  return(adata)
}


run_repspat_workflow_from_spe <- function(
    spe,
    assay_name = NULL,
    spatial_columns = NULL,
    sample_column = "sample_id",
    sample_name = "sample_1",
    cell_type_column = "cell_type",
    n_clusters = 7,
    n_neighs = 8,
    knn = 8,
    kernel = "IMQ",
    kernel_param = 1.0,
    nperm = 200,
    adj_p = "BH",
    thresholds = default_repspat_thresholds(),
    cluster_on_feature_mat = TRUE,
    python = NULL
) {
  repspat <- init_repspat(python)
  
  adata <- spe_to_anndata(
    spe = spe,
    assay_name = assay_name,
    spatial_columns = spatial_columns,
    sample_column = sample_column,
    sample_name = sample_name
  )
  
  data <- repspat$SampleData(
    sample_column = sample_column,
    sample_name = sample_name,
    adata_obj = adata,
    thresholds = as_python_thresholds(thresholds),
    cell_type_column = cell_type_column
  )
  
  clustering <- run_repspat_clustering(
    repspat_data = data,
    n_clusters = n_clusters,
    n_neighs = n_neighs,
    cluster_on_feature_mat = cluster_on_feature_mat,
    python = python
  )
  
  num_features <- ncol(clustering$feature_df) - 1
  
  mmd <- run_repspat_mmd(
    feature_df = clustering$feature_df,
    dist_matrix = data$dist_matrix,
    num_features = num_features,
    knn = knn,
    kernel = kernel,
    kernel_param = kernel_param,
    nperm = nperm,
    adj_p = adj_p,
    python = python
  )
  
  colData(spe)$repspat_region <- clustering$labels
  colData(spe)$repspat_polygon_id <- mmd$blocked_data$polygon_id
  
  metadata(spe)$repspat_mmd <- mmd$mmd_results
  metadata(spe)$repspat_feature_df <- clustering$feature_df
  
  return(spe)
}


# -----------------------------
# Plotting helper
# -----------------------------

plot_repspat_clusters <- function(
    spe,
    region_column = "repspat_region",
    point_size = 1,
    alpha = 0.8
) {
  coords <- spatialCoords(spe)
  labels <- colData(spe)[[region_column]]
  
  plot(
    coords[, 1],
    coords[, 2],
    col = as.integer(as.factor(labels)),
    pch = 16,
    cex = point_size,
    xlab = "X coordinate",
    ylab = "Y coordinate",
    main = "repSpat spatial regions"
  )
}


result <- run_repspat_workflow_from_h5ad(
  adata_path = "data/03_TNBC_2018_spe.h5ad",
  sample_column = "sample_id",
  sample_name = "Sample_04",
  cell_type_column = "mm",
  n_clusters = 7,
  n_neighs = 8,
  knn = 8,
  kernel = "IMQ",
  nperm = 10
)

head(result$mmd_results)

similar_regions <- result$mmd_results[result$mmd_results$adj_p >= 0.05, ]

print(similar_regions)
