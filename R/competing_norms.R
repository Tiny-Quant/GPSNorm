#' Apply upper-quartile (Q3) normalization to logcounts
#'
#' Computes Q3 size factors from raw counts and applies them as a
#' column-wise offset to a log2(count + 1) assay. Returns the same
#' SpatialExperiment with normalized logcounts.
#'
#' @param spe A SpatialExperiment with a counts assay and a logcounts assay.
#' @param counts_assay Character. Name of the raw counts assay (default "counts").
#' @param logcounts_assay Character. Name of the logcounts assay to normalize
#'   (default "logcounts").
#' @param zero.rm Logical. Whether to exclude zeros when computing Q3
#'   (recommended TRUE).
#'
#' @return A SpatialExperiment with Q3-normalized logcounts.
#' @export
apply_q3_logcounts <- function(
  spe,
  counts_assay = "counts",
  logcounts_assay = "logcounts",
  zero.rm = TRUE
) {
  for (pkg in c("SummarizedExperiment", "SpatialExperiment", "Matrix")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Missing package: ", pkg, call. = FALSE)
    }
  }

  if (!methods::is(spe, "SpatialExperiment")) {
    stop("`spe` must be a SpatialExperiment.", call. = FALSE)
  }

  if (!counts_assay %in% SummarizedExperiment::assayNames(spe)) {
    stop("Counts assay '", counts_assay, "' not found in spe.", call. = FALSE)
  }

  if (!logcounts_assay %in% SummarizedExperiment::assayNames(spe)) {
    stop("Logcounts assay '", logcounts_assay, "' not found in spe.", call. = FALSE)
  }

  counts <- SummarizedExperiment::assay(spe, counts_assay)
  logcounts <- SummarizedExperiment::assay(spe, logcounts_assay)

  if (any(counts < 0)) {
    stop("Counts must be non-negative.", call. = FALSE)
  }
  if (any(!is.finite(logcounts))) {
    stop("logcounts contains non-finite values.", call. = FALSE)
  }

  # ---- Compute Q3 size factors from counts ----
  q3 <- apply(counts, 2, function(x) {
    if (zero.rm) {
      x <- x[x > 0]
    }
    if (length(x) == 0L) return(NA_real_)
    stats::quantile(x, 0.75, names = FALSE, type = 7)
  })

  q3_med <- stats::median(q3, na.rm = TRUE)
  if (!is.finite(q3_med) || q3_med <= 0) {
    stop("Median Q3 is non-finite or zero; cannot normalize.", call. = FALSE)
  }

  size_factors <- q3 / q3_med
  size_factors[!is.finite(size_factors) | size_factors <= 0] <- 1

  # ---- Apply Q3 normalization on log2 scale ----
  # log2(count + 1) - log2(size_factor)
  logcounts_q3 <- sweep(logcounts, 2, log2(size_factors), FUN = "-")

  logcounts(spe) <- logcounts_q3

  spe
}


# Below taken from SpaNorm\SpaNorm_files\SpaNorm_files\codes\spatial_clustering_wrappers.R
#' Normalize with scran deconvolution size factors
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment` with `counts`.
#'
#' @return Object with `logcounts` added/updated using scran normalization.
#' @export
normalise_scran <- function(spe) {
  cl = scran::quickCluster(spe)
  spe = scran::computeSumFactors(spe, clusters = cl)
  spe$sizeFactor <- pmax(1e-08,spe$sizeFactor)
  spe = scater::logNormCounts(spe)

  return(spe)
}


#' Normalize with SCTransform (`sctransform::vst`)
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment` with `counts`.
#'
#' @return Object subset to modeled genes with `logcounts` set to SCT output.
#' @export
normalise_sct <- function(spe) {
  vst_res = sctransform::vst(SingleCellExperiment::counts(spe, vst.flavor = 'v2'),
                             verbosity = 0)
  #subset data
  spe = spe[rownames(vst_res$y), ]
  SingleCellExperiment::logcounts(spe) = vst_res$y

  return(spe)
}


#' Normalize with RUV-III-NB using NEG controls
#'
#' @param spe A `SpatialExperiment` containing `rowData(spe)$gene_type` with
#'   `"NEG"` controls.
#' @param K Integer number of unwanted factors for RUV-III.
#'
#' @return Object with RUV-adjusted values written to `logcounts`.
#' @export
normalise_ruv3nb <- function(spe, K = 1) {

  if (!requireNamespace("ruvIIInb", quietly = TRUE)) {
    stop("Package ruvIIInb required.")
  }
  if (!requireNamespace("tripack", quietly = TRUE)) {
    stop("Package tripack required.")
  }

  # --- Use simulated NEG controls ---
  if (!"gene_type" %in% colnames(rowData(spe))) {
    stop("gene_type column missing in rowData(spe).")
  }

  ctl <- rowData(spe)$gene_type == "NEG"
  if (sum(ctl, na.rm = TRUE) < 10) {
    stop("Not enough NEG control genes for RUV-III.")
  }

  # --- Build replicate matrix M ---
  coords <- spatialCoords(spe)
  tessel <- tripack::tri.mesh(x = coords[,1], y = coords[,2])
  nb <- tripack::neighbours(tessel)

  n.grid <- round(sqrt(0.1 * ncol(spe) / 20))
  x.pos <- quantile(coords[,1], prob = seq(0.2, 0.8, len = n.grid))
  y.pos <- quantile(coords[,2], prob = seq(0.2, 0.8, len = n.grid))

  seed.pos <- integer(0)
  for (x in x.pos) {
    for (y in y.pos) {
      seed.pos <- c(seed.pos,
        which.min((coords[,1] - x)^2 + (coords[,2] - y)^2)
      )
    }
  }

  M <- matrix(0, nrow(coords), length(seed.pos))
  for (j in seq_along(seed.pos)) {
    idx <- nb[[seed.pos[j]]]
    tmp <- unique(unlist(nb[idx]))
    tmp <- tmp[!rowSums(M) > 0]
    M[tmp, j] <- 1
  }

  # --- Run RUV-III NB ---
  out <- ruvIIInb::fastruvIII.nb(
    SingleCellExperiment::counts(spe),
    k = K,
    ctl = ctl,
    M = M,
    ncores = 1
  )

  res <- ruvIIInb::get.res(out, type = "pearson", batch = out$batch)

  spe <- spe[rownames(out$counts), ]
  SingleCellExperiment::logcounts(spe, withDimnames = FALSE) <- res

  spe
}


#' Normalize using Giotto standard normalization
#'
#' @param spe A `SpatialExperiment` with counts and spatial coordinates.
#'
#' @return `spe` with `logcounts` populated from Giotto normalized expression.
#' @export
normalise_Giotto <- function(spe) {
    library(Giotto)
    library(RcppZiggurat)
    #  giotto <- createGiottoObject(countMat=counts(spe),
    #                           spatial_locs=spatialCoords(spe),
    #                           cell_metadata=colData(spe), filter = FALSE,
    #                           normalize = TRUE )

    giotto <- createGiottoObject(
        raw_exprs=counts(spe),
        spatial_locs=spatialCoords(spe), 
        cell_metadata=colData(spe)
    )

    giotto <- normalizeGiotto(
        giotto, norm_methods = "standard"
    )

    assays(spe,withDimnames=FALSE)$logcounts <- giotto@norm_expr
    return(spe)
}


#' Normalize by library-size scaling then log-transform
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment` containing `counts`.
#'
#' @return Object with scaled `logcounts` assay.
#' @export
normalise_ls <- function(spe) {

  stopifnot("counts" %in% SummarizedExperiment::assayNames(spe))

  counts_mat <- SummarizedExperiment::assay(spe, "counts")

  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- Matrix::Matrix(counts_mat, sparse = TRUE)
  }

  # Library sizes
  libsize <- Matrix::colSums(counts_mat)
  sf <- median(libsize) / libsize

  SingleCellExperiment::sizeFactors(spe) <- sf

  # Scale columns (NO base t())
  emat <- counts_mat %*% Matrix::Diagonal(x = sf)

  log_mat <- log2(emat + 1)

  SummarizedExperiment::assay(spe, "logcounts", withDimnames = FALSE) <- log_mat

  spe
}
