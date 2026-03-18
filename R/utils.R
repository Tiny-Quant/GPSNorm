# ------------------------------------------------------------------------------
# Data structure adapters for GPSNorm
#   - long-form GeoMx tibble  <->  SpatialExperiment
#   - validation helpers for unit tests
#
# Conventions
#   - assays(spe)$counts: integer matrix [genes x AOIs]
#   - spatialCoords(spe): numeric matrix [AOIs x 2] with colnames c("x","y")
#   - rowData(spe)$gene_type: factor with levels c("NEG","HK","REG")
#   - colData(spe)$sample_id: slide/patient identifier (character)
#   - colData(spe)$Mix: binary group (0/1 numeric or integer)
# ------------------------------------------------------------------------------

#' Validate required columns in a GeoMx long-form tibble
#'
#' @param data A data.frame/tibble in long format.
#' @param required Character vector of required column names.
#'
#' @return Invisibly returns TRUE on success; errors otherwise.
#' @keywords internal
validate_geomx_long <- function(data, required) {
    if (!is.data.frame(data)) {
        stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!is.character(required) || anyNA(required)) {
        stop("`required` must be a character vector without NA.", call. = FALSE)
    }
    miss <- setdiff(required, names(data))
    if (length(miss) > 0L) {
        stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

#' Validate a SpatialExperiment contract used by GPSNorm
#'
#' @param spe A SpatialExperiment.
#' @param assay Name of the assay that must exist (default "counts").
#'
#' @return Invisibly TRUE on success.
#' @keywords internal
validate_spe_contract <- function(spe, assay = "counts") {
    if (!inherits(spe, "SpatialExperiment")) {
        stop("`spe` must be a SpatialExperiment.", call. = FALSE)
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Missing package: SummarizedExperiment.", call. = FALSE)
    }
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
        stop("Missing package: SpatialExperiment.", call. = FALSE)
    }

    a <- SummarizedExperiment::assayNames(spe)
    if (!assay %in% a) {
        stop("Assay '", assay, "' not found. Available: ", paste(a, collapse = ", "), call. = FALSE)
    }

    # spatial coords
    sc <- SpatialExperiment::spatialCoords(spe)
    if (!is.matrix(sc) || ncol(sc) != 2L) {
        stop("`spatialCoords(spe)` must be a numeric matrix with 2 columns.", call. = FALSE)
    }

    # rowData gene_type (recommended but we enforce if present)
    #   rd <- SummarizedExperiment::rowData(spe)
    #   if (!"gene_type" %in% colnames(rd)) {
    #     stop("`rowData(spe)$gene_type` is required (NEG/HK/REG).", call. = FALSE)
    #   }
    #   gt <- as.character(rd$gene_type)
    #   if (!all(gt %in% c("NEG", "HK", "REG"))) {
    #     stop("`rowData(spe)$gene_type` must be in {NEG, HK, REG}.", call. = FALSE)
    #   }

    invisible(TRUE)
}

#' Load a GeoMx DSP export for INLA-style modeling
#'
#' Reads a NanoString GeoMx DSP Excel export (xlsx) containing the sheets
#' `"SegmentProperties"`, `"TargetProperties"`, and `"TargetCountMatrix"`, and
#' returns:
#' (i) an AOI-level table, (ii) a target/gene-level table, and (iii) a long-form
#' AOI-by-gene count table with numeric IDs and common derived fields.
#'
#' This function is intended as a preprocessing step for downstream models
#' (e.g., INLA), but does not require INLA itself.
#'
#' @param path_xlsx Character scalar. Path to the GeoMx Excel file.
#' @param seg_covariates Character vector. Names of segment/AOI-level covariates
#'   to carry forward from the `SegmentProperties` sheet. These columns must be
#'   present in the sheet. Defaults to common sequencing/QC covariates.
#'
#' @return Long-form AOI-by-gene data (`tibble`) including counts, IDs,
#'   coordinates, requested covariates, and derived indicators.
#'
#' @details
#' Required columns:
#' \itemize{
#'   \item `SegmentProperties`: `SlideName`, `ROILabel`, `SegmentDisplayName`,
#'   `ROICoordinateX`, `ROICoordinateY`, plus any `seg_covariates`.
#'   \item `TargetProperties`: `TargetName`, `CodeClass`.
#'   \item `TargetCountMatrix`: a `TargetName` column; remaining columns
#'   correspond to `SegmentDisplayName`.
#' }
#'
#' `CodeClass_clean` is mapped as:
#' \itemize{
#'   \item `"Negative"` \eqn{\rightarrow} `"NEG"`
#'   \item `"Control"` \eqn{\rightarrow} `"HK"` (housekeeping)
#'   \item `"Endogenous"` \eqn{\rightarrow} `"REG"`
#' }
#'
#' @examples
#' \dontrun{
#' dat <- load_geomx_for_inla(
#'     "geomx.xlsx",
#'     seg_covariates = c("AOISurfaceArea", "AOINucleiCount")
#' )
#' head(dat$long)
#' }
#'
#' @export
load_geomx_for_inla <- function(path_xlsx,
                                seg_covariates = c(
                                    "AOISurfaceArea",
                                    "AOINucleiCount",
                                    "DeduplicatedReads",
                                    "SequencingSaturation",
                                    "Mix"
                                )) {
    if (!is.character(path_xlsx) ||
        length(path_xlsx) != 1L ||
        is.na(path_xlsx) ||
        !nzchar(path_xlsx)) {
        stop("`path_xlsx` must be a non-empty character scalar.", call. = FALSE)
    }
    if (!file.exists(path_xlsx)) {
        stop("File not found: ", path_xlsx, call. = FALSE)
    }

    if (is.null(seg_covariates)) {
        seg_covariates <- character(0)
    }
    if (!is.character(seg_covariates)) {
        stop("`seg_covariates` must be a character vector.", call. = FALSE)
    }
    if (anyNA(seg_covariates)) {
        stop("`seg_covariates` cannot contain NA.", call. = FALSE)
    }
    seg_covariates <- unique(seg_covariates)
    seg_covariates <- seg_covariates[nzchar(seg_covariates)]

    for (pkg in c("readxl", "dplyr", "tidyr")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            msg <- paste0(
                "Package '", pkg, "' is required. Install it with ",
                "install.packages('", pkg, "')."
            )
            stop(msg, call. = FALSE)
        }
    }

    seg <- readxl::read_excel(path_xlsx, sheet = "SegmentProperties")
    tgt <- readxl::read_excel(path_xlsx, sheet = "TargetProperties")
    mat <- readxl::read_excel(path_xlsx, sheet = "TargetCountMatrix")

    seg <- dplyr::as_tibble(seg)
    tgt <- dplyr::as_tibble(tgt)
    mat <- dplyr::as_tibble(mat)

    if (!"TargetName" %in% names(mat)) {
        if (ncol(mat) < 1L) {
            stop("`TargetCountMatrix` sheet has no columns.", call. = FALSE)
        }
        names(mat)[1] <- "TargetName"
    }

    base_seg <- c(
        "SlideName",
        "ROILabel",
        "SegmentDisplayName",
        "ROICoordinateX",
        "ROICoordinateY"
    )
    req_seg <- unique(c(base_seg, seg_covariates))
    missing_seg <- setdiff(req_seg, names(seg))
    if (length(missing_seg) > 0L) {
        stop(
            "`SegmentProperties` is missing required columns: ",
            paste(missing_seg, collapse = ", "),
            call. = FALSE
        )
    }

    req_tgt <- c("TargetName", "CodeClass")
    missing_tgt <- setdiff(req_tgt, names(tgt))
    if (length(missing_tgt) > 0L) {
        stop(
            "`TargetProperties` is missing required columns: ",
            paste(missing_tgt, collapse = ", "),
            call. = FALSE
        )
    }

    aoi_keep <- c(
        "AOI_ID",
        "Slide_ID",
        base_seg,
        seg_covariates
    )
    aoi_keep <- unique(aoi_keep)

    aoi_df <-
        seg |>
        dplyr::mutate(
            AOI_ID = as.integer(factor(SegmentDisplayName)),
            Slide_ID = as.integer(factor(SlideName))
        ) |>
        dplyr::select(dplyr::all_of(aoi_keep)) |>
        dplyr::distinct(AOI_ID, .keep_all = TRUE)

    target_df <-
        tgt |>
        dplyr::mutate(
            Gene_ID = as.integer(factor(TargetName)),
            CodeClass_clean = dplyr::case_when(
                CodeClass == "Negative" ~ "NEG",
                CodeClass == "Control" ~ "HK",
                CodeClass == "Endogenous" ~ "REG",
                TRUE ~ as.character(CodeClass)
            ),
            CodeClass_clean = factor(CodeClass_clean,
                levels = c("NEG", "HK", "REG")
            )
        ) |>
        dplyr::select(
            Gene_ID,
            TargetName,
            CodeClass,
            CodeClass_clean
        )

    mat_long <-
        mat |>
        tidyr::pivot_longer(
            cols = -TargetName,
            names_to = "SegmentDisplayName",
            values_to = "Count"
        )

    join_keep <- c(
        "AOI_ID",
        "Slide_ID",
        "SegmentDisplayName",
        base_seg,
        seg_covariates
    )
    join_keep <- unique(join_keep)

    aoi_join <- aoi_df |>
        dplyr::select(dplyr::all_of(join_keep))

    long_df <-
        mat_long |>
        dplyr::left_join(target_df, by = "TargetName") |>
        dplyr::left_join(aoi_join, by = "SegmentDisplayName") |>
        dplyr::mutate(
            aoi_id = AOI_ID,
            gene_id = Gene_ID,
            y_log1p = log1p(Count),
            is_NEG = as.integer(CodeClass_clean == "NEG"),
            is_HK = as.integer(CodeClass_clean == "HK"),
            is_REG = as.integer(CodeClass_clean == "REG"),
            Count_int = as.integer(Count)
        ) |>
        dplyr::arrange(aoi_id, gene_id)

    return(long_df)
}

#' Convert GeoMx long-form tibble to a SpatialExperiment
#'
#' @param data Long-form tibble from `load_geomx_for_inla()`.
#' @param count_col Name of count column to use for assays(counts).
#' @param aoi_col AOI identifier column name.
#' @param gene_col Gene identifier column name.
#' @param gene_name_col Optional gene name column for rownames (e.g., "TargetName").
#' @param gene_type_col Column giving gene type (default "CodeClass_clean": NEG/HK/REG).
#' @param slide_col Column to store in colData as `sample_id` (default "SlideName").
#' @param group_col Optional binary group variable to include in colData (default "Mix").
#' @param coord_cols Length-2 character vector for x/y columns.
#' @param extra_coldata Optional AOI-level columns to carry into colData.
#' @param extra_rowdata Optional gene-level columns to carry into rowData.
#'
#' @return A SpatialExperiment with assays(counts), rowData(gene_type), colData, spatialCoords.
#' @export
geomx_long_to_spe <- function(
    data,
    count_col = "Count_int",
    aoi_col = "AOI_ID",
    gene_col = "Gene_ID",
    gene_name_col = "TargetName",
    gene_type_col = "CodeClass_clean",
    slide_col = "SlideName",
    group_col = "Mix",
    coord_cols = c("ROICoordinateX", "ROICoordinateY"),
    extra_coldata = c(),
    extra_rowdata = c()) {
    for (pkg in c("dplyr", "tidyr", "SpatialExperiment", "SummarizedExperiment", "S4Vectors")) {
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
    }

    stopifnot(is.character(coord_cols), length(coord_cols) == 2L)
    if (is.null(extra_coldata)) extra_coldata <- character(0)
    if (is.null(extra_rowdata)) extra_rowdata <- character(0)

    required <- c(count_col, aoi_col, gene_col, gene_type_col, coord_cols, slide_col)
    if (!is.null(group_col) && nzchar(group_col)) required <- c(required, group_col)

    validate_geomx_long(data, unique(required))

    df <- dplyr::as_tibble(data)

    # ---- AOI-level table ----
    col_keep <- unique(c(aoi_col, slide_col, group_col, coord_cols, extra_coldata))
    col_keep <- col_keep[nzchar(col_keep)]

    coldata <- df |>
        dplyr::select(dplyr::all_of(col_keep)) |>
        dplyr::distinct(.data[[aoi_col]], .keep_all = TRUE) |>
        dplyr::arrange(.data[[aoi_col]])

    # enforce integer AOI IDs for stable ordering
    aoi_ids <- coldata[[aoi_col]]
    aoi_ids <- as.integer(aoi_ids)
    coldata[[aoi_col]] <- aoi_ids

    # sample_id for SpaNorm-style code
    coldata <- coldata |>
        dplyr::mutate(sample_id = as.character(.data[[slide_col]]))

    if (!"Slide_ID" %in% colnames(coldata)) {
        if ("slide_id" %in% colnames(coldata)) {
            coldata$Slide_ID <- as.integer(coldata$slide_id)
        } else {
            coldata$Slide_ID <- 1L
        }
    }

    # Create binary group_col factor.
    if (!is.null(group_col) && nzchar(group_col)) {
        # Coerce to factor
        g <- as.factor(coldata[[group_col]])

        # Enforce binary
        lv <- levels(g)
        if (length(lv) != 2L) {
            stop(
                "group_col must have exactly 2 levels; found: ",
                paste(lv, collapse = ", "),
                call. = FALSE
            )
        }

        # Keep current level order (deterministic)
        g <- factor(g, levels = lv)

        # Binary encoding: first level -> 0, second level -> 1
        g_bin <- as.integer(g) - 1L

        # ---- Print mapping (explicit + unambiguous) ----
        message(
            sprintf(
                "Group encoding: '%s' -> 0 ; '%s' -> 1",
                lv[1], lv[2]
            )
        )

        # Store results
        coldata[[group_col]] <- g
        coldata[[paste0(group_col, "_bin")]] <- g_bin
    }

    # spatial coordinates matrix aligned to coldata rows
    coords <- coldata |>
        dplyr::select(dplyr::all_of(coord_cols)) |>
        dplyr::mutate(
            x = as.numeric(.data[[coord_cols[1]]]),
            y = as.numeric(.data[[coord_cols[2]]])
        ) |>
        dplyr::select(x, y) |>
        as.matrix()
    colnames(coords) <- c("x", "y")

    # ---- gene-level table ----
    row_keep <- unique(c(gene_col, gene_type_col, gene_name_col, extra_rowdata))
    row_keep <- row_keep[nzchar(row_keep)]

    rowdata <- df |>
        dplyr::select(dplyr::all_of(row_keep)) |>
        dplyr::distinct(.data[[gene_col]], .keep_all = TRUE) |>
        dplyr::arrange(.data[[gene_col]])

    gene_ids <- as.integer(rowdata[[gene_col]])
    rowdata[[gene_col]] <- gene_ids

    # gene type
    gene_type <- as.character(rowdata[[gene_type_col]])
    if (!all(gene_type %in% c("NEG", "HK", "REG"))) {
        stop("`", gene_type_col, "` must contain only {NEG, HK, REG}.", call. = FALSE)
    }
    # rogene_type <- factor(gene_type, levels = c("NEG", "HK", "REG"))
    rowdata$gene_type <- factor(gene_type, levels = c("NEG", "HK", "REG"))
    if (gene_type_col != "gene_type") {
        rowdata[[gene_type_col]] <- NULL
    }

    # rownames: prefer gene_name_col if present and non-missing; else Gene_ID
    if (!is.null(gene_name_col) && gene_name_col %in% names(rowdata)) {
        rn <- as.character(rowdata[[gene_name_col]])
        rn[is.na(rn) | rn == ""] <- paste0("gene_", rowdata[[gene_col]][is.na(rn) | rn == ""])
        rownames_use <- make.unique(rn)
    } else {
        rownames_use <- paste0("gene_", rowdata[[gene_col]])
    }

    # ---- counts matrix [genes x AOIs] ----
    mat <- df |>
        dplyr::select(
            dplyr::all_of(c(gene_col, aoi_col, count_col))
        ) |>
        dplyr::mutate(
            Gene_ID = as.integer(.data[[gene_col]]),
            AOI_ID  = as.integer(.data[[aoi_col]]),
            y       = as.integer(.data[[count_col]])
        ) |>
        dplyr::select(Gene_ID, AOI_ID, y) |>
        tidyr::pivot_wider(
            names_from = AOI_ID,
            values_from = y,
            values_fill = 0L
        )

    gene_order <- rowdata[[gene_col]]
    aoi_order <- coldata[[aoi_col]]

    # pivot_wider creates a Gene_ID column + AOI columns
    mat_gene <- mat$Gene_ID
    mat_num <- mat |>
        dplyr::select(-Gene_ID) |>
        as.matrix()

    # reorder columns to match coldata AOI order
    colnames(mat_num) <- as.character(colnames(mat_num))
    want_cols <- as.character(aoi_order)
    if (!all(want_cols %in% colnames(mat_num))) {
        stop("Counts matrix is missing AOIs present in AOI table.", call. = FALSE)
    }
    mat_num <- mat_num[, want_cols, drop = FALSE]

    # reorder rows to match rowdata gene order
    idx <- match(gene_order, mat_gene)
    if (anyNA(idx)) stop("Counts matrix is missing genes present in gene table.", call. = FALSE)
    mat_num <- mat_num[idx, , drop = FALSE]

    # Dim check.
    if (!is.matrix(mat_num) || length(dim(mat_num)) != 2L) {
        stop("Internal error: assay matrix is not 2D.", call. = FALSE)
    }
    storage.mode(mat_num) <- "integer"
    rownames(mat_num) <- rownames_use
    colnames(mat_num) <- paste0("AOI_", aoi_order)

    # ---- build SpatialExperiment ----
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = mat_num),
        rowData = S4Vectors::DataFrame(rowdata),
        colData = S4Vectors::DataFrame(coldata),
        spatialCoords = coords
    )

    spe
}

#' Convert SpatialExperiment to GeoMx long-form tibble
#'
#' @param spe SpatialExperiment.
#' @param assay Name of assay to extract (default "counts").
#' @param group_col Name of group variable in colData to include (default "Mix").
#'
#' @return Long-form tibble with AOI_ID, Gene_ID, Count_int and metadata.
#' @export
spe_to_geomx_long <- function(
    spe,
    assay = "counts",
    group_col = "Mix") {
    for (pkg in c("dplyr", "tidyr", "SpatialExperiment", "SummarizedExperiment")) {
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
    }

    validate_spe_contract(spe, assay = assay)

    mat <- SummarizedExperiment::assay(spe, assay)
    rd <- SummarizedExperiment::rowData(spe)
    cd <- SummarizedExperiment::colData(spe)
    coords <- SpatialExperiment::spatialCoords(spe)

    # Extract if present
    Gene_ID <- if ("Gene_ID" %in% colnames(rd)) as.integer(rd$Gene_ID) else seq_len(nrow(mat))
    TargetName <- if ("gene_name" %in% colnames(rd)) rd$gene_name else seq_len(nrow(mat))
    CodeClass_clean <- if ("gene_type" %in% colnames(rd)) as.factor(rd$gene_type) else rep("REG", nrow(mat))
    AOI_ID <- if ("AOI_ID" %in% colnames(cd)) as.integer(cd$AOI_ID) else seq_len(ncol(mat))

    gene_meta <- tibble::tibble(
        Gene_ID = Gene_ID,
        TargetName = TargetName,
        CodeClass_clean = CodeClass_clean
    )

    # Extract and format counts.
    if (inherits(mat, "dgCMatrix")) {
        long <- Matrix::summary(mat)
        colnames(long) <- c("Gene_ID", "AOI_index", "Count_int")
        long$AOI_col <- colnames(mat)[long$AOI_index]
    } else {
        long <- as.data.frame(mat)
        long$Gene_ID <- Gene_ID
        long <- tidyr::pivot_longer(
            long,
            cols = -Gene_ID,
            names_to = "AOI_col",
            values_to = "Count_int"
        )
    }

    long <- long |>
        dplyr::mutate(
            # map AOI_col back to index
            AOI_index = match(AOI_col, colnames(mat)),
            AOI_ID = AOI_ID[AOI_index]
        ) |>
        dplyr::left_join(gene_meta, by = "Gene_ID") |>
        dplyr::select(AOI_ID, Gene_ID, Count_int, TargetName, CodeClass_clean)

    # attach AOI-level metadata
    cd_df <- as.data.frame(cd)
    cd_df$AOI_ID <- AOI_ID

    # Ensure coordinates only come from spatialCoords()
    coord_cols <- c("ROICoordinateX", "ROICoordinateY")
    cd_df <- cd_df[, setdiff(colnames(cd_df), coord_cols), drop = FALSE]

    coords_df <- data.frame(
        AOI_ID = AOI_ID,
        ROICoordinateX = coords[, 1],
        ROICoordinateY = coords[, 2]
    )

    out <- long |>
        dplyr::left_join(cd_df, by = "AOI_ID") |>
        dplyr::left_join(coords_df, by = "AOI_ID")

    # attach gene-level metadata
    rd_df <- as.data.frame(rd)
    rd_df$Gene_ID <- Gene_ID
    out <- out |>
        dplyr::left_join(rd_df, by = "Gene_ID")

    # keep key columns first
    first <- c("Gene_ID", "AOI_ID", "Count_int", group_col, "sample_id", "ROICoordinateX", "ROICoordinateY", "TargetName", "CodeClass_clean")
    first <- first[first %in% names(out)]
    out <- out |>
        dplyr::select(dplyr::all_of(first), dplyr::everything())

    out <- out |>
        dplyr::mutate(
            is_NEG = as.integer(CodeClass_clean == "NEG"),
            is_HK = as.integer(CodeClass_clean == "HK"),
            is_REG = as.integer(CodeClass_clean == "REG"),
        )

    dplyr::as_tibble(out)
}

#' Split a simulated SpatialExperiment into a list of per-slide objects
#'
#' @param spe A SpatialExperiment produced by simulate_spe_benchmark().
#' @param slide_col Character. Column in colData(spe) holding slide IDs.
#'   Default "slide_id".
#'
#' @return A named list of SpatialExperiment objects, one per slide.
#' @export
split_spe_by_slide <- function(spe, slide_col = "slide_id") {
    if (!inherits(spe, "SpatialExperiment")) {
        stop("`spe` must be a SpatialExperiment.", call. = FALSE)
    }
    cd <- SummarizedExperiment::colData(spe)
    if (!slide_col %in% colnames(cd)) {
        stop("`slide_col` not found in colData(spe).", call. = FALSE)
    }

    sid <- cd[[slide_col]]
    # Preserve ordering, tolerate integer/character slide ids
    u <- unique(as.character(sid))

    out <- lapply(u, function(k) {
        spe[, as.character(sid) == k, drop = FALSE]
    })
    names(out) <- paste0("slide_", u)
    out
}


#' Convert sim_domain to a factor that includes NA as an explicit level
#'
#' SpaNorm-style plotting/evaluation often expects domain labels to be a factor.
#' This helper converts `sim_domain` into a factor and replaces missing labels
#' with a chosen level name (default "Unlabeled").
#'
#' @param spe A SpatialExperiment (e.g., a per-slide split).
#' @param domain_col Character. Domain column in colData(spe).
#'   Default "sim_domain".
#' @param na_level Character. Label used for missing domains.
#'   Default "Unlabeled".
#'
#' @return The same SpatialExperiment with colData(spe)[[domain_col]] updated.
#' @export
domain_na_as_level <- function(
    spe,
    domain_col = "sim_domain",
    na_level = "Unlabeled") {
    if (!inherits(spe, "SpatialExperiment")) {
        stop("`spe` must be a SpatialExperiment.", call. = FALSE)
    }
    cd <- SummarizedExperiment::colData(spe)
    if (!domain_col %in% colnames(cd)) {
        stop("`domain_col` not found in colData(spe).", call. = FALSE)
    }

    d <- cd[[domain_col]]

    # Convert to character so we can replace NA deterministically
    d_chr <- as.character(d)
    d_chr[is.na(d_chr) | !nzchar(d_chr)] <- na_level

    # Define levels in a stable way: sorted unique + ensure na_level exists
    lev <- sort(unique(d_chr))
    if (!na_level %in% lev) lev <- c(lev, na_level)

    SummarizedExperiment::colData(spe)[[domain_col]] <- factor(d_chr, levels = lev)
    spe
}



.broadcast_param <- function(x, n, name) {
    if (length(x) == 1L) {
        rep(x, n)
    } else if (length(x) == n) {
        x
    } else {
        stop(
            sprintf(
                "Parameter '%s' must have length 1 or length(y) = %d; got %d",
                name, n, length(x)
            ),
            call. = FALSE
        )
    }
}

.safe_pull_re_mean <- function(re_list, name, idx_vec) {
    if (!is.list(re_list) || is.null(re_list[[name]])) {
        return(rep(NA_real_, length(idx_vec)))
    }

    re_df <- re_list[[name]]

    if (!("ID" %in% colnames(re_df)) || !("mean" %in% colnames(re_df))) {
        return(rep(NA_real_, length(idx_vec)))
    }

    ids <- re_df$ID
    m <- match(idx_vec, ids)
    out <- re_df$mean[m]
    out[is.na(out)] <- NA_real_
    out
}


#' Build block-diagonal kNN adjacency matrix (AOI-level)
#'
#' @param aoi AOI-level tibble from prep_inla_long_df()$aoi
#' @param k Number of nearest neighbors per AOI (within slide)
#' @param symmetric Whether to symmetrize adjacency
#'
#' @return Sparse symmetric adjacency matrix (dgCMatrix)
#' @export
build_block_diag_knn_adjacency <- function(
    aoi,
    k = 4L,
    symmetric = TRUE) {
    required <- c("AOI_ID", "AOI_x", "AOI_y", "Slide_ID")
    miss <- setdiff(required, names(aoi))
    if (length(miss) > 0L) {
        stop("Missing required AOI columns: ", paste(miss, collapse = ", "),
            call. = FALSE
        )
    }

    for (pkg in c("FNN", "Matrix", "dplyr")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    # Stable AOI ordering
    aoi <- dplyr::distinct(aoi, AOI_ID, AOI_x, AOI_y, Slide_ID) |>
        dplyr::arrange(AOI_ID)

    coords <- as.matrix(aoi[, c("AOI_x", "AOI_y")])
    slide <- aoi$Slide_ID
    slide_levels <- sort(unique(slide))

    adj_blocks <- vector("list", length(slide_levels))

    print("Building knn adjacency...")
    for (i in seq_along(slide_levels)) {
        s <- slide_levels[i]
        keep <- slide == s
        coords_s <- coords[keep, , drop = FALSE]
        n_s <- nrow(coords_s)

        if (n_s <= 1L) {
            adj_blocks[[i]] <- Matrix::Matrix(0, 1, 1, sparse = TRUE)
            next
        }

        k_eff <- min(k, n_s - 1L)
        nn <- FNN::get.knn(coords_s, k = k_eff)$nn.index

        adj <- matrix(0L, n_s, n_s)
        for (r in seq_len(n_s)) {
            adj[r, nn[r, ]] <- 1L
        }

        if (isTRUE(symmetric)) {
            adj <- ((adj + t(adj)) > 0L) * 1L
        }

        diag(adj) <- 0L
        adj_blocks[[i]] <- Matrix::Matrix(adj, sparse = TRUE)
    }

    Matrix::bdiag(adj_blocks)
}


#' Create AOI clusters for aggregation (slide-aware, knn merge)
#'
#' @param aoi AOI-level tibble (must include AOI_ID, Slide_ID, AOI_x_std, AOI_y_std)
#' @param target_per_slide Target number of clusters per slide
#' @param k Number of neighbors for merging
#' @param seed RNG seed
#' @param verbose Logical
#'
#' @return list(cluster_map, cluster_meta, meta)
#' @export
make_aoi_clusters <- function(
    aoi,
    target_per_slide = 600L,
    k = 5L,
    seed = 1L,
    verbose = FALSE) {
    for (pkg in c("dplyr", "FNN")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    required <- c("AOI_ID", "Slide_ID", "AOI_x_std", "AOI_y_std")
    miss <- setdiff(required, names(aoi))
    if (length(miss) > 0L) {
        stop("Missing required AOI columns: ", paste(miss, collapse = ", "),
            call. = FALSE
        )
    }

    set.seed(seed)

    slides <- split(aoi, aoi$Slide_ID)

    cluster_maps <- list()
    cluster_metas <- list()
    meta_slides <- list()

    for (sid in names(slides)) {
        aoi_s <- slides[[sid]]
        n_aoi <- nrow(aoi_s)

        if (verbose) {
            message("Slide ", sid, ": ", n_aoi, " AOIs")
        }

        if (n_aoi <= target_per_slide) {
            # no merging needed
            cluster_id <- seq_len(n_aoi)
            names(cluster_id) <- aoi_s$AOI_ID
        } else {
            # ---- Corrected knn-based clustering (size-limited) ----

            coords <- as.matrix(aoi_s[, c("AOI_x_std", "AOI_y_std")])

            knn <- FNN::get.knn(coords, k = min(k + 1L, n_aoi))
            nn <- knn$nn.index[, -1, drop = FALSE]

            target_size <- ceiling(n_aoi / target_per_slide)

            cluster_id <- integer(n_aoi)
            cluster_id[] <- NA_integer_

            current_cluster <- 0L
            order <- sample(seq_len(n_aoi)) # random traversal

            for (i in order) {
                if (!is.na(cluster_id[i])) next

                current_cluster <- current_cluster + 1L
                cluster_id[i] <- current_cluster
                members <- i

                for (j in nn[i, ]) {
                    if (length(members) >= target_size) break
                    if (is.na(cluster_id[j])) {
                        cluster_id[j] <- current_cluster
                        members <- c(members, j)
                    }
                }
            }

            # Any leftover unassigned AOIs (rare)
            na_idx <- which(is.na(cluster_id))
            if (length(na_idx) > 0L) {
                cluster_id[na_idx] <- current_cluster
            }

            names(cluster_id) <- aoi_s$AOI_ID
        }

        # ---- build cluster map ----
        map <- tibble::tibble(
            AOI_ID = aoi_s$AOI_ID,
            Slide_ID = as.integer(sid),
            CLUSTER_ID = as.integer(cluster_id)
        )

        # ---- cluster metadata ----
        meta <- aoi_s |>
            dplyr::mutate(CLUSTER_ID = as.integer(cluster_id)) |>
            dplyr::group_by(CLUSTER_ID) |>
            dplyr::summarise(
                Slide_ID = dplyr::first(Slide_ID),
                AOI_x_std = mean(AOI_x_std),
                AOI_y_std = mean(AOI_y_std),
                n_aoi = dplyr::n(),
                .groups = "drop"
            )

        cluster_maps[[sid]] <- map
        cluster_metas[[sid]] <- meta
        meta_slides[[sid]] <- list(
            n_aoi = n_aoi,
            n_cluster = nrow(meta)
        )

        if (verbose) {
            message(
                "  → ", nrow(meta), " clusters (target = ", target_per_slide, ")"
            )
        }
    }

    list(
        cluster_map = dplyr::bind_rows(cluster_maps),
        cluster_meta = dplyr::bind_rows(cluster_metas),
        meta = list(
            method = "knn_merge",
            target_per_slide = target_per_slide,
            k = k,
            slides = meta_slides
        )
    )
}


#' Aggregate AOI-level data to cluster-level for Stage B
#'
#' @param prep Output of prep_inla_long_df()
#' @param clusters Output of make_aoi_clusters()
#' @param tech_spatial Numeric vector of length n_AOI (Stage A spatial effect)
#' @param diff_mode How to aggregate diff covariate ("mean" or "majority")
#' @param verbose Logical
#'
#' @return list(df = aggregated_df, aoi = aggregated_aoi)
#' @export
aggregate_to_clusters <- function(
    prep,
    clusters,
    tech_spatial,
    group_col = "sim_group",
    aggregate_tech = TRUE,
    verbose = FALSE) {
    for (pkg in c("dplyr", "tidyr")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    df <- prep$df
    aoi <- prep$aoi

    if (length(tech_spatial) != nrow(aoi)) {
        stop("`tech_spatial` must be length nrow(prep$aoi).", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Attach cluster IDs to AOIs
    # ------------------------------------------------------------
    aoi2 <- aoi |>
        dplyr::left_join(
            clusters$cluster_map,
            by = c("AOI_ID", "Slide_ID")
        )

    if (anyNA(aoi2$CLUSTER_ID)) {
        stop("Some AOIs were not assigned to clusters.", call. = FALSE)
    }

    aoi2 <- aoi2 |>
        dplyr::mutate(
            tech_spatial = tech_spatial,
            lib_size = pmax(exp(log_lib), 1)
        )

    # ------------------------------------------------------------
    # Aggregate AOI metadata to CLUSTER level (per slide)
    # ------------------------------------------------------------
    aoi_agg <- aoi2 |>
        dplyr::group_by(Slide_ID, CLUSTER_ID) |>
        dplyr::summarise(
            AOI_x_std = mean(AOI_x_std),
            AOI_y_std = mean(AOI_y_std),
            AOI_x = mean(AOI_x_std),
            AOI_y = mean(AOI_y_std),
            tech_spatial = if (aggregate_tech) {
                weighted.mean(tech_spatial, w = lib_size, na.rm = TRUE)
            } else {
                0
            },
            lib_size = sum(lib_size),
            n_aoi = dplyr::n(),
            dplyr::across(
                dplyr::all_of(group_col),
                ~ as.numeric(mean(.x, na.rm = TRUE) > 0.5)
            ),
            .groups = "drop"
        ) |>
        dplyr::mutate(
            lib_size = pmax(lib_size, 1),
            log_lib  = log(lib_size)
        )

    # ------------------------------------------------------------
    # CRITICAL FIX: assign a *globally unique* AOI_ID
    # ------------------------------------------------------------
    aoi_agg <- aoi_agg |>
        dplyr::arrange(Slide_ID, CLUSTER_ID) |>
        dplyr::mutate(
            AOI_ID = dplyr::row_number()
        )

    # ------------------------------------------------------------
    # Attach cluster IDs to long df
    # ------------------------------------------------------------
    df2 <- df |>
        dplyr::left_join(
            clusters$cluster_map,
            by = c("AOI_ID", "Slide_ID")
        )

    if (anyNA(df2$CLUSTER_ID)) {
        stop("Some long-format rows were not assigned to clusters.", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Aggregate counts to cluster level
    # ------------------------------------------------------------
    df_agg <- df2 |>
        dplyr::group_by(
            Gene_ID, Gene_ID_de, Slide_ID, CLUSTER_ID
        ) |>
        dplyr::summarise(
            y = sum(y, na.rm = TRUE),
            diff = as.numeric(mean(diff, na.rm = TRUE) > 0.5),
            is_NEG = max(is_NEG),
            is_HK = max(is_HK),
            is_REG = max(is_REG),
            is_nonNEG = is_REG + is_HK,
            .groups = "drop"
        )

    # ------------------------------------------------------------
    # Join cluster-level offsets + global AOI_ID
    # ------------------------------------------------------------
    df_agg <- df_agg |>
        dplyr::left_join(
            aoi_agg |>
                dplyr::select(
                    Slide_ID, CLUSTER_ID, AOI_ID,
                    lib_size, log_lib, tech_spatial
                ),
            by = c("Slide_ID", "CLUSTER_ID")
        )

    # ------------------------------------------------------------
    # HARD GUARDRAILS (never remove these)
    # ------------------------------------------------------------
    stopifnot(
        !anyNA(df_agg$AOI_ID),
        !anyNA(df_agg$log_lib),
        all(is.finite(df_agg$log_lib)),
        !any(duplicated(df_agg[, c("Gene_ID", "AOI_ID")])),
        !any(duplicated(aoi_agg$AOI_ID))
    )

    if (verbose) {
        message(
            "Aggregated from ",
            nrow(df), " rows to ",
            nrow(df_agg), " rows (",
            nrow(aoi), " AOIs → ",
            nrow(aoi_agg), " clusters)"
        )
    }

    list(
        df  = df_agg,
        aoi = aoi_agg
    )
}


#'
#' @param fit_obj Output of normalize_inla_pac()
#' @param ref Reference anchor:
#'   "median","mean","median_bio","mean_bio",
#'   "median_bio_noHK","mean_bio_noHK"
#' @param method Projection method:
#'   "PAC","postmean","postmedian","eta"
#' @param housekeeping_col Optional column in prep_full$df for HK exclusion
#' @param assay_name Name of assay to add to spe
#'
#' @return SpatialExperiment with an added assay
#' @export
adjust_counts <- function(
    fit_obj,
    ref = c(
        "median", "mean",
        "median_bio", "mean_bio",
        "median_bio_noHK", "mean_bio_noHK"
    ),
    method = c("PAC", "postmean", "postmedian", "eta"),
    housekeeping_col = "is_HK",
    assay_out = "logcounts_readjusted") {
    ref <- match.arg(ref)
    method <- match.arg(method)

    spe <- fit_obj$spe
    fit <- fit_obj$fit

    if (is.null(fit_obj$prep_full) ||
        is.null(fit_obj$prep_full$df)) {
        stop("fit_obj$prep_full$df is required.", call. = FALSE)
    }

    df_long <- fit_obj$prep_full$df

    # ------------------------------------------------------------
    # Choose latent representation (cluster-safe)
    # ------------------------------------------------------------
    latent <- if (!is.null(fit$latent_full)) {
        fit$latent_full
    } else {
        fit$latent
    }

    need_latent <- c(
        "intercept", "gene_eff",
        "gene_de_eff", "log_lib"
    )

    miss <- setdiff(need_latent, names(latent))
    if (length(miss) > 0L) {
        stop(
            "Missing latent components: ",
            paste(miss, collapse = ", "),
            call. = FALSE
        )
    }

    # ------------------------------------------------------------
    # Define reference population
    # ------------------------------------------------------------
    keep <- rep(TRUE, length(latent$log_lib))

    if (grepl("_bio", ref)) {
        if (!"is_nonNEG" %in% names(df_long)) {
            stop(
                "prep_full$df must contain `is_nonNEG` ",
                "for *_bio references.",
                call. = FALSE
            )
        }
        keep <- keep & (df_long$is_nonNEG == 1)
    }

    if (grepl("noHK", ref)) {
        if (is.null(housekeeping_col) ||
            !housekeeping_col %in% names(df_long)) {
            stop(
                "Provide `housekeeping_col` present in ",
                "prep_full$df.",
                call. = FALSE
            )
        }
        keep <- keep &
            !(as.logical(df_long[[housekeeping_col]]))
    }

    if (!any(keep)) {
        stop(
            "No rows available for selected reference.",
            call. = FALSE
        )
    }

    ref_log_lib <- switch(ref,
        median = median(
            latent$log_lib,
            na.rm = TRUE
        ),
        mean = mean(
            latent$log_lib,
            na.rm = TRUE
        ),
        median_bio = median(
            latent$log_lib[keep],
            na.rm = TRUE
        ),
        mean_bio = mean(
            latent$log_lib[keep],
            na.rm = TRUE
        ),
        median_bio_noHK = median(
            latent$log_lib[keep],
            na.rm = TRUE
        ),
        mean_bio_noHK = mean(
            latent$log_lib[keep],
            na.rm = TRUE
        )
    )

    # ------------------------------------------------------------
    # Construct eta
    # ------------------------------------------------------------
    eta_full <-
        latent$intercept +
        latent$log_lib +
        latent$gene_eff +
        latent$gene_de_eff +
        latent$spatial_eff +
        latent$slide_eff

    eta_target <-
        latent$intercept +
        latent$gene_eff +
        latent$gene_de_eff +
        ref_log_lib

    # ------------------------------------------------------------
    # Project into expression space
    # ------------------------------------------------------------
    out_vec <- switch(method,
        PAC = {
            if (!"y" %in% names(df_long)) {
                stop(
                    "prep_full$df must contain y for PAC.",
                    call. = FALSE
                )
            }

            normalize_counts_pac_nb(
                y          = df_long$y,
                eta_full   = eta_full,
                eta_target = eta_target,
                size       = fit$size
            )
        },
        postmean = {
            log2(exp(eta_target) + 1)
        },
        postmedian = {
            if (is.null(fit$size)) {
                stop(
                    "fit must contain size for postmedian.",
                    call. = FALSE
                )
            }

            mu <- exp(eta_target)

            log2(
                stats::qnbinom(
                    0.5,
                    mu   = mu,
                    size = fit$size
                ) + 1
            )
        },
        eta = {
            eta_target
        }
    )

    # ------------------------------------------------------------
    # Map back to SPE geometry
    # ------------------------------------------------------------
    if (!all(c("Gene_ID", "AOI_ID") %in% names(df_long))) {
        stop(
            "prep_full$df must contain Gene_ID and AOI_ID.",
            call. = FALSE
        )
    }

    rd <- SummarizedExperiment::rowData(spe)
    cd <- SummarizedExperiment::colData(spe)

    gene_ids <- if ("Gene_ID" %in% colnames(rd)) {
        as.integer(rd$Gene_ID)
    } else {
        seq_len(nrow(spe))
    }

    aoi_ids <- if ("AOI_ID" %in% colnames(cd)) {
        as.integer(cd$AOI_ID)
    } else {
        seq_len(ncol(spe))
    }

    mat <- matrix(
        # NA_real_,
        0.0,
        nrow = nrow(spe),
        ncol = ncol(spe),
        dimnames = list(
            rownames(spe),
            colnames(spe)
        )
    )

    ok <- is.finite(out_vec)

    gi <- match(
        as.integer(df_long$Gene_ID[ok]),
        gene_ids
    )

    ai <- match(
        as.integer(df_long$AOI_ID[ok]),
        aoi_ids
    )

    if (anyNA(gi) || anyNA(ai)) {
        stop(
            "Failed to map Gene_ID/AOI_ID.",
            call. = FALSE
        )
    }

    mat[cbind(gi, ai)] <- out_vec[ok]

    assays <- SummarizedExperiment::assays(spe)
    assays[[assay_out]] <- mat
    SummarizedExperiment::assays(spe) <- assays

    spe
}
