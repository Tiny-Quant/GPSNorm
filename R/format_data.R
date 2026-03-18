# R/format_data.R
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
  extra_rowdata = c()
) {
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
    #rogene_type <- factor(gene_type, levels = c("NEG", "HK", "REG"))
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
    group_col = "Mix"
) {
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
    CodeClass_clean <- if("gene_type" %in% colnames(rd)) as.factor(rd$gene_type) else rep("REG", nrow(mat))
    AOI_ID  <- if ("AOI_ID" %in% colnames(cd)) as.integer(cd$AOI_ID) else seq_len(ncol(mat))

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
