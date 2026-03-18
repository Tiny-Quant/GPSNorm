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
    symmetric = TRUE
) {
    required <- c("AOI_ID", "AOI_x", "AOI_y", "Slide_ID")
    miss <- setdiff(required, names(aoi))
    if (length(miss) > 0L) {
        stop("Missing required AOI columns: ", paste(miss, collapse = ", "),
             call. = FALSE)
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
    slide  <- aoi$Slide_ID
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
    verbose = FALSE
) {

    for (pkg in c("dplyr", "FNN")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    required <- c("AOI_ID", "Slide_ID", "AOI_x_std", "AOI_y_std")
    miss <- setdiff(required, names(aoi))
    if (length(miss) > 0L) {
        stop("Missing required AOI columns: ", paste(miss, collapse = ", "),
             call. = FALSE)
    }

    set.seed(seed)

    slides <- split(aoi, aoi$Slide_ID)

    cluster_maps  <- list()
    cluster_metas <- list()
    meta_slides   <- list()

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
            nn  <- knn$nn.index[, -1, drop = FALSE]

            target_size <- ceiling(n_aoi / target_per_slide)

            cluster_id <- integer(n_aoi)
            cluster_id[] <- NA_integer_

            current_cluster <- 0L
            order <- sample(seq_len(n_aoi))  # random traversal

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

        cluster_maps[[sid]]  <- map
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
    verbose = FALSE
) {

    for (pkg in c("dplyr", "tidyr")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    df  <- prep$df
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
            is_HK  = max(is_HK),
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
    assay_out = "logcounts_readjusted"
) {
    ref    <- match.arg(ref)
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

    ref_log_lib <- switch(
        ref,
        median          = median(
            latent$log_lib,
            na.rm = TRUE
        ),
        mean            = mean(
            latent$log_lib,
            na.rm = TRUE
        ),
        median_bio      = median(
            latent$log_lib[keep],
            na.rm = TRUE
        ),
        mean_bio        = mean(
            latent$log_lib[keep],
            na.rm = TRUE
        ),
        median_bio_noHK = median(
            latent$log_lib[keep],
            na.rm = TRUE
        ),
        mean_bio_noHK   = mean(
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
    out_vec <- switch(
        method,

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
        #NA_real_,
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