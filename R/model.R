#' Randomized PIT for Negative Binomial observations
#'
#' Computes a randomized probability integral transform (PIT) for observations
#' from a Negative Binomial distribution. Randomization is required because the
#' NB distribution is discrete.
#'
#' @param y Integer vector of observed counts (>= 0).
#' @param mu Numeric vector of means (> 0).
#' @param size Numeric vector of NB size parameters (> 0).
#' @param randomized Logical; if TRUE (default), apply uniform randomization
#'   within the probability mass at y.
#' @param rng Optional numeric vector in (0,1) used for randomization.
#'   If NULL and randomized=TRUE, values are drawn via runif().
#'
#' @return Numeric vector u in (0, 1), same length as y.
#'
#' @details
#' The PIT is computed as:
#'   u = F(y - 1) + v * P(Y = y),
#' where v ~ Uniform(0,1).
#'
#' This ensures u is continuous even for discrete distributions.
#'
#' @export
compute_pit_nb <- function(
    y,
    mu,
    size,
    randomized = TRUE,
    rng = NULL
) {
    # Data formatting and input validation. 
    mu   <- .broadcast_param(mu, length(y), "mu")
    size <- .broadcast_param(size, length(y), "size")

    y <- as.integer(y)
    mu <- as.numeric(mu)
    size <- as.numeric(size)

    if (any(y < 0, na.rm = TRUE)) {
        stop("`y` must be non-negative.", call. = FALSE)
    }
    if (any(mu <= 0, na.rm = TRUE)) {
        stop("`mu` must be positive.", call. = FALSE)
    }
    if (any(size <= 0, na.rm = TRUE)) {
        stop("`size` must be positive.", call. = FALSE)
    }

    # Only compute quantiles where nothing is NA. 
    out <- rep(NA_real_, length(y))
    ok <- !(is.na(y) | is.na(mu) | is.na(size))
    if (!any(ok)) return(out) # Early exit. 

    y_ok <- y[ok]
    mu_ok <- mu[ok]
    size_ok <- size[ok]

    # F(y - 1)
    Fy_minus <- stats::pnbinom(
        q = pmax(y_ok - 1L, 0L),
        size = size_ok,
        mu = mu_ok
    )

    if (!randomized) {
        out[ok] <- Fy_minus
        return(out)
    }

    # P(Y = y)
    Py <- stats::dnbinom(
        x = y_ok,
        size = size_ok,
        mu = mu_ok
    )

    # Randomization
    if (is.null(rng)) {
        v <- stats::runif(length(y_ok))
    } else {
        if (length(rng) != length(y_ok)) {
            stop("`rng` must have same length as non-missing y.", call. = FALSE)
        }
        if (any(rng <= 0 | rng >= 1)) {
            stop("`rng` must lie strictly in (0,1).", call. = FALSE)
        }
        v <- rng
    }

    out[ok] <- Fy_minus + v * Py

    # Numerical guardrails
    out[out <= 0] <- .Machine$double.eps
    out[out >= 1] <- 1 - .Machine$double.eps

    out
}


#' Percentile-adjusted counts (PAC) for Negative Binomial models
#'
#' Maps PIT values into a target Negative Binomial distribution to obtain
#' percentile-invariant adjusted counts.
#'
#' @param u Numeric vector in (0,1), typically from `compute_pit_nb()`.
#' @param mu_target Numeric vector of target means (> 0).
#' @param size Numeric vector of NB size parameters (> 0).
#'
#' @return Integer vector of adjusted counts.
#'
#' @details
#' The adjusted count is computed as:
#'   y_tilde = F^{-1}(u | mu_target, size)
#'
#' This preserves ranks and local variability while changing the mean structure.
#'
#' @export
compute_pac_nb <- function(
    u,
    mu_target,
    size, 
    jitter = TRUE 
) {
    # Data formatting and input validation. 
    mu_target <- .broadcast_param(mu_target, length(u), "mu_target")
    size <- .broadcast_param(size, length(u), "size")

    u <- as.numeric(u)
    mu_target <- as.numeric(mu_target)
    size <- as.numeric(size)

    if (any(mu_target <= 0, na.rm = TRUE)) {
        stop("`mu_target` must be positive.", call. = FALSE)
    }
    if (any(size <= 0, na.rm = TRUE)) {
        stop("`size` must be positive.", call. = FALSE)
    }

    # Only compute where nothing is NA. 
    out <- rep(NA_real_, length(u))
    ok <- !(is.na(u) | is.na(mu_target) | is.na(size))
    if (!any(ok)) return(out) # Early exit. 

    u_ok <- u[ok]

    # Guard against boundary issues
    u_ok[u_ok <= 0] <- .Machine$double.eps
    u_ok[u_ok >= 1] <- 1 - .Machine$double.eps

    out[ok] <- stats::qnbinom(
        p = u_ok,
        size = size[ok],
        mu = mu_target[ok]
    )

    if (jitter) {
        out[ok] <- out[ok] + u_ok
    }

    out
}


#' PAC-based normalization under Negative Binomial
#'
#' Applies percentile-invariant transformation using a full NB model
#' and a target NB model (with selected effects removed).
#'
#' @param y Integer vector of observed counts
#' @param eta_full Linear predictor including all effects
#' @param eta_target Linear predictor with only retained effects
#' @param size NB size (1 / overdispersion), scalar or vector
#'
#' @return Numeric vector: log2(PAC + 1)
#'
#' @details
#' This implements:
#'   1. PIT under full model
#'   2. Inverse CDF under target model
#'   3. log2(PAC + 1)
#'
#' This mirrors SpaNorm's PAC, but with INLA-based parameters.
#'
#' @export
normalize_counts_pac_nb <- function(
    y,
    eta_full,
    eta_target,
    size, 
    jitter = TRUE
) {

    # Data formatting and input validation. 
    eta_full   <- .broadcast_param(eta_full, length(y), "eta_full")
    eta_target <- .broadcast_param(eta_target, length(y), "eta_target")
    size       <- .broadcast_param(size, length(y), "size")

    mu_full   <- exp(eta_full)
    mu_target <- exp(eta_target)

    ## Step 1: PIT under full model
    u <- compute_pit_nb(
        y    = y,
        mu   = mu_full,
        size = size
    )

    ## Step 2: PAC under target model
    y_pac <- compute_pac_nb(
        u = u,
        mu_target = mu_target,
        size = size,
        jitter = jitter
    )

    ## Step 3: log2(PAC + 1)
    log2(y_pac + 1)
}


#' Prepare long-form GeoMx data for INLA normalization models
#'
#' This function performs *only* data preparation and returns:
#'   (1) a long-form gene × AOI dataframe for likelihood modeling
#'   (2) an AOI-level metadata table for spatial models (BESAG / SPDE / SVG)
#'
#' No model fitting or spatial objects are created here.
#'
#' @param data Long-form data.frame (gene × AOI)
#' @param response Name of count column
#' @param diff_group Name of DE group column
#' @param tech_covariates Optional technical covariates
#' @param center_scale_covariates Logical
#' @param subsample_frac Optional AOI-level subsampling fraction
#'
#' @return A list with elements:
#'   \describe{
#'     \item{df}{Long-form dataframe for INLA likelihood}
#'     \item{aoi}{AOI-level spatial metadata (one row per AOI)}
#'   }
#' @export
prep_inla_long_df <- function(
    data,
    response = "Count_int",
    diff_group,
    tech_covariates = NULL,
    center_scale_covariates = FALSE,
    subsample_frac = NULL
) {

    for (pkg in c("dplyr")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    stopifnot(is.data.frame(data))
    stopifnot(is.character(response), length(response) == 1L)
    stopifnot(is.character(diff_group), length(diff_group) == 1L)

    required <- c(
        response,
        diff_group,
        "AOI_ID",
        "Slide_ID",
        "Gene_ID",
        "CodeClass_clean",
        "ROICoordinateX",
        "ROICoordinateY"
    )

    miss <- setdiff(required, names(data))
    if (length(miss) > 0L) {
        stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }

    if (is.null(tech_covariates)) tech_covariates <- character(0)
    tech_covariates <- unique(tech_covariates)

    miss_tc <- setdiff(tech_covariates, names(data))
    if (length(miss_tc) > 0L) {
        stop("Missing technical covariates: ", paste(miss_tc, collapse = ", "), call. = FALSE)
    }

    df <- dplyr::as_tibble(data)

    # ---------------------------------------------------------
    # AOI-level metadata (FIRST-CLASS OBJECT)
    # ---------------------------------------------------------
    aoi <- df |>
        dplyr::distinct(
            AOI_ID,
            Slide_ID,
            AOI_x = ROICoordinateX,
            AOI_y = ROICoordinateY
        ) |>
        dplyr::arrange(AOI_ID)

    if (anyDuplicated(aoi$AOI_ID)) {
        stop("AOI_ID must be unique at the AOI level.", call. = FALSE)
    }

    # ---------------------------------------------------------
    # AOI-level DE grouping (FIRST-CLASS)
    # ---------------------------------------------------------
    aoi_group <- df |>
        dplyr::distinct(
            AOI_ID,
            Slide_ID,
            .data[[diff_group]]
        )

    # sanity: group must be constant within AOI
    chk <- aoi_group |>
        dplyr::count(AOI_ID) |>
        dplyr::filter(n > 1)

    if (nrow(chk) > 0L) {
        stop(
            "DE group varies within AOI. ",
            "diff_group must be AOI-level.",
            call. = FALSE
        )
    }

    aoi <- aoi |>
        dplyr::left_join(
            aoi_group,
            by = c("AOI_ID", "Slide_ID")
        )

    # ---------------------------------------------------------
    # Core row-level variables
    # ---------------------------------------------------------
    df <- df |>
        dplyr::mutate(
            y = as.integer(.data[[response]]),

            diff_raw = as.numeric(.data[[diff_group]]),
            diff     = diff_raw * (.data$CodeClass_clean == "REG"),

            is_NEG    = as.integer(.data$CodeClass_clean == "NEG"),
            is_HK     = as.integer(.data$CodeClass_clean == "HK"),
            is_REG    = as.integer(.data$CodeClass_clean == "REG"),
            is_nonNEG = as.integer(.data$CodeClass_clean != "NEG"),

            Gene_ID    = as.integer(Gene_ID),
            Gene_ID_de = as.integer(Gene_ID),

            AOI_ID   = as.integer(AOI_ID),
            Slide_ID = as.integer(Slide_ID)
        )

    # ---------------------------------------------------------
    # Library size (AOI-level, joined back)
    # ---------------------------------------------------------
    lib_tbl <- df |>
        dplyr::group_by(AOI_ID) |>
        dplyr::summarise(
            lib_size = sum(y, na.rm = TRUE),
            .groups = "drop"
        )

    df <- df |>
        dplyr::left_join(lib_tbl, by = "AOI_ID") |>
        dplyr::mutate(
            lib_size = pmax(lib_size, 1),
            log_lib  = log(lib_size)
        )

    aoi <- aoi |>
        dplyr::left_join(
            df |>
            dplyr::distinct(AOI_ID, log_lib),
            by = "AOI_ID"
        )

    # ---------------------------------------------------------
    # AOI-level NEG / HK summaries (technical-only signals)
    # ---------------------------------------------------------
    tech_sum_tbl <- df |>
        dplyr::group_by(AOI_ID) |>
        dplyr::summarise(
            neg_sum = sum(y[is_NEG == 1L], na.rm = TRUE),
            hk_sum  = sum(y[is_HK  == 1L], na.rm = TRUE),
            .groups = "drop"
        )

    aoi <- aoi |>
        dplyr::left_join(tech_sum_tbl, by = "AOI_ID") |>
        dplyr::mutate(
            neg_sum = dplyr::coalesce(neg_sum, 0L),
            hk_sum  = dplyr::coalesce(hk_sum, 0L)
        )

    # ---------------------------------------------------------
    # Standardize to median nearest-neighbor distance per slide
    # ---------------------------------------------------------
    aoi <- aoi |>
        dplyr::group_by(Slide_ID) |>
        dplyr::group_modify(function(.x, .g) {
            coords <- as.matrix(.x[, c("AOI_x", "AOI_y")])

            # k=2: first neighbor is itself, second is true NN
            nn <- FNN::get.knn(coords, k = 2)$nn.dist[, 2]
            s  <- stats::median(nn, na.rm = TRUE)

            if (!is.finite(s) || s <= 0) {
                # fallback if degenerate
                s <- stats::median(abs(coords[, 1] - stats::median(coords[, 1]))) +
                        stats::median(abs(coords[, 2] - stats::median(coords[, 2])))
                s <- ifelse(!is.finite(s) || s <= 0, 1, s)
            }

            .x$AOI_x_std <- .x$AOI_x / s
            .x$AOI_y_std <- .x$AOI_y / s

            # store per-slide NN summaries (for tuning diagnostics)
            .x$nn_med <- s
            .x$nn_p25 <- stats::quantile(nn, 0.25, na.rm = TRUE, names = FALSE)
            .x$nn_p75 <- stats::quantile(nn, 0.75, na.rm = TRUE, names = FALSE)

            .x
        }) |>
        dplyr::ungroup()

    # ---------------------------------------------------------
    # Technical covariates
    # ---------------------------------------------------------
    if (length(tech_covariates) > 0L) {
        for (v in tech_covariates) {
            df[[v]] <- as.numeric(df[[v]])
            if (center_scale_covariates) {
                df[[v]] <- as.numeric(scale(df[[v]]))
            }
        }
    }

    # ---------------------------------------------------------
    # AOI-level likelihood subsampling
    # ---------------------------------------------------------
    if (!is.null(subsample_frac)) {

        if (!is.numeric(subsample_frac) ||
            length(subsample_frac) != 1L ||
            subsample_frac <= 0 || subsample_frac > 1) {
            stop("`subsample_frac` must be in (0,1].", call. = FALSE)
        }

        keep_aoi <- aoi |>
            dplyr::group_by(Slide_ID) |>
            dplyr::slice_sample(prop = subsample_frac) |>
            dplyr::ungroup() |>
            dplyr::pull(AOI_ID)

        df$y[!df$AOI_ID %in% keep_aoi] <- NA_integer_
    }

    # ---------------------------------------------------------
    # Final row-level selection (NO COORDINATES)
    # ---------------------------------------------------------
    keep <- c(
        "y",
        "AOI_ID",
        "Slide_ID",
        "Gene_ID",
        "Gene_ID_de",
        "diff",
        "diff_raw",
        "is_NEG",
        "is_HK",
        "is_REG",
        "is_nonNEG",
        "lib_size",
        "log_lib",
        tech_covariates
    )

    df <- df |>
        dplyr::select(dplyr::all_of(keep))

    prep <- list(
        df  = df,
        aoi = aoi
    )

    prep
}


#' Normalize a SpatialExperiment using an INLA + PAC model
#'
#' @return List with SpatialExperiment (`spe`) and fitted model (`fit`)
#' @export
normalize_spe_inla_pac <- function(
    spe,
    model_fun,
    group_col = "sim_group",
    counts_assay = "counts",
    tech_covariates = NULL,
    center_scale_covariates = FALSE,
    subsample_frac = NULL,
    model_args = list(),
    cluster_aoi = FALSE,
    cluster_args = list()
) {

    validate_spe_contract(spe, assay = counts_assay)

    # ------------------------------------------------------------
    # Long format (FULL RESOLUTION)  ---- keep this for projection!
    # ------------------------------------------------------------
    long_full <- spe_to_geomx_long(
        spe = spe,
        assay = counts_assay,
        group_col = group_col
    )

    prep_full <- prep_inla_long_df(
        data = long_full,
        diff_group = group_col,
        tech_covariates = tech_covariates,
        center_scale_covariates = center_scale_covariates,
        subsample_frac = subsample_frac
    )

    df_full  <- prep_full$df
    aoi_full <- prep_full$aoi

    # ------------------------------------------------------------
    # Optional AOI clustering for FITTING ONLY
    #   - We fit on clustered data
    #   - But we return normalized counts at FULL resolution
    # ------------------------------------------------------------
    clusters <- NULL
    prep_fit <- prep_full

    if (isTRUE(cluster_aoi)) {

        clusters <- do.call(
            make_aoi_clusters,
            c(list(aoi = prep_full$aoi), cluster_args)
        )

        agg <- aggregate_to_clusters(
            prep = prep_full,
            clusters = clusters,
            tech_spatial = rep(0, nrow(prep_full$aoi)),  # placeholder
            group_col = group_col
        )

        prep_fit$df  <- agg$df
        prep_fit$aoi <- agg$aoi

        # Guardrails: clustered df must be unique per (Gene_ID, AOI_ID)
        if (any(duplicated(prep_fit$df[, c("Gene_ID", "AOI_ID")]))) {
            stop(
                "Duplicate (Gene_ID, AOI_ID) pairs after clustering. ",
                "aggregate_to_clusters() must produce unique rows.",
                call. = FALSE
            )
        }

        if (any(!is.finite(prep_fit$aoi$log_lib))) {
            stop(
                "Non-finite log_lib after clustering. ",
                "Check aggregate_to_clusters().",
                call. = FALSE
            )
        }
    }

    # ------------------------------------------------------------
    # Fit INLA model (possibly on clusters)
    # ------------------------------------------------------------
    model_out <- do.call(
        model_fun,
        c(list(prep = prep_fit), model_args)
    )

    required <- c("eta_full", "size")
    missing <- setdiff(required, names(model_out))
    if (length(missing) > 0L) {
        stop(
            "Model function did not return: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }

    # -----------------------------
    # Construct desired target mean
    # -----------------------------
    eta_target_default <- (
        model_out$eta_full
        - model_out$latent$spatial_eff
        - model_out$latent$slide_eff
        - model_out$latent$log_lib 
        + median(model_out$latent$log_lib, na.rm = TRUE)
    )

    # We REQUIRE the INLA fit object for projection
    if (isTRUE(cluster_aoi)) {
        if (is.null(model_out$fit)) {
            stop(
                "When cluster_aoi=TRUE, model_fun must return `fit` (INLA object) ",
                "so we can project effects back to full AOI resolution.",
                call. = FALSE
            )
        }
    }

    # ------------------------------------------------------------
    # Projection: clustered fit -> full-resolution latent + eta
    # ------------------------------------------------------------
    eta_full_use   <- NULL
    eta_target_use <- NULL
    y_use          <- NULL
    latent_full    <- NULL

    if (isTRUE(cluster_aoi)) {

        fit <- model_out$fit

        # --- helper to pull iid/besag means by ID ---
        pull_re <- function(name, ids) {
            .safe_pull_re_mean(fit$summary.random, name, ids)
        }

        # --- fixed intercept ---
        if (is.null(fit$summary.fixed) || !"(Intercept)" %in% rownames(fit$summary.fixed)) {
            stop("INLA fit missing fixed intercept.", call. = FALSE)
        }
        intercept <- fit$summary.fixed["(Intercept)", "mean"]

        # --- map each AOI to its cluster (slide-aware) ---
        cmap <- clusters$cluster_map
        cmap$AOI_ID     <- as.integer(cmap$AOI_ID)
        cmap$Slide_ID   <- as.integer(cmap$Slide_ID)
        cmap$CLUSTER_ID <- as.integer(cmap$CLUSTER_ID)

        key_full <- paste(df_full$AOI_ID, df_full$Slide_ID, sep = "_")
        key_map  <- paste(cmap$AOI_ID,   cmap$Slide_ID,   sep = "_")
        cluster_id_full <- cmap$CLUSTER_ID[match(key_full, key_map)]
        if (anyNA(cluster_id_full)) {
            stop("Failed to map some full-resolution AOIs to clusters.", call. = FALSE)
        }

        spatial_eff <- pull_re("AOI_ID", cluster_id_full)
        slide_eff   <- pull_re("Slide_ID", df_full$Slide_ID)

        gene_eff_raw <- pull_re("Gene_ID", df_full$Gene_ID)
        gene_eff     <- gene_eff_raw * df_full$is_nonNEG

        gene_de_raw  <- pull_re("Gene_ID_de", df_full$Gene_ID_de)
        gene_de_eff  <- gene_de_raw * df_full$diff

        log_lib <- df_full$log_lib

        latent_full <- list(
            intercept    = intercept,
            spatial_eff  = spatial_eff,
            slide_eff    = slide_eff,
            gene_eff     = gene_eff,
            gene_de_eff  = gene_de_eff,
            log_lib      = log_lib
        )

        model_out$latent_full <- latent_full

        eta_full_use <-
            latent_full$intercept +
            latent_full$log_lib +
            latent_full$gene_eff +
            latent_full$gene_de_eff +
            latent_full$spatial_eff +
            latent_full$slide_eff

        ref_log_lib <- median(df_full$log_lib, na.rm = TRUE)

        eta_target_use <-
            eta_full_use -
            latent_full$spatial_eff -
            latent_full$slide_eff -
            latent_full$log_lib +
            ref_log_lib

        y_use <- df_full$y

    } else {

        eta_full_use   <- model_out$eta_full
        eta_target_use <- eta_target_default
        y_use          <- prep_fit$df$y
    }

    # ------------------------------------------------------------
    # PAC normalization (at the resolution of y_use)
    # ------------------------------------------------------------
    logcounts_vec <- normalize_counts_pac_nb(
        y          = y_use,
        eta_full   = eta_full_use,
        eta_target = eta_target_use,
        size       = model_out$size
    )

    # ------------------------------------------------------------
    # Map back to FULL-resolution SPE
    #   - ALWAYS return original SPE with original colData (includes group_col)
    # ------------------------------------------------------------
    rd <- SummarizedExperiment::rowData(spe)
    cd <- SummarizedExperiment::colData(spe)

    gene_ids <- if ("Gene_ID" %in% colnames(rd)) as.integer(rd$Gene_ID) else seq_len(nrow(spe))
    aoi_ids  <- if ("AOI_ID"  %in% colnames(cd)) as.integer(cd$AOI_ID)  else seq_len(ncol(spe))

    mat_norm <- matrix(
        #NA_real_,
        0.0, 
        nrow = nrow(spe),
        ncol = ncol(spe),
        dimnames = list(rownames(spe), colnames(spe))
    )

    # Use the FULL long mapping for write-back
    ok <- is.finite(logcounts_vec) & !is.na(long_full$Gene_ID) & !is.na(long_full$AOI_ID)
    if (sum(ok) == 0L) {
        stop("All PAC-normalized values are non-finite.", call. = FALSE)
    }

    gi <- match(as.integer(long_full$Gene_ID[ok]), gene_ids)
    ai <- match(as.integer(long_full$AOI_ID[ok]),  aoi_ids)

    if (anyNA(gi) || anyNA(ai)) {
        stop("Failed to map Gene_ID or AOI_ID back to original SPE.", call. = FALSE)
    }

    mat_norm[cbind(gi, ai)] <- logcounts_vec[ok]

    # Final guardrail: no AOI should be entirely missing
    bad_cols <- which(colSums(is.finite(mat_norm)) == 0L)
    if (length(bad_cols) > 0L) {
        stop(
            "Some AOIs have no finite normalized values (",
            length(bad_cols), " columns). This will break downstream DE.",
            call. = FALSE
        )
    }

    assays <- SummarizedExperiment::assays(spe)
    assays$logcounts <- mat_norm
    SummarizedExperiment::assays(spe) <- assays

    # ------------------------------------------------------------
    # Return: full-resolution SPE + fit + optional clustering metadata
    # ------------------------------------------------------------
    out <- list(
        spe = spe,
        fit = model_out, 
        prep_full = prep_full
    )

    if (isTRUE(cluster_aoi)) {
        out$clusters <- clusters
        out$prep_fit <- prep_fit   # clustered prep used for fitting
    }

    out
}


#' Baseline INLA normalization model: gene baseline + spatial AOI
#'
#' @export
fit_inla_baseline_model <- function(
    prep,
    k = 4, 
    symmetric = TRUE, 
    family = "nbinomial",
    verbose = FALSE
) {
    df <- prep$df
    required <- c(
        "y", "log_lib", "Gene_ID", "Gene_ID_de", "diff", "is_nonNEG", 
        "AOI_ID", "Slide_ID" 
    )
    miss <- setdiff(required, names(df))
    if (length(miss) > 0L) {
        stop(
            "Missing required columns: ", 
            paste(miss, collapse = ", "), call. = FALSE
        )
    }

    adj <- build_block_diag_knn_adjacency(
        prep$aoi, k = k, symmetric = symmetric
    )
    graph <- INLA::inla.read.graph(adj) 
    if (is.null(graph)) {
        stop("Missing INLA graph.")
    }

    formula <- (
        y ~ 1 
            + offset(log_lib) 
            + f(Gene_ID, is_nonNEG, model = "iid") 
            + f(Gene_ID_de, diff, model = "iid")
            + f(AOI_ID, model = "besag", graph = graph, scale.model = TRUE)
            + f(Slide_ID, model = "iid")
    )

    fit <- INLA::inla(
        formula,
        data = df,
        family = family,
        control.predictor = list(compute = TRUE, link = 1),
        verbose = verbose
    )

    eta_full <- fit$summary.linear.predictor$mean

    # ---- Extract components ----
    intercept <- fit$summary.fixed["(Intercept)", "mean"]

    spatial_eff <- .safe_pull_re_mean(
        fit$summary.random,
        "AOI_ID",
        df$AOI_ID
    )

    slide_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Slide_ID",
        df$Slide_ID
    )

    gene_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Gene_ID",
        df$Gene_ID
    ) * df$is_nonNEG

    gene_de_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Gene_ID_de",
        df$Gene_ID_de
    ) * df$diff

    latent_effs <- list(
        intercept = intercept, 
        spatial_eff = spatial_eff, 
        slide_eff = slide_eff, 
        gene_eff = gene_eff,
        gene_de_eff = gene_de_eff,
        log_lib = df$log_lib
    )

    size <- fit$summary.hyperpar[
        "size for the nbinomial observations (1/overdispersion)",
        "mean"
    ]

    list(
        fit = fit,
        eta_full = eta_full,
        size = size, 
        latent = latent_effs
    )
}


#'
#' @export
fit_inla_tuned_pc <- function(
    prep,
    k = 4, 
    symmetric = TRUE, 
    family = "nbinomial",
    de_U = 0.3,  
    de_alpha = 0.1,  
    verbose = FALSE
) {
    df <- prep$df
    required <- c(
        "y", "log_lib", "Gene_ID", "Gene_ID_de", "diff", "is_nonNEG", 
        "AOI_ID", "Slide_ID" 
    )
    miss <- setdiff(required, names(df))
    if (length(miss) > 0L) {
        stop(
            "Missing required columns: ", 
            paste(miss, collapse = ", "), call. = FALSE
        )
    }

    adj <- build_block_diag_knn_adjacency(
        prep$aoi, k = k, symmetric = symmetric
    )
    graph <- INLA::inla.read.graph(adj) 
    if (is.null(graph)) {
        stop("Missing INLA graph.")
    }

    formula <- (
        y ~ 1 
            + offset(log_lib) 
            + f(Gene_ID, is_nonNEG, model = "iid") 
            + f(
                Gene_ID_de, diff, model = "iid", 
                hyper = list(
                    prec = list(
                        prior = "pc.prec", 
                        param = c(de_U, de_alpha)
                    )
                )
            )
            + f(AOI_ID, model = "besag", graph = graph, scale.model = TRUE)
            + f(Slide_ID, model = "iid")
    )

    fit <- INLA::inla(
        formula,
        data = df,
        family = family,
        control.predictor = list(compute = TRUE, link = 1),
        verbose = verbose
    )

    eta_full <- fit$summary.linear.predictor$mean

    # ---- Extract components ----
    intercept <- fit$summary.fixed["(Intercept)", "mean"]

    spatial_eff <- .safe_pull_re_mean(
        fit$summary.random,
        "AOI_ID",
        df$AOI_ID
    )

    slide_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Slide_ID",
        df$Slide_ID
    )

    gene_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Gene_ID",
        df$Gene_ID
    ) * df$is_nonNEG

    gene_de_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Gene_ID_de",
        df$Gene_ID_de
    ) * df$diff

    latent_effs <- list(
        intercept = intercept, 
        spatial_eff = spatial_eff, 
        slide_eff = slide_eff, 
        gene_eff = gene_eff,
        gene_de_eff = gene_de_eff,
        log_lib = df$log_lib
    )

    size <- fit$summary.hyperpar[
        "size for the nbinomial observations (1/overdispersion)",
        "mean"
    ]

    list(
        fit = fit,
        eta_full = eta_full,
        size = size, 
        latent = latent_effs
    )
}


#'
#' @export 
fit_inla_heavy_tail <- function(
    prep,
    k = 4, 
    symmetric = TRUE, 
    family = "nbinomial",
    de_lgamma = c(a = 0.5, b = 0.15),
    verbose = FALSE
) {
    df <- prep$df
    required <- c(
        "y", "log_lib", "Gene_ID", "Gene_ID_de", "diff", "is_nonNEG", 
        "AOI_ID", "Slide_ID" 
    )
    miss <- setdiff(required, names(df))
    if (length(miss) > 0L) {
        stop(
            "Missing required columns: ", 
            paste(miss, collapse = ", "), call. = FALSE
        )
    }

    adj <- build_block_diag_knn_adjacency(
        prep$aoi, k = k, symmetric = symmetric
    )
    graph <- INLA::inla.read.graph(adj) 
    if (is.null(graph)) {
        stop("Missing INLA graph.")
    }

    formula <- (
        y ~ 1 
            + offset(log_lib) 
            + f(Gene_ID, is_nonNEG, model = "iid") 
            + f(
                Gene_ID_de, diff, model = "iid",
                hyper = list(
                    prec = list(
                        prior = "loggamma", param = de_lgamma
                    )
                )
            )
            + f(AOI_ID, model = "besag", graph = graph, scale.model = TRUE)
            + f(Slide_ID, model = "iid")
    )

    fit <- INLA::inla(
        formula,
        data = df,
        family = family,
        control.predictor = list(compute = TRUE, link = 1),
        verbose = verbose
    )

    eta_full <- fit$summary.linear.predictor$mean

    # ---- Extract components ----
    intercept <- fit$summary.fixed["(Intercept)", "mean"]

    spatial_eff <- .safe_pull_re_mean(
        fit$summary.random,
        "AOI_ID",
        df$AOI_ID
    )

    slide_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Slide_ID",
        df$Slide_ID
    )

    gene_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Gene_ID",
        df$Gene_ID
    ) * df$is_nonNEG

    gene_de_eff <- .safe_pull_re_mean(
        fit$summary.random, 
        "Gene_ID_de",
        df$Gene_ID_de
    ) * df$diff

    latent_effs <- list(
        intercept = intercept, 
        spatial_eff = spatial_eff, 
        slide_eff = slide_eff, 
        gene_eff = gene_eff,
        gene_de_eff = gene_de_eff,
        log_lib = df$log_lib
    )

    size <- fit$summary.hyperpar[
        "size for the nbinomial observations (1/overdispersion)",
        "mean"
    ]

    list(
        fit = fit,
        eta_full = eta_full,
        size = size, 
        latent = latent_effs
    )
}


#'
#' @export
fit_inla_bin <- function(
    prep,
    k = 4,
    symmetric = TRUE,
    family = "nbinomial",
    n_de_bins = 5,
    de_U = 0.30,
    de_alpha = 0.10,
    verbose = FALSE
) {

    df <- prep$df

    required <- c(
        "y", "log_lib", "Gene_ID",
        "Gene_ID_de", "diff", "is_nonNEG",
        "AOI_ID", "Slide_ID"
    )

    miss <- setdiff(required, names(df))
    if (length(miss) > 0L) {
        stop("Missing required columns: ",
             paste(miss, collapse = ", "), call. = FALSE)
    }

    # ------------------------------------------------------------
    # 1️⃣ Compute gene-level mean expression internally
    # ------------------------------------------------------------

    gene_mean <- df |>
        dplyr::group_by(Gene_ID_de) |>
        dplyr::summarise(
            mean_expr = mean(y),
            .groups = "drop"
        )

    # ------------------------------------------------------------
    # 2️⃣ Create equal-frequency bins
    # ------------------------------------------------------------

    gene_mean$de_bin <- dplyr::ntile(
        gene_mean$mean_expr,
        n_de_bins
    )

    # map bins back to long df
    df$de_bin <- gene_mean$de_bin[
        match(df$Gene_ID_de, gene_mean$Gene_ID_de)
    ]

    # ------------------------------------------------------------
    # 3️⃣ Create bin-specific DE index
    # ------------------------------------------------------------

    df$Gene_ID_de_bin <- interaction(
        df$Gene_ID_de,
        df$de_bin,
        drop = TRUE
    )

    # ------------------------------------------------------------
    # 4️⃣ Build adjacency graph
    # ------------------------------------------------------------

    adj <- build_block_diag_knn_adjacency(
        prep$aoi,
        k = k,
        symmetric = symmetric
    )

    graph <- INLA::inla.read.graph(adj)
    if (is.null(graph)) stop("Missing INLA graph.")

    # ------------------------------------------------------------
    # 5️⃣ Model formula
    # ------------------------------------------------------------

    formula <- (
        y ~ 1 +
            offset(log_lib) +

            # baseline gene effect
            f(Gene_ID, is_nonNEG, model = "iid") +

            # expression-stratified DE prior
            f(
                Gene_ID_de_bin,
                diff,
                model = "iid",
                hyper = list(
                    prec = list(
                        prior = "pc.prec",
                        param = c(de_U, de_alpha)
                    )
                )
            ) +

            f(
                AOI_ID,
                model = "besag",
                graph = graph,
                scale.model = TRUE
            ) +

            f(Slide_ID, model = "iid")
    )

    fit <- INLA::inla(
        formula,
        data = df,
        family = family,
        control.predictor = list(compute = TRUE, link = 1),
        verbose = verbose
    )

    # ------------------------------------------------------------
    # 6️⃣ Extract components (same structure as baseline)
    # ------------------------------------------------------------

    eta_full <- fit$summary.linear.predictor$mean
    intercept <- fit$summary.fixed["(Intercept)", "mean"]

    spatial_eff <- .safe_pull_re_mean(
        fit$summary.random, "AOI_ID", df$AOI_ID
    )

    slide_eff <- .safe_pull_re_mean(
        fit$summary.random, "Slide_ID", df$Slide_ID
    )

    gene_eff <- .safe_pull_re_mean(
        fit$summary.random, "Gene_ID", df$Gene_ID
    ) * df$is_nonNEG

    de_eff_bin <- .safe_pull_re_mean(
        fit$summary.random,
        "Gene_ID_de_bin",
        df$Gene_ID_de_bin
    )

    gene_de_eff <- de_eff_bin * df$diff

    size <- fit$summary.hyperpar[
        "size for the nbinomial observations (1/overdispersion)",
        "mean"
    ]

    list(
        fit = fit,
        eta_full = eta_full,
        size = size,
        latent = list(
            intercept = intercept,
            spatial_eff = spatial_eff,
            slide_eff = slide_eff,
            gene_eff = gene_eff,
            gene_de_eff = gene_de_eff,
            log_lib = df$log_lib
        )
    )
}


#' Baseline INLA model with two-component DE mixture prior
#'
#' Heavy-tailed approximation:
#'   beta_de = beta_small + beta_large
#'
#' @export
fit_inla_mix <- function(
    prep,
    k = 4,
    symmetric = TRUE,
    family = "nbinomial",
    verbose = FALSE,

    # ---- mixture prior controls ----
    small_U = 0.10,
    small_alpha = 0.50,
    large_U = 0.75,
    large_alpha = 0.01
) {

    df <- prep$df

    required <- c(
        "y", "log_lib", "Gene_ID", "Gene_ID_de",
        "diff", "is_nonNEG", "AOI_ID", "Slide_ID"
    )

    miss <- setdiff(required, names(df))
    if (length(miss) > 0L) {
        stop(
            "Missing required columns: ",
            paste(miss, collapse = ", "),
            call. = FALSE
        )
    }

    # ------------------------------------------------------------
    # NEW: duplicate DE indices for mixture components
    # ------------------------------------------------------------
    df$Gene_ID_de_small <- df$Gene_ID_de
    df$Gene_ID_de_large <- df$Gene_ID_de

    # ------------------------------------------------------------
    # Graph construction
    # ------------------------------------------------------------
    adj <- build_block_diag_knn_adjacency(
        prep$aoi,
        k = k,
        symmetric = symmetric
    )

    graph <- INLA::inla.read.graph(adj)

    if (is.null(graph)) {
        stop("Missing INLA graph.", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Mixture DE prior formula
    # ------------------------------------------------------------
    formula <- (
        y ~ 1 +
            offset(log_lib) +

            # gene baseline
            f(Gene_ID, is_nonNEG, model = "iid") +

            # ---- NEW mixture DE prior ----
            f(
                Gene_ID_de_small,
                diff,
                model = "iid",
                hyper = list(
                    prec = list(
                        prior = "pc.prec",
                        param = c(small_U, small_alpha)
                    )
                )
            ) +
            f(
                Gene_ID_de_large,
                diff,
                model = "iid",
                hyper = list(
                    prec = list(
                        prior = "pc.prec",
                        param = c(large_U, large_alpha)
                    )
                )
            ) +

            # shared spatial technical field
            f(
                AOI_ID,
                model = "besag",
                graph = graph,
                scale.model = TRUE
            ) +

            # slide effect
            f(Slide_ID, model = "iid")
    )

    # ------------------------------------------------------------
    # Fit model
    # ------------------------------------------------------------
    fit <- INLA::inla(
        formula,
        data = df,
        family = family,
        control.predictor = list(
            compute = TRUE,
            link = 1
        ),
        verbose = verbose
    )

    eta_full <- fit$summary.linear.predictor$mean

    # ------------------------------------------------------------
    # Extract latent effects
    # ------------------------------------------------------------
    intercept <- fit$summary.fixed["(Intercept)", "mean"]

    spatial_eff <- .safe_pull_re_mean(
        fit$summary.random,
        "AOI_ID",
        df$AOI_ID
    )

    slide_eff <- .safe_pull_re_mean(
        fit$summary.random,
        "Slide_ID",
        df$Slide_ID
    )

    gene_eff <- .safe_pull_re_mean(
        fit$summary.random,
        "Gene_ID",
        df$Gene_ID
    ) * df$is_nonNEG

    # ---- mixture DE effect = sum of components ----
    de_small <- .safe_pull_re_mean(
        fit$summary.random,
        "Gene_ID_de_small",
        df$Gene_ID_de_small
    )

    de_large <- .safe_pull_re_mean(
        fit$summary.random,
        "Gene_ID_de_large",
        df$Gene_ID_de_large
    )

    gene_de_eff <- (de_small + de_large) * df$diff

    latent_effs <- list(
        intercept   = intercept,
        spatial_eff = spatial_eff,
        slide_eff   = slide_eff,
        gene_eff    = gene_eff,
        gene_de_eff = gene_de_eff,
        log_lib     = df$log_lib
    )

    size <- fit$summary.hyperpar[
        "size for the nbinomial observations (1/overdispersion)",
        "mean"
    ]

    list(
        fit     = fit,
        eta_full = eta_full,
        size     = size,
        latent   = latent_effs
    )
}