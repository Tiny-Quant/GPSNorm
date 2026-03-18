library(dplyr)
library(tidyr)


#' Canonicalize a SpatialExperiment for benchmarking simulation.
#'
#' Ensures the baseline `spe` has:
#'   (i) a per-gene type label in {NEG, HK, REG},
#'  (ii) a per-AOI binary group label (for DE),
#' (iii) a per-AOI spatial domain/cluster label (for ARI benchmarking).
#'
#' If any are missing, they are constructed from the observed counts and/or
#' spatial coordinates using deterministic rules (given `seed`).
#'
#' @param spe A SpatialExperiment with a 'counts' assay.
#' @param gene_type_col Optional character. If present in rowData(spe), used as gene type.
#' @param group_col Optional character. If present in colData(spe), used as group label.
#' @param domain_col Optional character. If present in colData(spe), used as domain label.
#' @param counts_assay Character. Assay name for counts. Default "counts".
#' @param prop_neg,prop_hk Fractions of genes to label NEG/HK if missing.
#' @param k_domains Integer number of domains (only used if domain labels missing).
#' @param seed Integer seed for reproducibility.
#'
#' @return A list with:
#'   - `spe`: SpatialExperiment with added `rowData(spe)$gene_type` and/or
#'            `colData(spe)$sim_group`, `colData(spe)$sim_domain` as needed.
#'   - `meta`: list describing which fields were created and key summaries.
canonicalize_spe_inputs <- function(
    spe,
    gene_type_col = NULL,
    group_col = NULL,
    domain_col = NULL,
    counts_assay = "counts",
    prop_neg = 0.05,
    prop_hk = 0.05,
    k_domains = 7,
    seed = 1
) {
    # Check for dependencies. 
    for (
        pkg in 
        c("SpatialExperiment", "SummarizedExperiment", "S4Vectors", "Matrix", "stats")
    ) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    set.seed(seed)

    # Validate input data types. 
    if (!methods::is(spe, "SpatialExperiment")) {
        stop("`spe` must be a SpatialExperiment.", call. = FALSE)
    }
    if (!counts_assay %in% SummarizedExperiment::assayNames(spe)) {
        stop("Assay '", counts_assay, "' not found in `spe`.", call. = FALSE)
    }
    mat <- SummarizedExperiment::assay(spe, counts_assay)

    if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
        stop("Counts assay must be a matrix or Matrix.", call. = FALSE)
    }
    if (nrow(mat) < 10 || ncol(mat) < 10) {
        stop(
            "`spe` too small for benchmark simulation (need >=10 genes and >=10 AOIs).", 
            call. = FALSE
        )
    }
    coords <- SpatialExperiment::spatialCoords(spe)

    if (is.null(coords) || ncol(coords) < 2) {
        stop(
            "`spe` must contain spatialCoords with at least 2 columns.", 
            call. = FALSE
        )
    }
    created <- list(gene_type = FALSE, group = FALSE, domain = FALSE)

    # ---- 1) gene types ----
    spe <- .ensure_gene_types(
        spe = spe,
        mat = mat,
        gene_type_col = gene_type_col,
        prop_neg = prop_neg,
        prop_hk = prop_hk
    )
    if (
        is.null(gene_type_col) || 
        !gene_type_col %in% names(SummarizedExperiment::rowData(spe))
    ) {
        created$gene_type <- TRUE
    }

    # ---- 2) group label ----
    spe <- .ensure_group_label(
        spe = spe,
        group_col = group_col
    )
    if (is.null(group_col) || 
        !group_col %in% names(SummarizedExperiment::colData(spe))
    ) {
        created$group <- TRUE
    }

    # ---- 3) domain label ----
    spe <- .ensure_domain_label(
        spe = spe,
        coords = coords,
        domain_col = domain_col,
        k_domains = k_domains
    )
    if (is.null(domain_col) || 
        !domain_col %in% names(SummarizedExperiment::colData(spe))) {
        created$domain <- TRUE
    }

    meta <- list(
        created = created,
        n_gene = nrow(spe),
        n_aoi = ncol(spe),
        gene_type_table = table(
            SummarizedExperiment::rowData(spe)$gene_type, useNA = "ifany"
        ),
        group_table = table(
            SummarizedExperiment::colData(spe)$sim_group, useNA = "ifany"
        ),
        domain_table = table(
            SummarizedExperiment::colData(spe)$sim_domain, useNA = "ifany"
        )
    )

    list(spe = spe, meta = meta)
} 
# Helper functions. 
.ensure_gene_types <- function(
    spe, mat, gene_type_col, prop_neg, prop_hk 
) {
    rd <- SummarizedExperiment::rowData(spe)

    # If provided and present: normalize names to {NEG,HK,REG} and store as gene_type
    if (!is.null(gene_type_col) && gene_type_col %in% names(rd)) {
        gt <- as.character(rd[[gene_type_col]])
        gt[is.na(gt) | gt == ""] <- "REG"
        gt <- toupper(gt)
        # allow common variants; map to canonical
        gt <- dplyr::case_when(
            gt %in% c("NEG", "NEGATIVE") ~ "NEG",
            gt %in% c("HK", "HOUSEKEEPING", "CONTROL") ~ "HK",
            gt %in% c("REG", "ENDOGENOUS") ~ "REG",
            TRUE ~ gt
        )
        if (!all(gt %in% c("NEG", "HK", "REG"))) {
            stop("Existing gene types must map to {NEG, HK, REG}.", call. = FALSE)
        }
        rd$gene_type <- factor(gt, levels = c("NEG", "HK", "REG"))
        SummarizedExperiment::rowData(spe) <- rd
        return(spe)
    }

    # Otherwise: construct strong-control types from counts
    if (!is.numeric(prop_neg) || prop_neg <= 0 || prop_neg >= 1) {
        stop("`prop_neg` must be in (0,1).", call. = FALSE)
    }
    if (!is.numeric(prop_hk)  || prop_hk  <= 0 || prop_hk  >= 1) {
        stop("`prop_hk` must be in (0,1).", call. = FALSE)
    }
    if (prop_neg + prop_hk >= 0.8) {
        stop("`prop_neg + prop_hk` too large.", call. = FALSE)
    }

    # Use library-size adjusted log1p counts for stability
    lib <- Matrix::colSums(mat)
    lib <- pmax(as.numeric(lib), 1)
    sf <- mean(lib) / lib
    mat_norm <- mat %*% Matrix::Diagonal(x = sf)
    #mat_norm <- t(t(mat) / lib) * mean(lib)  # simple global scaling
    logmat <- log1p(mat_norm)

    gene_mean <- Matrix::rowMeans(logmat)
    gene_var  <- .row_var_sparse(logmat) 

    n_gene <- length(gene_mean)
    n_neg <- max(1L, as.integer(round(n_gene * prop_neg)))
    n_hk  <- max(1L, as.integer(round(n_gene * prop_hk)))

    # NEG: lowest mean expression (strongly low biology)
    neg_idx <- order(gene_mean, decreasing = FALSE)[seq_len(n_neg)]

    # HK: stable across AOIs (low variance), among remaining genes
    remain <- setdiff(seq_len(n_gene), neg_idx)
    hk_idx <- remain[order(gene_var[remain], decreasing = FALSE)][
        seq_len(min(n_hk, length(remain)))
    ]
    gt <- rep("REG", n_gene)
    gt[neg_idx] <- "NEG"
    gt[hk_idx]  <- "HK"

    rd$gene_type <- factor(gt, levels = c("NEG", "HK", "REG"))
    SummarizedExperiment::rowData(spe) <- rd

    spe
}
# .row_var_dense <- function(x) {
#   # Returns row-wise variance for matrix or dgCMatrix via dense ops in blocks.
#   # Replace with a sparse-native version later if needed.
#   if (inherits(x, "Matrix")) {
#     x <- as.matrix(x)
#   }
#   m <- rowMeans(x)

#   rowMeans((x - m)^2)
# }
.row_var_sparse <- function(mat) {
  # mat: dgCMatrix (genes x AOIs)
  stopifnot(inherits(mat, "Matrix"))

  n <- ncol(mat)
  if (n <= 1) {
    return(rep(NA_real_, nrow(mat)))
  }

  rs1 <- Matrix::rowSums(mat)           # Σ x
  rs2 <- Matrix::rowSums(mat ^ 2)       # Σ x^2

  (rs2 - (rs1^2) / n) / (n - 1)
}
.ensure_group_label <- function(
    spe, group_col
) {
    cd <- SummarizedExperiment::colData(spe)

    # Format to binary if grouping is available. 
    if (!is.null(group_col) && group_col %in% names(cd)) {
        g <- cd[[group_col]]
        # coerce to 0/1 numeric, keep as sim_group for downstream consistency
        if (is.logical(g)) g <- as.integer(g)
        if (is.factor(g) || is.character(g)) g <- as.integer(as.factor(g)) - 1L
        g <- as.integer(g)
        if (!all(g %in% c(0L, 1L))) {
            stop("`group_col` must be binary.", call. = FALSE)
        }
        cd$sim_group <- g
        SummarizedExperiment::colData(spe) <- cd
        return(spe)
    }

    # Simulate if missing. 
    cd$sim_group <- stats::rbinom(ncol(spe), 1, 0.5) # TODO: Class imbalance. 
    SummarizedExperiment::colData(spe) <- cd
    spe
}
.ensure_domain_label <- function(
    spe, coords, domain_col, k_domains
) {
    cd <- SummarizedExperiment::colData(spe)

    # Format if column provided. 
    if (!is.null(domain_col) && domain_col %in% names(cd)) {
        z <- cd[[domain_col]]
        if (is.factor(z) || is.character(z)) z <- as.integer(as.factor(z))
        z <- as.integer(z)
        #if (anyNA(z)) stop("`domain_col` contains NA.", call. = FALSE)
        cd$sim_domain <- z
        SummarizedExperiment::colData(spe) <- cd
        return(spe)
    }

    # Simulate if missing. 
    # TODO: Simulate NA/Unlabeled domains. 
    if (!is.numeric(k_domains) || length(k_domains) != 1L || is.na(k_domains) || k_domains < 2) {
        stop("`k_domains` must be an integer >= 2.", call. = FALSE)
    }
    k_domains <- as.integer(k_domains)
    k_domains <- min(k_domains, ncol(spe) - 1L)

    xy <- coords[, 1:2, drop = FALSE]
    xy <- apply(xy, 2, as.numeric)
    xy <- scale(xy)

    # TODO: Allow other domain simulation methods. 
    km <- stats::kmeans(xy, centers = k_domains, nstart = 10) 
    cd$sim_domain <- as.integer(km$cluster)
    SummarizedExperiment::colData(spe) <- cd
    spe
}

#' 
#' @export 
simulate_spe_benchmark <- function(
    spe,
    n_slide = 3L,
    counts_assay = "counts",
    gene_type_col = NULL,
    group_col = NULL,
    domain_col = NULL,
    prop_neg = 0.05,
    prop_hk = 0.05,
    k_domains = 5L,
    de_prop = 0.10,
    de_logfc_log_mean = 0.5,
    dom_prop = 0.10,
    dom_sd = 0.35,
    svg_prop = 0.10,
    svg_sd = 0.35,
    de_svg_overlap = 0.30,
    de_dom_overlap = 0.30,
    tech_spatial_sd = 0.4,
    tech_rbf_knots = 10L,
    tech_rbf_scale = 0.16,
    slide_rho = 0.6,
    slide_lib_sd = 0.6, 
    slide_group_probs = NULL, 
    zi_prob = 0.01,
    spot_prop = NULL,
    gene_prop = NULL,
    n_spot_sub = NULL,
    n_gene_sub = NULL,
    bg_prop = 0.05,
    bg_only_controls = FALSE, 

    group_confounded = FALSE,
    group_conf_rho = 0.0,            
    group_conf_mode = c("logit", "hard"),
    group_conf_center_by_slide = TRUE,     
    group_conf_target_p = 0.5,       
    seed = 1
) {

    for (pkg in c("Matrix", "SummarizedExperiment",
                  "SpatialExperiment", "S4Vectors", "stats")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    set.seed(seed)

    if (!is.numeric(bg_prop) || length(bg_prop) != 1L || !is.finite(bg_prop) ||
        bg_prop < 0 || bg_prop >= 1) {
        stop("`bg_prop` must be in [0, 1).", call. = FALSE)
    }

    if (!is.numeric(slide_lib_sd) || length(slide_lib_sd) != 1L ||
        !is.finite(slide_lib_sd) || slide_lib_sd < 0) {
        stop("`slide_lib_sd` must be a non-negative finite scalar.", call. = FALSE)
    }

    if (is.null(slide_group_probs)) {
        #slide_group_probs <- rep(NA_real_, n_slide)
        slide_group_probs <- rep(1 / n_slide, n_slide)
    } else {
        if (!is.numeric(slide_group_probs) ||
            length(slide_group_probs) != n_slide ||
            any(!is.finite(slide_group_probs)) ||
            any(slide_group_probs <= 0 | slide_group_probs >= 1)) {
            stop("`slide_group_probs` must be numeric in (0,1) with length n_slide.",
                call. = FALSE)
        }
    }

    if (!is.logical(group_confounded) || length(group_confounded) != 1L) {
        stop("`group_confounded` must be TRUE/FALSE.", call. = FALSE)
    }

    if (!is.numeric(group_conf_rho) || length(group_conf_rho) != 1L ||
        !is.finite(group_conf_rho) || group_conf_rho < 0) {
        stop("`group_conf_rho` must be a finite non-negative scalar.", call. = FALSE)
    }

    if (!is.numeric(group_conf_target_p) || length(group_conf_target_p) != 1L ||
        !is.finite(group_conf_target_p) || group_conf_target_p <= 0 || group_conf_target_p >= 1) {
        stop("`group_conf_target_p` must be in (0,1).", call. = FALSE)
    }

    if (!is.logical(group_conf_center_by_slide) || length(group_conf_center_by_slide) != 1L) {
        stop("`group_conf_center_by_slide` must be TRUE/FALSE.", call. = FALSE)
    }

    # ------------------------------------------------------------------
    # Canonicalize inputs
    # ------------------------------------------------------------------
    canon <- canonicalize_spe_inputs(
        spe = spe,
        gene_type_col = gene_type_col,
        group_col = group_col,
        domain_col = domain_col,
        counts_assay = counts_assay,
        prop_neg = prop_neg,
        prop_hk = prop_hk,
        k_domains = k_domains,
        seed = seed
    )
    spe0 <- canon$spe

    # ------------------------------------------------------------------
    # Subsample spots
    # ------------------------------------------------------------------
    if (is.null(n_spot_sub) && !is.null(spot_prop)) {
        n_spot_sub <- max(1L, as.integer(round(ncol(spe0) * spot_prop)))
    }
    if (!is.null(n_spot_sub)) {
        n_spot_sub <- min(as.integer(n_spot_sub), ncol(spe0))
        keep_cols <- sort(sample(seq_len(ncol(spe0)), n_spot_sub))
        spe0 <- spe0[, keep_cols]
    }

    # ------------------------------------------------------------------
    # Subsample genes (stratified by gene_type) — FIXED VERSION
    # ------------------------------------------------------------------
    if (is.null(n_gene_sub) && !is.null(gene_prop)) {
        n_gene_sub <- max(1L, as.integer(round(nrow(spe0) * gene_prop)))
    }

    if (!is.null(n_gene_sub)) {
        n_gene_sub <- min(as.integer(n_gene_sub), nrow(spe0))

        gt <- as.character(SummarizedExperiment::rowData(spe0)$gene_type)
        tab <- table(gt)
        prop <- as.numeric(tab) / sum(tab)

        tgt <- floor(prop * n_gene_sub)
        names(tgt) <- names(tab)

        rem <- n_gene_sub - sum(tgt)
        if (rem > 0L) {
            frac <- prop * n_gene_sub - tgt
            add <- order(frac, decreasing = TRUE)[seq_len(rem)]
            tgt[add] <- tgt[add] + 1L
        }

        keep_rows <- integer(0)
        for (lvl in names(tgt)) {
            idx <- which(gt == lvl)
            nk <- min(length(idx), tgt[lvl])
            if (nk > 0L) {
                keep_rows <- c(keep_rows, sample(idx, nk))
            }
        }

        if (length(keep_rows) < n_gene_sub) {
            pool <- setdiff(seq_len(nrow(spe0)), keep_rows)
            keep_rows <- c(
                keep_rows,
                sample(pool, n_gene_sub - length(keep_rows))
            )
        }

        spe0 <- spe0[sort(keep_rows), ]
    }

    # ------------------------------------------------------------------
    # NEW: scale factor for gene subsampling
    # ------------------------------------------------------------------
    gene_scale <- 1
    if (!is.null(gene_prop)) {
        gene_scale <- 1 / gene_prop
    } else if (!is.null(n_gene_sub)) {
        gene_scale <- nrow(SummarizedExperiment::assay(spe, counts_assay)) /
                      n_gene_sub
    }

    # ------------------------------------------------------------------
    # Extract baseline data
    # ------------------------------------------------------------------
    mat0 <- SummarizedExperiment::assay(spe0, counts_assay)
    if (!inherits(mat0, "dgCMatrix")) {
        mat0 <- Matrix::Matrix(mat0, sparse = TRUE)
    }

    coords0 <- .as_numeric_matrix(
        SpatialExperiment::spatialCoords(spe0)
    )

    n_gene <- nrow(mat0)
    n_spot <- ncol(mat0)

    rd0 <- SummarizedExperiment::rowData(spe0)
    gene_type <- as.character(rd0$gene_type)
    is_neg <- gene_type == "NEG"
    is_hk  <- gene_type == "HK"
    is_reg <- gene_type == "REG"
    reg_ids <- which(is_reg)

    if (length(reg_ids) < 2L) {
        stop("Need at least 2 REG genes to simulate.", call. = FALSE)
    }

    # ------------------------------------------------------------------
    # Baseline gene composition and dispersion
    # ------------------------------------------------------------------
    lambda_g <- Matrix::rowMeans(mat0)
    lambda_g <- pmax(as.numeric(lambda_g), 1e-8)
    lambda_g <- .winsorize_log(lambda_g, 0.001, 0.999)
    p_g <- lambda_g / sum(lambda_g)

    size_g <- .nb_size_from_sparse(mat0)
    size_g <- pmin(pmax(size_g, 0.5), 200)

    # ---- NEW: background distribution over genes ----
    # Default: uniform background across genes (strongest anti-sparsity)
    if (!isTRUE(bg_only_controls)) {
        u_g <- rep(1 / n_gene, n_gene)
    } else {
        ctrl <- which(is_neg | is_hk)
        if (length(ctrl) == 0L) {
            u_g <- rep(1 / n_gene, n_gene)
        } else {
            u_g <- numeric(n_gene)
            u_g[ctrl] <- 1 / length(ctrl)
        }
    }

    # ------------------------------------------------------------------
    # Choose DE / DOM / SVG genes (ALLOW OVERLAP WITH DE)
    # ------------------------------------------------------------------
    if (!is.numeric(de_svg_overlap) || length(de_svg_overlap) != 1L ||
        !is.finite(de_svg_overlap) || de_svg_overlap < 0 || de_svg_overlap > 1) {
    stop("`de_svg_overlap` must be in [0, 1].", call. = FALSE)
    }
    if (!is.numeric(de_dom_overlap) || length(de_dom_overlap) != 1L ||
        !is.finite(de_dom_overlap) || de_dom_overlap < 0 || de_dom_overlap > 1) {
    stop("`de_dom_overlap` must be in [0, 1].", call. = FALSE)
    }

    # 1) DE genes
    de_ids <- .sample_subset(reg_ids, de_prop)

    # Targets (marginal) among REG genes
    n_reg <- length(reg_ids)
    tgt_svg <- max(0L, as.integer(round(svg_prop * n_reg)))
    tgt_dom <- max(0L, as.integer(round(dom_prop * n_reg)))

    # 2) Forced overlap subsets from DE
    n_de <- length(de_ids)
    n_de_svg <- min(n_de, max(0L, as.integer(round(de_svg_overlap * tgt_svg))))
    n_de_dom <- min(n_de, max(0L, as.integer(round(de_dom_overlap * tgt_dom))))

    de_svg_ids <- if (n_de_svg > 0L) sample(de_ids, n_de_svg) else integer(0)
    de_dom_ids <- if (n_de_dom > 0L) sample(de_ids, n_de_dom) else integer(0)

    # 3) Fill remaining SVG/DOM genes from non-DE pool to hit marginal targets
    pool_nonde <- setdiff(reg_ids, de_ids)

    need_svg <- max(0L, tgt_svg - length(de_svg_ids))
    need_dom <- max(0L, tgt_dom - length(de_dom_ids))

    extra_svg_ids <- if (need_svg > 0L) {
    sample(pool_nonde, min(need_svg, length(pool_nonde)))
    } else integer(0)

    pool_nonde2 <- setdiff(pool_nonde, extra_svg_ids)

    extra_dom_ids <- if (need_dom > 0L) {
    sample(pool_nonde2, min(need_dom, length(pool_nonde2)))
    } else integer(0)

    # 4) Final sets (now can overlap with DE)
    svg_ids <- sort(unique(c(de_svg_ids, extra_svg_ids)))
    dom_ids <- sort(unique(c(de_dom_ids, extra_dom_ids)))

    true_logfc <- numeric(n_gene)
    if (length(de_ids) > 0L) {
        sign <- sample(c(-1, 1), length(de_ids), replace = TRUE)
        mag  <- abs(
            rnorm(length(de_ids), mean = de_logfc_log_mean, sd = 0.25)
        )
        true_logfc[de_ids] <- sign * mag
    }


    cd0 <- SummarizedExperiment::colData(spe0)
    dom0 <- cd0$sim_domain
    k_dom <- max(as.integer(stats::na.omit(dom0)), na.rm = TRUE)
    if (!is.finite(k_dom) || k_dom < 1L) k_dom <- k_domains

    gamma_gd <- matrix(0, n_gene, k_dom)
    for (g in dom_ids) {
        v <- rnorm(k_dom, 0, dom_sd)
        gamma_gd[g, ] <- v - mean(v)
    }


    truth_gene <- S4Vectors::DataFrame(
        gene_type = factor(gene_type, levels = c("NEG", "HK", "REG")),
        is_de = as.integer(seq_len(n_gene) %in% de_ids),
        true_logfc = true_logfc,
        is_svg = as.integer(seq_len(n_gene) %in% svg_ids),
        is_dom = as.integer(seq_len(n_gene) %in% dom_ids),
        lambda_g = lambda_g,
        size_g = size_g
    )
    rownames(truth_gene) <- rownames(spe0)

    # ------------------------------------------------------------------
    # Build slides / spot metadata
    # ------------------------------------------------------------------
    base_group <- as.integer(cd0$sim_group)
    p_group <- mean(base_group == 1, na.rm = TRUE)
    p_group <- ifelse(is.finite(p_group), p_group, 0.5)

    slides <- vector("list", n_slide)
    for (s in seq_len(n_slide)) {
        tr <- .coord_transform_params(seed + 1000L + s)
        xy <- .apply_coord_transform(coords0, tr)

        p_s <- slide_group_probs[s]
        if (!is.finite(p_s)) p_s <- p_group
        grp <- rbinom(n_spot, 1, p_s)

        slides[[s]] <- S4Vectors::DataFrame(
            spot_id = paste0("spot_", seq_len(n_spot)),
            Slide_ID = s,
            sim_group = grp,
            sim_domain = dom0,
            x = xy[, 1],
            y = xy[, 2]
        )
    }

    truth_spot <- do.call(S4Vectors::rbind, slides)
    coords_all <- cbind(truth_spot$x, truth_spot$y)
    n_spot_all <- nrow(truth_spot)

    tech_field_raw <- rep(0, n_spot_all)
    tech_field_ctr <- rep(0, n_spot_all)

    grp_all <- as.integer(truth_spot$sim_group)
    dom_all <- as.integer(truth_spot$sim_domain)
    ok_dom <- !is.na(dom_all) & dom_all >= 1 & dom_all <= k_dom

    # ------------------------------------------------------------------
    # Library size calibration
    # ------------------------------------------------------------------
    lib0 <- as.numeric(Matrix::colSums(mat0))
    lib0 <- lib0[lib0 > 0 & is.finite(lib0)]
    lo <- quantile(lib0, 0.05, names = FALSE)
    lib0_trim <- lib0[lib0 >= lo]
    if (length(lib0_trim) < 10L) lib0_trim <- lib0

    L_min <- max(1, quantile(lib0_trim, 0.01, names = FALSE))
    L_all <- pmax(sample(lib0_trim, n_spot_all, TRUE), L_min)

    # ------------------------------------------------------------------
    # NEW: compensate library size for gene subsampling
    # ------------------------------------------------------------------
    L_all <- L_all * gene_scale
    L_all <- pmax(L_all, L_min)
    slide_alpha <- rnorm(n_slide, mean = 0, sd = slide_lib_sd)

    idx0 <- 0L
    for (s in seq_len(n_slide)) {
        rows <- idx0 + seq_len(n_spot)
        idx0 <- idx0 + n_spot
        L_all[rows] <- L_all[rows] * exp(slide_alpha[s])
    }
    L_all <- pmax(L_all, L_min)

    # ------------------------------------------------------------------
    # Spatial exposure field (library size only)
    # ------------------------------------------------------------------
    if (tech_spatial_sd > 0) {
        w0 <- rnorm(tech_rbf_knots)
        idx0 <- 0L
        for (s in seq_len(n_slide)) {
            rows <- idx0 + seq_len(n_spot)
            idx0 <- idx0 + n_spot
            B <- .rbf_basis(coords_all[rows, ], tech_rbf_knots,
                            tech_rbf_scale, seed + 2000L + s)
            e <- rnorm(length(w0))
            ws <- slide_rho * w0 + sqrt(1 - slide_rho^2) * e
            f_raw <- .zscore(as.numeric(B %*% ws)) * tech_spatial_sd
            f_ctr <- f_raw - mean(f_raw)
            L_all[rows] <- L_all[rows] * exp(f_ctr)

            tech_field_raw[rows] <- f_raw
            tech_field_ctr[rows] <- f_ctr
        }
        L_all <- pmax(L_all, L_min)
    }

    truth_spot$tech_exposure_raw <- tech_field_raw
    truth_spot$tech_exposure_ctr <- tech_field_ctr
    truth_spot$L_all <- L_all 
    truth_spot$log_L_all <- log(L_all)

    # ------------------------------------------------------------
    # Guardrail: technical spatial field must be centered per slide
    # ------------------------------------------------------------
    slide_means <- tapply(
        tech_field_ctr,
        truth_spot$Slide_ID,
        mean
    )

    tol <- 1e-8

    if (any(abs(slide_means) > tol)) {
        stop(
            "Technical spatial field is not centered within slide.\n",
            "Max |mean| = ",
            signif(max(abs(slide_means)), 4), "\n",
            "This indicates a bug in spatial field indexing or assignment.",
            call. = FALSE
        )
    }

    # ------------------------------------------------------------------
    # NEW (Scenario G): override sim_group using tech field
    # ------------------------------------------------------------------
    if (isTRUE(group_confounded) && group_conf_rho > 0) {

        # choose confounder: centered-by-slide technical field is safest
        z <- if (isTRUE(group_conf_center_by_slide)) {
            truth_spot$tech_exposure_ctr
        } else {
            truth_spot$tech_exposure_raw
        }

        # per-slide intercept calibration (logit mode) so marginal p ~= target
        if (group_conf_mode == "logit") {
            # For each slide, pick alpha_s so mean(sigmoid(alpha_s + rho*z)) ~= target.
            # We do a quick 1D root-find per slide on [-20, 20].
            slide_ids <- unique(truth_spot$Slide_ID)

            alpha_s <- numeric(max(slide_ids))
            for (s in slide_ids) {
                idx <- which(truth_spot$Slide_ID == s)
                zz <- z[idx]

                f_obj <- function(a) mean(stats::plogis(a + group_conf_rho * zz)) - group_conf_target_p
                # robust bracket; if tech field is extreme, this still works
                alpha_s[s] <- tryCatch(
                    stats::uniroot(f_obj, interval = c(-20, 20))$root,
                    error = function(e) stats::qlogis(group_conf_target_p) # fallback
                )
            }

            p_i <- stats::plogis(alpha_s[truth_spot$Slide_ID] + group_conf_rho * z)
            truth_spot$sim_group <- stats::rbinom(n_spot_all, 1, p_i)

        } else { # hard mode
            # Group = 1 if tech field above its within-slide median (keeps ~50/50 per slide)
            truth_spot$sim_group <- 0L
            for (s in unique(truth_spot$Slide_ID)) {
                idx <- which(truth_spot$Slide_ID == s)
                thr <- stats::median(z[idx], na.rm = TRUE)
                truth_spot$sim_group[idx] <- as.integer(z[idx] > thr)
            }
        }
    }

    # refresh derived vectors after possible override
    grp_all <- as.integer(truth_spot$sim_group)
    dom_all <- as.integer(truth_spot$sim_domain)
    ok_dom <- !is.na(dom_all) & dom_all >= 1 & dom_all <= k_dom

    # ------------------------------------------------------------------
    # SVG fields
    # ------------------------------------------------------------------
    svg_term <- matrix(0, n_spot_all, length(svg_ids))
    svg_map <- rep(NA_integer_, n_gene)

    if (length(svg_ids) > 0L) {
        svg_map[svg_ids] <- seq_along(svg_ids)
        w_svg0 <- matrix(rnorm(tech_rbf_knots * length(svg_ids)),
                         tech_rbf_knots)

        idx0 <- 0L
        for (s in seq_len(n_slide)) {
            rows <- idx0 + seq_len(n_spot)
            idx0 <- idx0 + n_spot
            B <- .rbf_basis(coords_all[rows, ], tech_rbf_knots,
                            tech_rbf_scale, seed + 3000L + s)
            e <- matrix(rnorm(length(w_svg0)), nrow(w_svg0))
            ws <- slide_rho * w_svg0 + sqrt(1 - slide_rho^2) * e
            svg_term[rows, ] <- apply(B %*% ws, 2, .zscore) * svg_sd
        }
    }
    # Extract SVG effect sizes. 
    true_svg_var <- numeric(n_gene)

    if (length(svg_ids) > 0L) {
        for (g in svg_ids) {
            k <- svg_map[g]
            true_svg_var[g] <- var(svg_term[, k])
        }
    }

    truth_gene$true_svg_var <- true_svg_var

    # ------------------------------------------------------------------
    # OPTIMIZED NORMALIZATION (build expB_reg and denom once)
    # ------------------------------------------------------------------
    reg_ids <- which(is_reg)
    n_reg <- length(reg_ids)

    B_reg <- matrix(0, n_reg, n_spot_all)

    for (k in seq_len(n_reg)) {
        g <- reg_ids[k]
        b <- numeric(n_spot_all)

        if (true_logfc[g] != 0) {
            b <- b + true_logfc[g] * grp_all
        }
        if (g %in% dom_ids) {
            add <- numeric(n_spot_all)
            add[ok_dom] <- gamma_gd[g, dom_all[ok_dom]]
            b <- b + add
        }
        if (!is.na(svg_map[g])) {
            b <- b + svg_term[, svg_map[g]]
        }
        B_reg[k, ] <- b
    }

    expB_reg <- exp(B_reg)
    base_mass <- sum(p_g[!is_reg])
    denom <- base_mass + colSums(expB_reg * p_g[reg_ids])
    denom <- pmax(denom, 1e-12)

    # ------------------------------------------------------------------
    # Sample counts (with background mixture)
    # ------------------------------------------------------------------
    ii <- integer(0)
    jj <- integer(0)
    xx <- numeric(0)

    eps <- bg_prop  # naming consistent with math description

    for (g in seq_len(n_gene)) {

        # current composition-based proportion pi_gi
        if (!is_reg[g]) {
            pi <- (p_g[g] / denom)
        } else {
            k <- match(g, reg_ids)
            pi <- (p_g[g] * expB_reg[k, ] / denom)
        }

        # ---- NEW: mix in uniform technical background ----
        prop_gi <- (1 - eps) * pi + eps * u_g[g]
        prop_gi <- pmax(prop_gi, 1e-15)

        mu <- L_all * prop_gi
        mu <- pmax(mu, 1e-12)

        y <- rnbinom(n_spot_all, size = size_g[g], mu = mu)

        if (zi_prob > 0) {
            z <- rbinom(n_spot_all, 1, zi_prob)
            y[z == 1] <- 0L
        }

        nz <- which(y > 0L)
        if (length(nz) > 0L) {
            ii <- c(ii, rep.int(g, length(nz)))
            jj <- c(jj, nz)
            xx <- c(xx, as.numeric(y[nz]))
        }
    }

    mat_sim <- Matrix::sparseMatrix(
        i = ii, j = jj, x = xx,
        dims = c(n_gene, n_spot_all)
    )
    rownames(mat_sim) <- rownames(spe0)
    colnames(mat_sim) <- paste0("spot_", seq_len(n_spot_all))

    spe_out <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = mat_sim),
        rowData = truth_gene,
        colData = truth_spot,
        spatialCoords = coords_all
    )

    # Hot Patch -- fix later. 
    cd <- SummarizedExperiment::colData(spe_out)
    rd <- SummarizedExperiment::rowData(spe_out)

    # Ensure identifiers exist for the INLA model
    if (!"AOI_ID" %in% colnames(cd)) {
        cd$AOI_ID <- seq_len(ncol(spe_out))
    }
    if (!"Slide_ID" %in% colnames(cd)) {
        if ("slide_id" %in% colnames(cd)) {
            cd$Slide_ID <- as.integer(cd$slide_id)
        } else {
            cd$Slide_ID <- 1L
        }
    }
    SummarizedExperiment::colData(spe_out) <- cd

    if (!"Gene_ID" %in% colnames(rd)) {
        rd$Gene_ID <- seq_len(nrow(spe_out))
    }
    if (!"gene_type" %in% colnames(rd)) {
        stop("`rowData(spe_out)` must contain `gene_type` for NEG/HK/REG labels.",
            call. = FALSE)
    }
    SummarizedExperiment::rowData(spe_out) <- rd

    list(
        spe = spe_out,
        truth_gene = truth_gene,
        truth_spot = truth_spot, 
        gamma_gd = gamma_gd,
        svg_term = svg_term,
        svg_map = svg_map
    )
}
# Helper functions.
.as_numeric_matrix <- function(x) {
    x <- as.matrix(x)
    x[, 1] <- as.numeric(x[, 1])
    x[, 2] <- as.numeric(x[, 2])
    x
}
.nb_size_from_sparse <- function(mat) {
    # Estimates per-gene NB "size" using Var(Y)=mu + mu^2/size.
    if (!inherits(mat, "dgCMatrix")) {
        mat <- Matrix::Matrix(mat, sparse = TRUE)
    }
    mu <- Matrix::rowMeans(mat)

    mat2 <- mat
    mat2@x <- mat2@x^2
    mu2 <- Matrix::rowMeans(mat2)

    v <- pmax(as.numeric(mu2 - mu^2), 0)
    mu <- pmax(as.numeric(mu), 1e-8)

    size <- (mu^2) / pmax(v - mu, 1e-6)
    as.numeric(size)
}
.sample_subset <- function(ids, prop) {
    ids <- as.integer(ids)
    if (length(ids) == 0L) return(integer(0))
    if (!is.finite(prop) || prop <= 0) return(integer(0))
    n <- max(1L, as.integer(round(length(ids) * prop)))
    sample(ids, size = min(n, length(ids)), replace = FALSE)
}
.coord_transform_params <- function(seed) {
    set.seed(seed)
    list(
        theta = stats::runif(1, 0, 2 * pi),
        refl_x = stats::rbinom(1, 1, 0.5) == 1,
        refl_y = stats::rbinom(1, 1, 0.5) == 1,
        shift = stats::rnorm(2, mean = 0, sd = 0.25)
    )
}
.apply_coord_transform <- function(xy, tr) {
    xy <- .as_numeric_matrix(xy)
    c0 <- colMeans(xy)
    z <- sweep(xy, 2, c0, "-")

    if (isTRUE(tr$refl_x)) z[, 1] <- -z[, 1]
    if (isTRUE(tr$refl_y)) z[, 2] <- -z[, 2]

    ct <- cos(tr$theta)
    st <- sin(tr$theta)
    R <- matrix(c(ct, -st, st, ct), nrow = 2, byrow = TRUE)

    z <- z %*% t(R)
    z <- sweep(z, 2, c0, "+")
    z <- sweep(z, 2, tr$shift, "+")
    z
}
.rbf_basis <- function(xy, m = 12L, scale = 0.2, seed = 1) {
    #xy <- .as_numeric_matrix(xy)
    xy <- .scale_coords01(xy)
    n <- nrow(xy)
    if (n < 2L) {
        return(matrix(1, nrow = n, ncol = 1))
    }

    m <- as.integer(m)
    m <- max(1L, min(m, n))

    set.seed(seed)
    idx <- sample(seq_len(n), size = m, replace = FALSE)
    knots <- xy[idx, , drop = FALSE]

    # squared distances to knots
    d2 <- matrix(0, nrow = n, ncol = m)
    for (j in seq_len(m)) {
        dx <- xy[, 1] - knots[j, 1]
        dy <- xy[, 2] - knots[j, 2]
        d2[, j] <- dx * dx + dy * dy
    }

    B <- exp(-d2 / (2 * scale * scale))
    B
}
.zscore <- function(x) {
    x <- as.numeric(x)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
    (x - mean(x)) / s
}
.winsorize_log <- function(x, p_lo = 0.001, p_hi = 0.999) {
    lx <- log(pmax(as.numeric(x), 1e-12))
    qs <- stats::quantile(
        lx, probs = c(p_lo, p_hi), na.rm = TRUE, names = FALSE
    )
    exp(pmin(pmax(lx, qs[1]), qs[2]))
}
.scale_coords01 <- function(xy) {
    xy <- .as_numeric_matrix(xy)

    for (j in 1:2) {
        r <- range(xy[, j], finite = TRUE)
        den <- r[2] - r[1]
        if (!is.finite(den) || den <= 0) den <- 1
        xy[, j] <- (xy[, j] - r[1]) / den
    }
    xy
}
