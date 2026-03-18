#'
#' @export
evaluate_de_performance_limma <- function(
    spe,
    norm_assay = "logcounts",
    group_col = "sim_group",
    fdr_level = 0.05,
    top_k = NULL,
    return_table = FALSE
) {
    for (pkg in c("SummarizedExperiment", "limma")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    validate_spe_contract(spe, assay = norm_assay)
    rd <- SummarizedExperiment::rowData(spe)
    cd <- SummarizedExperiment::colData(spe)

    required_gene <- c("is_de", "true_logfc")
    miss_gene <- setdiff(required_gene, colnames(rd))
    if (length(miss_gene) > 0L) {
        stop(
        "rowData(spe) is missing: ",
        paste(miss_gene, collapse = ", "),
        call. = FALSE
        )
    }

    if (!group_col %in% colnames(cd)) {
        stop("colData(spe) is missing group column '",
            group_col, "'.", call. = FALSE)
    }

    y <- SummarizedExperiment::assay(spe, norm_assay)
    if (any(!is.finite(y))) {
        stop("`norm_assay` contains non-finite values.", call. = FALSE)
    }

    group_vec <- as.factor(cd[[group_col]])
    if (!all(levels(group_vec) %in% c("0", "1"))) {
        stop("`group_col` must be binary 0/1.", call. = FALSE)
    }

    if (!"lambda_g" %in% colnames(rd)) {
        stop("rowData(spe) is missing 'lambda_g' (needed for calibration strata).",
            call. = FALSE)
    }

    design <- stats::model.matrix(~ group_vec)

    fit <- limma::lmFit(y, design)
    fit <- limma::eBayes(fit)

    stats_tbl <- limma::topTable(
        fit,
        coef = 2,
        number = nrow(y),
        sort.by = "none",
        confint = 0.95
    )

    logfc_hat <- stats_tbl$logFC
    pval      <- stats_tbl$P.Value
    logfc_025_quant <- stats_tbl$CI.L
    logfc_975_quant <- stats_tbl$CI.R

    # -----------------------------
    # Restrict to REG genes
    # -----------------------------
    gene_type <- if ("gene_type" %in% colnames(rd)) {
        as.character(rd$gene_type)
    } else {
        rep("REG", length(logfc_hat))
    }
    keep <- gene_type == "REG"

    truth_de     <- as.numeric(rd$is_de)[keep]
    truth_logfc  <- as.numeric(rd$true_logfc)[keep]
    logfc_eval   <- as.numeric(logfc_hat[keep])
    pval_eval    <- as.numeric(pval[keep])
    logfc_025_quant_eval <- as.numeric(logfc_025_quant[keep])
    logfc_975_quant_eval <- as.numeric(logfc_975_quant[keep])

    # Evidence score (limma-native)
    score <- -log10(pval_eval + .Machine$double.xmin)

    # Check for classes. 
    n_pos <- sum(truth_de == 1, na.rm = TRUE)
    n_neg <- sum(truth_de == 0, na.rm = TRUE)
    has_two_classes <- n_pos > 0L && n_neg > 0L

    # -----------------------------
    # Effect size metrics
    # -----------------------------
    mse  <- mean((logfc_eval - truth_logfc)^2)
    bias <- mean(logfc_eval - truth_logfc)
    spearman_r <- suppressWarnings(
        stats::cor(logfc_eval, truth_logfc, method = "spearman")
    )
    coverage <- mean(
        logfc_025_quant_eval <= truth_logfc &
        truth_logfc <= logfc_975_quant_eval
    )
    width = mean(
        logfc_975_quant_eval - logfc_025_quant_eval
    )

    # -----------------------------
    # Ranking metrics
    # -----------------------------
    auc <- NA_real_
    ap  <- NA_real_

    if (has_two_classes && requireNamespace("pROC", quietly = TRUE)) {
            roc_obj <- pROC::roc(
                response = truth_de,
                predictor = score,
                direction = "auto",
                quiet = TRUE
            )
            auc <- as.numeric(pROC::auc(roc_obj))
        }

    if (has_two_classes && requireNamespace("PRROC", quietly = TRUE)) {
        fg <- score[truth_de == 1]
        bg <- score[truth_de == 0]
        if (length(fg) > 0L && length(bg) > 0L) {
            ap <- PRROC::pr.curve(
                scores.class0 = fg,
                scores.class1 = bg,
                curve = FALSE
            )$auc.integral
        }
    }

    # -----------------------------
    # Threshold-based metrics
    # -----------------------------
    adj_p <- stats::p.adjust(pval_eval, method = "BH")

    ord <- order(adj_p)
    ranked_truth <- truth_de[ord]

    cum_fp <- cumsum(ranked_truth == 0)
    cum_tp <- cumsum(ranked_truth == 1)

    fdr <- cum_fp / pmax(cum_fp + cum_tp, 1)
    tpr <- if (n_pos > 0L) cum_tp / n_pos else rep(NA_real_, length(cum_tp))

    tpr_at_fdr <- NA_real_
    idx <- which(fdr <= fdr_level)
    if (length(idx) > 0L) {
        tpr_at_fdr <- max(tpr[idx])
    }

    # -----------------------------
    # Top-K metrics
    # -----------------------------
    if (is.null(top_k)) {
        # top_k <- sum(truth_de == 1)
        top_k <- max(10L, as.integer(ceiling(0.05 * length(truth_de))))
    }

    k <- min(top_k, length(truth_de))
    topk_truth <- ranked_truth[seq_len(k)]

    precision_k <- mean(topk_truth == 1)
    recall_k <- if (n_pos > 0L) {
        sum(topk_truth == 1) / n_pos
        } else {
            NA_real_
        }

    diagnostics <- list(
        logfc_hat = logfc_hat,
        p_value = pval, 
        n_genes = length(truth_de),
        n_de = sum(truth_de),
        score_sd = sd(score, na.rm = TRUE),
        unique_scores = length(unique(round(score, 6))),
        min_fdr = min(fdr, na.rm = TRUE)
    )

    # ---------------------
    # Gene-wise metrics
    # ---------------------
    de_table <- NULL
    if (isTRUE(return_table)) {
        de_table <- data.frame(
            truth_de    = truth_de,
            truth_logfc = truth_logfc,
            logfc_hat   = logfc_eval,
            lower       = logfc_025_quant_eval,
            upper       = logfc_975_quant_eval,
            lambda_g    = as.numeric(rd$lambda_g)[keep],
            gene_type   = gene_type[keep]
        )
    }

    # -----------------------------
    # Output
    # -----------------------------
    list(
        # effect recovery
        mse = mse,
        bias = bias,
        spearman_r = spearman_r,
        coverage = coverage, 
        width = width, 

        # ranking
        auc = auc,
        average_precision = ap,

        # decision quality
        tpr_at_fdr = tpr_at_fdr,

        # top-k
        precision_k = precision_k,
        recall_k = recall_k,

        diagnostics = diagnostics, 
        
        de_table = de_table
    )
}

#'
#' @export 
evaluate_de_performance_bayes <- function(
    spe,
    fit_inla,
    gene_id_col   = "Gene_ID",
    gene_type_col = "gene_type",
    fdr_level     = 0.05,
    top_k         = NULL,
    return_table = FALSE
) {
    for (pkg in c("SummarizedExperiment")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    validate_spe_contract(spe)

    rd <- SummarizedExperiment::rowData(spe)

    required_gene <- c("is_de", "true_logfc")
    miss_gene <- setdiff(required_gene, colnames(rd))
    if (length(miss_gene) > 0L) {
        stop("rowData(spe) is missing: ",
            paste(miss_gene, collapse = ", "),
            call. = FALSE)
    }

    if (!gene_id_col %in% colnames(rd)) {
        stop("rowData(spe) is missing gene id column '",
            gene_id_col, "'.", call. = FALSE)
    }

    if (!"lambda_g" %in% colnames(rd)) {
        stop("rowData(spe) is missing 'lambda_g' (needed for calibration strata).",
            call. = FALSE)
    }

    # -----------------------------
    # Extract DE random effects
    # Supports:
    #   - single iid DE term
    #   - mixture model (small + large)
    # -----------------------------

    if (!is.null(fit_inla$summary.random$Gene_ID_de_bin)) {
        re_single = fit_inla$summary.random$Gene_ID_de_bin
    } else {
        re_single = fit_inla$summary.random$Gene_ID_de
    }

    re_small  <- fit_inla$summary.random$Gene_ID_de_small
    re_large  <- fit_inla$summary.random$Gene_ID_de_large

    if (!is.null(re_single)) {

        # ---- OLD MODEL ----
        beta_mean <- re_single$mean
        beta_sd   <- re_single$sd
        beta_025_quant <- re_single$'0.025quant'
        beta_975_quant <- re_single$'0.975quant'

    } else if (!is.null(re_small) && !is.null(re_large)) {

        # ---- MIXTURE MODEL ----

        # Ensure alignment (same gene ordering)
        if (length(re_small$mean) != length(re_large$mean)) {
            stop("Mixture DE components have different lengths.", call. = FALSE)
        }

        beta_mean <- re_small$mean + re_large$mean

        # variance additivity approximation
        beta_sd <- sqrt(
            re_small$sd^2 +
            re_large$sd^2
        )

        # approximate quantiles (normal approx)
        beta_025_quant <- beta_mean - 1.96 * beta_sd
        beta_975_quant <- beta_mean + 1.96 * beta_sd

    } else {

        stop(
            "INLA fit does not contain DE random effects.\n",
            "Expected one of:\n",
            "  Gene_ID_de\n",
            "  Gene_ID_de_bin\n",
            "  Gene_ID_de_small + Gene_ID_de_large",
            call. = FALSE
        )
    }

    # -----------------------------
    # Align with rowData(spe)
    # -----------------------------
    #gid <- as.character(rd[[gene_id_col]])

    logfc_hat <- beta_mean
    logfc_sd  <- beta_sd
    logfc_025_quant <- beta_025_quant
    logfc_975_quant <- beta_975_quant

    logfc_hat[is.na(logfc_hat)] <- 0
    logfc_sd[is.na(logfc_sd)]   <- Inf
    logfc_025_quant[is.na(logfc_025_quant)] <- 0
    logfc_975_quant[is.na(logfc_975_quant)] <- 0

    # -----------------------------
    # Restrict to REG genes
    # -----------------------------
    gene_type <- if (gene_type_col %in% colnames(rd)) {
        as.character(rd[[gene_type_col]])
    } else {
        rep("REG", length(logfc_hat))
    }

    truth_de_all    <- as.numeric(rd$is_de)
    truth_logfc_all <- as.numeric(rd$true_logfc)

    logfc_all <- as.numeric(logfc_hat)
    lower_all <- as.numeric(logfc_025_quant)
    upper_all <- as.numeric(logfc_975_quant)

    gene_type_all <- gene_type
    lambda_all <- as.numeric(rd$lambda_g)

    keep <- gene_type == "REG"

    truth_de    <- as.numeric(rd$is_de)[keep]
    truth_logfc <- as.numeric(rd$true_logfc)[keep]

    logfc_eval <- as.numeric(logfc_hat[keep])
    sd_eval    <- as.numeric(logfc_sd[keep])

    lower_eval <- as.numeric(logfc_025_quant[keep])
    upper_eval <- as.numeric(logfc_975_quant[keep])

    # -----------------------------
    # Bayesian evidence score
    # -----------------------------
    z <- abs(logfc_eval) / sd_eval
    post_prob <- 2 * stats::pnorm(z) - 1   # P(|β| > 0 | data)

    # Bayesian analogue of p-value
    p_bayes <- 1 - post_prob

    # -----------------------------
    # Effect size metrics
    # -----------------------------
    mse  <- mean((logfc_eval - truth_logfc)^2)
    bias <- mean(logfc_eval - truth_logfc)
    spearman_r <- suppressWarnings(
        stats::cor(logfc_eval, truth_logfc, method = "spearman")
    )
    coverage <- mean(
        lower_eval <= truth_logfc &
        truth_logfc <= upper_eval
    )
    width = mean(
        upper_eval - lower_eval
    )

    # -----------------------------
    # Ranking metrics
    # -----------------------------
    auc <- NA_real_
    ap  <- NA_real_

    if (requireNamespace("pROC", quietly = TRUE)) {
        roc_obj <- pROC::roc(
        response = truth_de,
        predictor = post_prob,
        direction = "auto",
        quiet = TRUE
        )
        auc <- as.numeric(pROC::auc(roc_obj))
    }

    if (requireNamespace("PRROC", quietly = TRUE)) {
        fg <- post_prob[truth_de == 1]
        bg <- post_prob[truth_de == 0]
        ap <- PRROC::pr.curve(
        scores.class0 = fg,
        scores.class1 = bg,
        curve = FALSE
        )$auc.integral
    }

    # -----------------------------
    # Threshold-based metrics
    # -----------------------------
    ord <- order(p_bayes)
    ranked_truth <- truth_de[ord]

    # TPR at fixed FDR
    cum_fp <- cumsum(ranked_truth == 0)
    cum_tp <- cumsum(ranked_truth == 1)

    fdr <- cum_fp / pmax(cum_fp + cum_tp, 1)
    tpr <- cum_tp / sum(truth_de == 1)

    tpr_at_fdr <- NA_real_
    idx <- which(fdr <= fdr_level)
    if (length(idx) > 0L) {
        tpr_at_fdr <- max(tpr[idx])
    }

    # -----------------------------
    # Top-K metrics
    # -----------------------------
    if (is.null(top_k)) {
        #top_k <- sum(truth_de == 1)
        top_k <- max(10L, as.integer(ceiling(0.05 * length(truth_de))))
    }

    k <- min(top_k, length(truth_de))
    topk_truth <- ranked_truth[seq_len(k)]

    precision_k <- mean(topk_truth == 1)
    recall_k    <- sum(topk_truth == 1) / sum(truth_de == 1)

    diagnostics <- list(
        logfc_hat = logfc_hat,
        post_prob = post_prob,
        p_bayes = p_bayes, 
        n_genes = length(truth_de),
        n_de = sum(truth_de),
        min_fdr = min(fdr, na.rm = TRUE)
    )

    # ---------------------
    # Gene-wise metrics
    # ---------------------
    de_table <- NULL
    if (isTRUE(return_table)) {

        de_table <- data.frame(
            truth_de    = truth_de_all,
            truth_logfc = truth_logfc_all,
            logfc_hat   = logfc_all,
            logfc_sd    = logfc_sd, 
            lower       = lower_all,
            upper       = upper_all,
            lambda_g    = lambda_all,
            gene_type   = gene_type_all
        )
    }


    # -----------------------------
    # Output
    # -----------------------------
    list(
        # effect recovery
        mse = mse,
        bias = bias,
        spearman_r = spearman_r,
        coverage = coverage, 
        width = width, 

        # ranking
        auc = auc,
        average_precision = ap,

        # decision quality
        tpr_at_fdr = tpr_at_fdr,

        # top-k
        precision_k = precision_k,
        recall_k = recall_k,

        # raw outputs (for diagnostics)
        diagnostics = diagnostics,

        de_table = de_table 
    )
}

compute_dual_marker_performance <- function(
    results_df,
    positive_markers,
    negative_markers,
    method_name = "method",
    effect_col = NULL,
    prob_col = NULL,
    top_k = 100L
) {

    df <- results_df

    # ----------------------------
    # Detect effect column
    # ----------------------------
    if (is.null(effect_col)) {
        if ("logFC" %in% colnames(df)) {
            effect_col <- "logFC"
        } else if ("beta" %in% colnames(df)) {
            effect_col <- "beta"
        } else {
            stop("No effect column found.")
        }
    }

    # ----------------------------
    # Detect probability column
    # ----------------------------
    if (is.null(prob_col)) {
        if ("adj_p" %in% colnames(df)) {
            prob_col <- "adj_p"
        } else if ("post_prob" %in% colnames(df)) {
            prob_col <- "post_prob"
        }
    }

    df <- df[!is.na(df[[effect_col]]), ]

    n_genes <- nrow(df)

    # ----------------------------
    # Evidence score (scale-free)
    # ----------------------------
    if (!is.null(prob_col)) {
        if (prob_col == "adj_p") {
            df$evidence <- -log10(pmax(df[[prob_col]], .Machine$double.xmin))
        } else {
            df$evidence <- df[[prob_col]]
        }
    } else {
        df$evidence <- NA_real_
    }

    # Evidence percentile (higher = stronger)
    if (any(is.finite(df$evidence))) {
        r <- rank(-df$evidence, ties.method = "average")
        df$evidence_pct <- 1 - (r - 1) / (n_genes - 1)
    } else {
        df$evidence_pct <- NA_real_
    }

    # Top-K flag
    df$in_topK <- FALSE
    if (any(is.finite(df$evidence))) {
        ord <- order(df$evidence, decreasing = TRUE)
        k <- min(top_k, n_genes)
        df$in_topK[ord[seq_len(k)]] <- TRUE
    }

    # Rank by effect
    df$rank_effect <- rank(-df[[effect_col]], ties.method = "average")

    # Subsets
    pos_df <- df[df$gene %in% positive_markers, ]
    neg_df <- df[df$gene %in% negative_markers, ]

    # Coverage
    pos_cov <- nrow(pos_df) / length(positive_markers)
    neg_cov <- nrow(neg_df) / length(negative_markers)

    # ----------------------------
    # Direction Concordance
    # ----------------------------
    pos_cds <- mean(pos_df[[effect_col]] > 0)
    neg_cds <- mean(neg_df[[effect_col]] < 0)
    combined_cds <- mean(c(
        pos_df[[effect_col]] > 0,
        neg_df[[effect_col]] < 0
    ))

    # ----------------------------
    # AUC-like separation
    # ----------------------------
    pos_auc <- mean(pos_df$rank_effect <= n_genes / 2)
    neg_auc <- mean(neg_df$rank_effect > n_genes / 2)
    combined_auc <- mean(c(
        pos_df$rank_effect <= n_genes / 2,
        neg_df$rank_effect > n_genes / 2
    ))

    # ----------------------------
    # Relative Effect Enrichment
    # ----------------------------
    background_effect <- df[[effect_col]]

    pos_rel_eff <- median(abs(pos_df[[effect_col]]), na.rm = TRUE) /
        median(abs(background_effect), na.rm = TRUE)

    neg_rel_eff <- median(abs(neg_df[[effect_col]]), na.rm = TRUE) /
        median(abs(background_effect), na.rm = TRUE)

    # ----------------------------
    # Robust Standardized Effect
    # ----------------------------
    mad_val <- mad(background_effect, constant = 1, na.rm = TRUE) + 1e-12
    df$robust_z <- df[[effect_col]] / mad_val

    pos_robust <- median(abs(df$robust_z[df$gene %in% positive_markers]), na.rm = TRUE)
    neg_robust <- median(abs(df$robust_z[df$gene %in% negative_markers]), na.rm = TRUE)

    # ----------------------------
    # Evidence Percentile + TopK
    # ----------------------------
    pos_evidence_pct <- mean(df$evidence_pct[df$gene %in% positive_markers], na.rm = TRUE)
    neg_evidence_pct <- mean(df$evidence_pct[df$gene %in% negative_markers], na.rm = TRUE)

    pos_topK <- mean(df$in_topK[df$gene %in% positive_markers], na.rm = TRUE)
    neg_topK <- mean(df$in_topK[df$gene %in% negative_markers], na.rm = TRUE)

    # Composite biological score
    bio_score <- 0.5 * combined_cds + 0.5 * combined_auc

    data.frame(
        method = method_name,
        n_genes = n_genes,

        pos_coverage = pos_cov,
        neg_coverage = neg_cov,

        pos_CDS = pos_cds,
        neg_CDS = neg_cds,
        combined_CDS = combined_cds,

        pos_AUC = pos_auc,
        neg_AUC = neg_auc,
        combined_AUC = combined_auc,

        pos_rel_effect = pos_rel_eff,
        neg_rel_effect = neg_rel_eff,

        pos_robust_effect = pos_robust,
        neg_robust_effect = neg_robust,

        pos_mean_evidence_pct = pos_evidence_pct,
        neg_mean_evidence_pct = neg_evidence_pct,

        pos_topK = pos_topK,
        neg_topK = neg_topK,

        BioScore = bio_score
    )
}


#' Plot posterior for differential expression of a single gene
#'
#' Extracts the posterior marginal of the DE coefficient
#' f(Gene_ID_de, diff, model = "iid") and plots it.
#'
#' @param fit_obj Output of normalize_inla_pac() (list with $fit, $prep_full)
#' @param gene Gene identifier (name or ID)
#' @param gene_name_col Column in prep_full$df with gene names
#' @param gene_id_col Column in prep_full$df with integer Gene_ID
#' @param transform Scale for display: "logFC" or "FC"
#' @param prob Credible interval probabilities (length 3, e.g. c(.025,.5,.975))
#' @param rope Optional ROPE half-width on logFC scale (e.g. 0.1); NULL to skip
#'
#' @return List with ggplot object and a data.frame of posterior summaries
#' @export
plot_de_posterior <- function(
    fit_obj,
    gene,
    gene_name_col = "gene",
    gene_id_col   = "Gene_ID",
    transform = c("logFC", "FC"),
    prob = c(0.025, 0.5, 0.975),
    rope = NULL
) {
    transform <- match.arg(transform)

    if (is.null(fit_obj$prep_full) ||
        is.null(fit_obj$prep_full$df)) {
        stop("fit_obj$prep_full$df is required.", call. = FALSE)
    }

    df <- fit_obj$prep_full$df
    fit <- fit_obj$fit$fit

    # ------------------------------------------------------------
    # Resolve gene -> Gene_ID_de index
    # ------------------------------------------------------------
    if (is.character(gene)) {
        if (!gene_name_col %in% names(df)) {
            stop("gene_name_col not found in prep_full$df.", call. = FALSE)
        }
        rows <- which(df[[gene_name_col]] == gene)
        if (length(rows) == 0L) {
            stop("Gene not found: ", gene, call. = FALSE)
        }
        gene_id <- df[[gene_id_col]][rows[1L]]
    } else {
        gene_id <- as.integer(gene)
    }

    if (!"Gene_ID_de" %in% names(df)) {
        stop("prep_full$df must contain Gene_ID_de.", call. = FALSE)
    }

    gene_de_idx <- unique(df$Gene_ID_de[df[[gene_id_col]] == gene_id])
    gene_de_idx <- gene_de_idx[!is.na(gene_de_idx)]

    if (length(gene_de_idx) != 1L) {
        stop(
            "Expected exactly one Gene_ID_de for gene; got ",
            length(gene_de_idx),
            call. = FALSE
        )
    }

    # ------------------------------------------------------------
    # Extract posterior marginal for DE effect
    # ------------------------------------------------------------
    if (is.null(fit$marginals.random$Gene_ID_de)) {
        stop("DE marginals not found in INLA fit.", call. = FALSE)
    }

    marg <- fit$marginals.random$Gene_ID_de[[gene_de_idx]]

    # ------------------------------------------------------------
    # Posterior summaries on logFC scale
    # ------------------------------------------------------------
    qs_log <- INLA::inla.qmarginal(prob, marg)
    mean_log <- INLA::inla.emarginal(identity, marg)
    p_gt0 <- 1 - INLA::inla.pmarginal(0, marg)
    p_lt0 <- INLA::inla.pmarginal(0, marg)

    rope_prob <- NA_real_
    if (!is.null(rope)) {
        rope_prob <- INLA::inla.pmarginal(rope, marg) -
            INLA::inla.pmarginal(-rope, marg)
    }

    # ------------------------------------------------------------
    # Optional transform to FC scale
    # ------------------------------------------------------------
    if (transform == "FC") {
        marg_plot <- INLA::inla.tmarginal(exp, marg)
        qs_plot   <- exp(qs_log)
        xlab      <- "Fold change"
        vline     <- 1
    } else {
        marg_plot <- marg
        qs_plot   <- qs_log
        xlab      <- "log fold change"
        vline     <- 0
    }

    # ------------------------------------------------------------
    # Build data for plotting
    # ------------------------------------------------------------
    dens <- as.data.frame(marg_plot)
    colnames(dens) <- c("x", "density")

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------
    p <- ggplot2::ggplot(dens, ggplot2::aes(x = x, y = density)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_vline(
            xintercept = vline,
            linetype = "dashed",
            color = "grey40"
        ) +
        ggplot2::geom_vline(
            xintercept = qs_plot[c(1, 3)],
            linetype = "dotted",
            color = "grey60"
        ) +
        ggplot2::labs(
            title = paste0("DE posterior: ", gene),
            x = xlab,
            y = "Posterior density"
        ) +
        ggplot2::theme_minimal(base_size = 12)

    # ------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------
    summary <- data.frame(
        gene      = gene,
        scale     = transform,
        mean      = if (transform == "FC") exp(mean_log) else mean_log,
        q_low     = qs_plot[1],
        q_med     = qs_plot[2],
        q_high    = qs_plot[3],
        p_gt0     = p_gt0,
        p_lt0     = p_lt0,
        rope_prob = rope_prob,
        stringsAsFactors = FALSE
    )

    list(
        plot = p,
        summary = summary
    )
}


#' Evaluate spatial domain preservation after normalization (ARI)
#'
#' Given a SpatialExperiment with *true* domain labels (simulation truth),
#' this function evaluates how well a normalized assay preserves spatial
#' domains by running one or more clustering backends and computing ARI
#' between predicted clusters and truth.
#'
#' This is the domain-analogue of `evaluate_de_performance()`: it treats the
#' normalization as a black box (already applied) and evaluates a downstream
#' task (spatial clustering) against known truth.
#'
#' @param spe A SpatialExperiment containing domain truth in colData.
#' @param norm_assay Character. Name of assay to evaluate (default "logcounts").
#'   The clustering backends in this code use `logcounts(spe)` as input; if
#'   `norm_assay != "logcounts"`, the function will temporarily set
#'   `logcounts(spe)` to `assay(spe, norm_assay)` for evaluation.
#' @param domain_col Character. Column in colData(spe) holding true domains.
#'   Default "sim_domain". If not present and `fallback_annotated = TRUE`,
#'   uses "AnnotatedCluster" instead.
#' @param group_col Character. Optional grouping column. Not used for ARI,
#'   but validated if present in downstream codebases.
#' @param methods Character vector of clustering backends to run. Supported:
#'   "SNN", "BayesSpace", "SpaGCN".
#' @param snn_k Integer. k for SNN graph construction (default 10).
#' @param q Numeric multiplier for expected clusters in BayesSpace/SpaGCN
#'   (passed through to your helpers; default 1).
#' @param bayespace_nrep Numeric. Passed to BayesSpace::spatialCluster via
#'   your helper (if you later expose it); currently your helper hard-codes.
#' @param spagcn_res_init,spagcn_res_step Numeric. Passed to cluster_SpaGCN().
#' @param fallback_annotated Logical. If TRUE and domain_col missing, use
#'   colData(spe)$AnnotatedCluster.
#' @param drop_na_truth Logical. If TRUE, drops spots with NA truth before ARI.
#'   If FALSE, errors when truth contains NA (default TRUE).
#' @param return_clusters Logical. If TRUE, return predicted labels per method.
#' @param verbose Logical. If TRUE, emits warnings for method failures.
#'
#' @return A list with:
#'   - `ari`: named numeric vector with ARI per method (NA if method failed)
#'   - `diagnostics`: list with truth summary, n_spots_used, etc.
#'   - `clusters` (optional): list of (Truth, Predicted) per method
#' @export
evaluate_domain_preservation <- function(
  spe,
  norm_assay = "logcounts",
  domain_col = "sim_domain",
  group_col = "sim_group",
  methods = c("SNN", "BayesSpace", "SpaGCN"),
  snn_k = 10L,
  q = 1,
  bayespace_nrep = 5e4,
  spagcn_res_init = 0.7,
  spagcn_res_step = 0.1,
  fallback_annotated = TRUE,
  drop_na_truth = TRUE,
  return_clusters = FALSE,
  verbose = FALSE,
  strict = TRUE
) {
  # ---- base deps used regardless of method ----
  for (pkg in c("SummarizedExperiment", "SpatialExperiment", "SingleCellExperiment", "Matrix", "mclust")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Missing package: ", pkg, call. = FALSE)
    }
  }

  # ---- validate inputs ----
  if (!methods::is(spe, "SpatialExperiment")) {
    stop("`spe` must be a SpatialExperiment.", call. = FALSE)
  }
  if (!norm_assay %in% SummarizedExperiment::assayNames(spe)) {
    stop("Assay '", norm_assay, "' not found in `spe`.", call. = FALSE)
  }
  coords <- SpatialExperiment::spatialCoords(spe)
  if (is.null(coords) || ncol(coords) < 2L) {
    stop("`spe` must contain spatialCoords with >= 2 columns.", call. = FALSE)
  }

  cd <- SummarizedExperiment::colData(spe)

  truth_col <- domain_col
  if (!truth_col %in% colnames(cd)) {
    if (isTRUE(fallback_annotated) && "AnnotatedCluster" %in% colnames(cd)) {
      truth_col <- "AnnotatedCluster"
    } else {
      stop(
        "colData(spe) is missing domain truth column '", domain_col, "'.",
        if (fallback_annotated) " Also checked 'AnnotatedCluster'." else "",
        call. = FALSE
      )
    }
  }

  truth_raw <- cd[[truth_col]]
  if (drop_na_truth) {
    keep_spots <- !is.na(truth_raw)
    if (!any(keep_spots)) {
      stop("All truth domain labels are NA in colData(spe)$", truth_col, ".", call. = FALSE)
    }
  } else {
    if (any(is.na(truth_raw))) {
      stop("Truth domain labels contain NA; set drop_na_truth=TRUE to drop them.", call. = FALSE)
    }
    keep_spots <- rep(TRUE, length(truth_raw))
  }

  spe_eval <- spe[, keep_spots, drop = FALSE]
  cd_eval <- SummarizedExperiment::colData(spe_eval)
  truth <- as.factor(cd_eval[[truth_col]])
  if (nlevels(truth) < 2L) {
    stop("Need at least 2 unique domains in truth labels to compute ARI.", call. = FALSE)
  }

  # ---- dependency checks per method (FAIL FAST by default) ----
  methods <- unique(as.character(methods))
  supported <- c("SNN", "BayesSpace", "SpaGCN")
  bad <- setdiff(methods, supported)
  if (length(bad) > 0L) {
    stop("Unsupported methods: ", paste(bad, collapse = ", "),
         ". Supported: ", paste(supported, collapse = ", "), call. = FALSE)
  }

  deps <- list(
    SNN = c("scater", "scran", "igraph"),
    BayesSpace = c("BayesSpace"),
    SpaGCN = c("reticulate") # Python modules checked at runtime inside cluster_SpaGCN
  )

  missing <- lapply(methods, function(m) deps[[m]][!vapply(deps[[m]], requireNamespace, logical(1), quietly = TRUE)])
  names(missing) <- methods

  any_missing <- any(vapply(missing, length, integer(1)) > 0L)
  if (any_missing) {
    msg_vec <- vapply(
        names(missing),
        function(m) {
            if (length(missing[[m]]) == 0L) {
                ""
            } else {
                paste0(m, ": ", paste(missing[[m]], collapse = ", "))
            }
    },
    character(1)
    )
    msg <- paste(msg_vec[nzchar(msg_vec)], collapse = "\n")
    msg <- msg[nzchar(msg)]
    if (isTRUE(strict)) {
      stop("Missing required package(s) for requested clustering method(s):\n", msg, call. = FALSE)
    } else if (isTRUE(verbose)) {
      warning("Missing required package(s) for requested clustering method(s):\n", msg, call. = FALSE)
    }
  }

  # ---- standardize logcounts for downstream helpers ----
  y <- SummarizedExperiment::assay(spe_eval, norm_assay)
  if (any(!is.finite(y))) stop("`norm_assay` contains non-finite values.", call. = FALSE)
#   if (any(y < 0)) stop("`norm_assay` values must be non-negative (expected logcounts-like).", call. = FALSE)

  SingleCellExperiment::logcounts(spe_eval) <- y

  # computeSNN() expects spe$AnnotatedCluster
  cd2 <- SummarizedExperiment::colData(spe_eval)
  cd2$AnnotatedCluster <- truth
  SummarizedExperiment::colData(spe_eval) <- cd2

  # ---- PCA ----
  spe_eval <- addPCA(spe_eval)

  # ---- run methods ----
  ari <- stats::setNames(rep(NA_real_, length(methods)), methods)
  clusters_out <- vector("list", length(methods))
  names(clusters_out) <- methods

  run_one <- function(meth) {
    if (meth == "SNN") {
      g <- computeSNN(spe_eval, snn = as.integer(snn_k))
      cluster_graph(g)
    } else if (meth == "BayesSpace") {
      cluster_BayesSpace(spe_eval, q = q)
    } else if (meth == "SpaGCN") {
      cluster_SpaGCN(spe_eval, q = q, res_init = spagcn_res_init, res_step = spagcn_res_step)
    } else {
      stop("Internal error: unknown method ", meth, call. = FALSE)
    }
  }

  for (meth in methods) {
    # if strict=FALSE and package missing, skip and leave NA
    if (!isTRUE(strict) && length(missing[[meth]]) > 0L) next

    cl <- run_one(meth)  # let errors propagate (fail fast)

    cl$Truth <- as.factor(cl$Truth)
    cl$Predicted <- as.factor(cl$Predicted)

    if (length(cl$Truth) != length(cl$Predicted)) {
      stop("Method '", meth, "' returned mismatched Truth/Predicted lengths.", call. = FALSE)
    }

    ari[meth] <- rand_index(cl)
    if (isTRUE(return_clusters)) clusters_out[[meth]] <- cl
  }

  diagnostics <- list(
    truth_col = truth_col,
    n_spots_total = ncol(spe),
    n_spots_used = ncol(spe_eval),
    n_domains = nlevels(truth),
    domain_sizes = table(truth)
  )

  out <- list(ari = ari, diagnostics = diagnostics)
  if (isTRUE(return_clusters)) out$clusters <- clusters_out
  out
}
# Helper Functions 
# Note: Taken from SpaNorm\SpaNorm_files\SpaNorm_files\codes\spatial_clustering_wrappers.R
#----intermmediate analyses----
addPCA <- function(spe) {
  stopifnot(is(spe, 'SpatialExperiment'))

  #cluster
  spe = suppressWarnings(scater::runPCA(spe, ntop = nrow(spe), name = 'PCA', exprs_values = 'logcounts'))

  return(spe)
}

computeSNN <- function(spe, snn = 10) {
  g = suppressWarnings(scran::buildSNNGraph(spe, k = snn, use.dimred = 'PCA'))
  igraph::V(g)$Truth = as.character(spe$AnnotatedCluster)

  return(g)
}

#----cluster----
cluster_graph <- function(g, alg = igraph::cluster_walktrap, ...) {
  stopifnot(is(g, 'igraph'))

  cl = alg(g, ...)$membership
  cl = as.factor(cl)
  cl = list(
    'Truth' = as.factor(igraph::V(g)$Truth),
    'Predicted' = cl
  )

  return(cl)
}

cluster_BayesSpace <- function(spe, q = 1) {
  # expected clusters
  ncl = getNClusters(spe)
  ncl = round(ncl * q)

  # add BayesSpace metadata
  metadata(spe)$BayesSpace.data <- list()
  metadata(spe)$BayesSpace.data$platform <- "Visium"
  metadata(spe)$BayesSpace.data$is.enhanced <- FALSE
  rowData(spe)[["is.HVG"]] <- TRUE
  colnames(SpatialExperiment::spatialCoords(spe)) = c("col", "row")
  SummarizedExperiment::colData(spe) = cbind(SummarizedExperiment::colData(spe), SpatialExperiment::spatialCoords(spe))

  # identify clusters
  spe = BayesSpace::spatialCluster(spe, q = ncl, use.dimred = "PCA", nrep = 5e4)
  cl = list(
    "Truth" = spe$AnnotatedCluster,
    "Predicted" = as.factor(spe$spatial.cluster)
  )

  return(cl)
}


#' 
#' Error in py_call_impl(callable, call_args$unnamed, call_args$named) : 
#' AttributeError: 'csc_matrix' object has no attribute 'A'
#' Run `reticulate::py_last_error()` for details.
#' Calls: lapply ... run_one -> cluster_SpaGCN -> <Anonymous> -> py_call_impl
cluster_SpaGCN <- function(spe, q = 1, res_init = 0.7, res_step = 0.1) {
  # expected clusters
  ncl = getNClusters(spe)
  ncl = round(ncl * q)

  np     <- reticulate::import("numpy")
  pd     <- reticulate::import("pandas")
  spg    <- reticulate::import("SpaGCN")
  torch  <- reticulate::import("torch")
  random <- reticulate::import("random")

  # Convert SCE to AnnData
  assays(spe) = list('logcounts' = logcounts(spe))
  adata = zellkonverter::SCE2AnnData(spe)
  # ---- FIX: SpaGCN expects dense matrix ----
  if (inherits(adata$X, "scipy.sparse.csc_matrix") ||
      inherits(adata$X, "scipy.sparse.csr_matrix")) {
    adata$X <- adata$X$toarray()
  }

  # SpaGCN - build graph
  x_array = spatialCoords(spe)[, 1]
  y_array = spatialCoords(spe)[, 2]
  adj = spg$calculate_adj_matrix(x = x_array, y = y_array, histology = FALSE)

  # SpaGCN - Hyper-parameters
  p = 0.5
  l = spg$search_l(
    p,
    adj,
    start = 0.01,
    end = 1000,
    tol = 0.01,
    max_run = 100L
  )
  r_seed = t_seed = n_seed = 100L
  # search for suitable resolution
  res = spg$search_res(
    adata,
    adj,
    l,
    ncl,
    start = res_init,
    step = res_step,
    # start = 0.3,
    # step = 0.1,
    tol = 5e-3,
    lr = 0.05,
    max_epochs = 20L,
    r_seed = r_seed,
    t_seed = t_seed,
    n_seed = n_seed
  )

  # SpaGCN - run
  clf = spg$SpaGCN()
  clf$set_l(l)
  # set seed
  random$seed(r_seed)
  torch$manual_seed(t_seed)
  np$random$seed(n_seed)
  #Run
  clf$train(
    adata,
    adj,
    init_spa = TRUE,
    init = "louvain",
    res = res,
    tol = 5e-3,
    lr = 0.05,
    max_epochs = 200L
  )
  pred_res = clf$predict()
  adata$obs["pred"] = as.factor(pred_res[[1]])
  #Do cluster refinement(optional)
  #shape="hexagon" for Visium data, "square" for ST data.
  refined_pred = spg$refine(
    sample_id = rownames(adata$obs),
    pred = adata$obs[["pred"]],
    dis = adj,
    shape = "hexagon"
  )

  # Convert back to an SCE:
  cl = list(
    "Truth" = spe$AnnotatedCluster,
    "Predicted" = as.factor(refined_pred)
  )

  return(cl)
}

#----performance metrics----
rand_index <- function(x) {
  mclust::adjustedRandIndex(x$Predicted, x$Truth)
}

variance_analysis <- function(spe) {
  # filter very lowly expressed genes
  spe = spe[rowSums(logcounts(spe) != 0) > 5, ]
  x = t(as.matrix(SingleCellExperiment::logcounts(spe)))
  m1 = lm(x ~ spe$AnnotatedCluster)
  aov1 = car::Anova(m1, type = "II")
  df = summary(aov1, univariate = TRUE, p.adjust.method = "fdr")$univaov$SS
  df = t(df / df[, 1])[-1, ]
  colnames(df) = c("BetweenSS", "WithinSS")

  return(df)
}

getNClusters <- function(spe, domain_col = "sim_domain") {
  cd <- SummarizedExperiment::colData(spe)

  if (domain_col %in% colnames(cd)) {
    return(nlevels(as.factor(cd[[domain_col]])))
  }
  if ("AnnotatedCluster" %in% colnames(cd)) {
    return(nlevels(as.factor(cd$AnnotatedCluster)))
  }

  stop("No domain labels found to infer cluster count.", call. = FALSE)
}

#' 
#' @export
evaluate_svg_performance <- function(
  spe,
  norm_assay = "logcounts",
  truth_col = "is_svg", # TODO: Add effect size col. 
  k_frac = 0.10, # For recall metrics
  compute_iso_mse = TRUE, 
  methods = c("SPARKX", "nnSVG", "Seurat_Moran", "MERINGUE"),
  n_cores = 1L,
  verbose = FALSE,
  strict = FALSE
) {

  ## ---------------------------
  ## Basic validation
  ## ---------------------------
  stopifnot(methods::is(spe, "SpatialExperiment"))

  rd <- rowData(spe)
  if (!truth_col %in% colnames(rd)) {
    stop("Missing SVG truth column in rowData(spe): ", truth_col, call. = FALSE)
  }

  truth <- as.integer(rd[[truth_col]])
  ### Allow null-SVG simulations
  n_pos <- sum(truth == 1, na.rm = TRUE)
  n_neg <- sum(truth == 0, na.rm = TRUE)
  has_two_classes <- n_pos > 0L && n_neg > 0L

  gene_ids <- rownames(spe)
  obj_std <- .standardize_spot_ids(spe)
  logmat <- obj_std$mat
  coords <- obj_std$coords
  #logmat <- SummarizedExperiment::assay(spe, norm_assay)
  #coords <- SpatialExperiment::spatialCoords(spe)

  out_auc       <- setNames(rep(NA_real_, length(methods)), methods)
  out_prauc     <- out_auc
  out_recall_k  <- out_auc
  out_spearman  <- out_auc
  out_mse       <- out_auc
  out_iso_mse   <- out_auc

  raw_scores <- list()

  ## ---------------------------
  ## Method dispatch
  ## ---------------------------
  for (m in methods) {

    if (verbose) message("Running SVG method: ", m)

    res <- tryCatch({

      if (m == "SPARKX") {
        tbl <- callSVG.SPARKX(spe, n_cores = n_cores)
        score <- -log10(tbl$pvalue)

      } else if (m == "nnSVG") {
        tbl <- callSVG.nnSVG(
            spe, gene_ids = gene_ids, 
            n_threads = n_cores
        )
        score <- tbl$nnsvg_score

      } else if (m == "Seurat_Moran") {
        seu <- createSeuObject(logmat, normalize = FALSE)
        svf <- FindSpatiallyVariableFeatures(
          seu,
          selection.method = "moransi",
          spatial.location = coords,
          nfeatures = nrow(seu)
        )
        score <- svf$moranI

      } else if (m == "MERINGUE") {
        I <- callSVG.MERINGUE(logmat, coords)
        score <- I[, "I"]

      } else {
        stop("Unknown SVG method: ", m)
      }

      .standardize_svg_output(score, gene_ids)

    }, error = function(e) {
      if (strict) stop(e)
      if (verbose) message("❌ SVG method failed: ", m, " — ", e$message)
      NULL
    })

    if (is.null(res)) next

    ## ---------------------------
    ## Metrics
    ## ---------------------------
    df <- merge(
      data.frame(gene_id = gene_ids, truth = truth),
      res,
      by = "gene_id"
    )

    ## Ensure numeric
    df$score <- as.numeric(df$score)
    df$truth <- as.integer(df$truth)

    ## ---------- AUC ----------
    if (has_two_classes && requireNamespace("pROC", quietly = TRUE)) {
        out_auc[m] <- as.numeric(
            pROC::auc(df$truth, df$score)
        )
    }

    ## ---------- PR-AUC ----------
    if (has_two_classes && requireNamespace("PRROC", quietly = TRUE)) {
        fg <- df$score[df$truth == 1]
        bg <- df$score[df$truth == 0]
        if (length(fg) > 0L && length(bg) > 0L) {
            pr <- PRROC::pr.curve(
            scores.class0 = fg,
            scores.class1 = bg,
            curve = FALSE
            )
            out_prauc[m] <- pr$auc.integral
        }
    }

    ## ---------- Recall@K ----------
    if (n_pos > 0L) {
      K <- ceiling(k_frac * nrow(df))
      ord <- order(df$score, decreasing = TRUE)
      topK <- df[ord[seq_len(min(K, nrow(df)))], ]
      out_recall_k[m] <- sum(topK$truth == 1) / n_pos
    }

    ## ---------- Effect-size metrics ----------
    if ("true_svg_var" %in% colnames(rd)) {

      true_eff <- rd[df$gene_id, "true_svg_var"]
      ok <- is.finite(true_eff) & is.finite(df$score)

      if (sum(ok) > 1L) {
        out_spearman[m] <- suppressWarnings(
          cor(df$score[ok], true_eff[ok], method = "spearman")
        )
      }

      if (sum(ok) > 1L) {
        out_mse[m] <- mean(
          (scale(df$score[ok]) - scale(true_eff[ok]))^2,
          na.rm = TRUE
        )
      }

      if (compute_iso_mse && sum(ok) > 2L) {
        iso <- isoreg(df$score[ok], true_eff[ok])
        pred <- iso$yf[match(df$score[ok], iso$x)]
        out_iso_mse[m] <- mean((pred - true_eff[ok])^2, na.rm = TRUE)
      }
    }

    raw_scores[[m]] <- df

  }

  list(
    auc        = out_auc,
    pr_auc     = out_prauc,
    recall_k   = out_recall_k,
    spearman   = out_spearman,
    mse        = out_mse,
    iso_mse    = if (compute_iso_mse) out_iso_mse else NULL,
    raw        = raw_scores,
    diagnostics = list(
      n_genes = nrow(spe),
      n_svg   = sum(truth == 1),
      k_frac  = k_frac
    )
  )

}
.standardize_svg_output <- function(score, gene_ids) {
  stopifnot(length(score) == length(gene_ids))
  data.frame(
    gene_id = gene_ids,
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )
}

# From: SpaNorm\SpaNorm_files\SpaNorm_files\codes\SVG\SVG_wrappers.R
createGobject <- function(countMat, 
                          spatial_locs, 
                          cell_metadata, 
                          filter = FALSE, 
                          normalize = TRUE) {
  
  require(Giotto)
  
  python_path = NULL
  if(is.null(python_path)) {
    installGiottoEnvironment()
  }
  
  giotto = createGiottoObject(raw_exprs = countMat,
                                    spatial_locs = spatial_locs,
                                    cell_metadata = cell_metadata)
  
  if (filter == TRUE) {
    message("performing filtering...............")
    giotto <- filterGiotto(gobject = giotto,
                                 expression_threshold = 1,
                                 gene_det_in_min_cells = 50,
                                 min_det_genes_per_cell = 1000,
                                 expression_values = c('raw'),
                                 verbose = T)
    
  }
  if (normalize == TRUE) {
    message("performing normalization...............")
    giotto <- normalizeGiotto(gobject = giotto, scalefactor = 6000, verbose = FALSE)
    giotto <- addStatistics(gobject = giotto)
  }
  
  return(giotto)
}

createSPEObject <- function(countMat, 
                            spatial_locs, 
                            cell_metadata,
                            normalize = TRUE) {
  
  require(SpatialExperiment)
  require(S4Vectors)
  require(scran)
  
  rd = S4Vectors::DataFrame(symbol = rownames(countMat))
  
  data_spe = SpatialExperiment(
    assays = list(counts = countMat),
    colData = S4Vectors::DataFrame(cell_metadata),
    rowData = rd,
    spatialCoords = as.matrix(spatial_locs)
  )
  
  if (normalize == TRUE) {
    message("performing normalization...............")
    set.seed(123)
    qclus <- quickCluster(data_spe)
    data_spe <- computeSumFactors(data_spe, cluster = qclus)
    data_spe$sizeFactor <- pmax(1e-08,data_spe$sizeFactor)
    data_spe <- scater::logNormCounts(data_spe)
  }
  
  return(data_spe)
}

createSeuObject <- function(countMat, normalize=TRUE) {

  require(Seurat)

  seu = CreateSeuratObject(counts=countMat, assay="Spatial")

  if (normalize==TRUE) {
    message("performing normalization...............")
    seu = SCTransform(seu, assay = "Spatial", verbose = FALSE, return.only.var.genes=FALSE)
  }

  return(seu)
}


callSVG.Seurat <- function(logMat, method = "moransi", spatial.loc=spatial.loc) {
  
  require(Seurat)

  if(method == "moransi") {
    svf.info <- FindSpatiallyVariableFeatures(logMat, selection.method = "moransi",
                                          spatial.location = spatial.loc,
                                          nfeatures = nrow(seu))
    
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
    return(svf.info)
  }
}

callSVG.MERINGUE <- function(normMat, spatial_locs, verbose = FALSE) {
  require(MERINGUE)
  require(spdep)

  # ---- transpose: genes x spots -> spots x genes ----
  expr <- t(normMat)

  if (is.null(rownames(expr))) {
    stop("Expression matrix must have spot names as rownames.", call. = FALSE)
  }

  coords <- as.matrix(spatial_locs)
  rownames(coords) <- rownames(expr)

  stopifnot(identical(rownames(expr), rownames(coords)))

  if (verbose) {
    message("MERINGUE dims (spots x genes): ", paste(dim(expr), collapse = " x "))
  }

  # ---- build neighbor list ----
  nb <- spdep::knn2nb(
    spdep::knearneigh(coords, k = 6)
  )

  # ---- convert to spatial weights (THIS IS THE KEY FIX) ----
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  # ---- run MERINGUE ----
  I <- MERINGUE::getSpatialPatterns(expr, lw)

  I
}

.standardize_spot_ids <- function(spe, norm_assay = "logcounts") {

  mat <- SummarizedExperiment::assay(spe, norm_assay)
  coords <- SpatialExperiment::spatialCoords(spe)

  if (is.null(colnames(mat))) {
    stop("Expression matrix has NULL colnames.", call. = FALSE)
  }
  if (nrow(coords) != ncol(mat)) {
    stop("Mismatch: nrow(coords) != ncol(expression).", call. = FALSE)
  }

  spot_ids <- colnames(mat)
  rownames(coords) <- spot_ids

  list(
    mat = as.matrix(mat),
    coords = coords
  )
}

callSVG.Giotto <- function(gobject,
                           svg_method = c("kmeans", "rank"),
                           runHVG = FALSE,
                           runPCA = TRUE, n_cores = 1L) {
  
  if (runHVG == TRUE) {
    gobject <- Giotto::calculateHVG(gobject = gobject)
  }
  
  if (runPCA == TRUE) {
    gobject <- runPCA(gobject =  gobject, center = TRUE, scale_unit = TRUE)
  }
  
  gobject = createSpatialNetwork(gobject = gobject, minimum_k = 0)
  
  if (svg_method == "kmeans") {
    kmtest = Giotto::binSpect(gobject, do_parallel=TRUE, cores=n_cores)
    return(kmtest)
  } else {
    ranktest = Giotto::binSpect(gobject, bin_method = "rank", do_parallel=TRUE, cores=n_cores)
    return(ranktest)
  }
  
}

callSVG.SPARKX <- function(spe,
                          n_cores = 1L) {
  
  require(SPARK)
  require(SpatialExperiment)
  
  sp_count = counts(spe)
  spatial_loc = spatialCoords(spe)
  
  res <- sparkx(sp_count,spatial_loc, 
                numCores = n_cores,
                option="mixture")
  
  return(res$res_mtest)
  
}

# callSVG.nnSVG <- function(spe,
#                           seed = 1,
#                           return_spe = FALSE,
#                           n_threads = 1L, k=10) {


#   set.seed(seed)
#   spe <- nnSVG::nnSVG(spe, n_threads = n_threads, n_neighbors=k)
  
#   if (return_spe == FALSE) {
#     return(rowData(spe))
#   } else {
#     return(spe)
#   }
  
# }
callSVG.nnSVG <- function(spe, gene_ids, seed = 1, n_threads = 1L, k = 10) {

    require(nnSVG)

    set.seed(seed)

    spe_out <- nnSVG::nnSVG(spe)

    rd <- SummarizedExperiment::rowData(spe_out)

    if (!"prop_sv" %in% colnames(rd)) {
        stop(
        "nnSVG ran but `prop_sv` not found. Available columns: ",
        paste(colnames(rd), collapse = ", "),
        call. = FALSE
        )
    }

    data.frame(
        gene_id = rownames(rd),
        nnsvg_score = rd$prop_sv,
        stringsAsFactors = FALSE
    )

}
