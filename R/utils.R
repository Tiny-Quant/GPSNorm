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
  na_level = "Unlabeled"
) {
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

#'
#' @export
de_volcano_bayes <- function(
    spe,
    fit_inla,
    gene_id_col = "Gene_ID",
    gene_type_col = "gene_type",
    delta = 0.1,
    label_top = 20,
    prob_cut = 0.95, 
    title = "Bayesian DE volcano (ROPE-based)",
    subtitle = NULL
) {
    rd <- SummarizedExperiment::rowData(spe)

    re <- fit_inla$summary.random$Gene_ID_de
    stopifnot(!is.null(re))

    beta <- re$mean
    sd   <- re$sd
    gid  <- as.character(re$ID)

    names(beta) <- gid
    names(sd)   <- gid

    beta_hat <- beta[as.character(rd[[gene_id_col]])]
    beta_sd  <- sd[as.character(rd[[gene_id_col]])]

    z <- (abs(beta_hat) - delta) / beta_sd
    post_prob <- pmax(0, 1 - 2 * pnorm(-z))

    df <- tibble::tibble(
        gene = rownames(spe),
        beta = beta_hat,
        post_prob = post_prob,
        gene_type = rd[[gene_type_col]]
    ) |>
        dplyr::filter(gene_type == "REG")

    df <- df |>
        dplyr::mutate(
            neglog_prob = -log10(pmax(1 - post_prob, 1e-12)),
            is_strong = post_prob >= prob_cut
        )

    # genes to label: top by posterior probability × effect size
    label_df <- df |>
        dplyr::filter(is_strong) |>
        dplyr::arrange(desc(post_prob), desc(abs(beta))) |>
        dplyr::slice_head(n = label_top)

    plot <- ggplot2::ggplot(df, ggplot2::aes(x = beta, y = neglog_prob)) +

        # main points
        ggplot2::geom_point(
        ggplot2::aes(color = is_strong),
            alpha = 0.7,
            size = 1.8
        ) +

        # ROPE vertical lines
        ggplot2::geom_vline(
            xintercept = c(-delta, delta),
            linetype = "dashed",
            linewidth = 0.6,
            color = "grey40"
        ) +

        # evidence threshold
        ggplot2::geom_hline(
            yintercept = -log10(1 - prob_cut),
            linetype = "dotted",
            linewidth = 0.6,
            color = "grey40"
        ) +

        # gene labels
        ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(label = gene),
            size = 3,
            max.overlaps = Inf,
            box.padding = 0.4,
            min.segment.length = 0
        ) +

        # scales
        ggplot2::scale_color_manual(
            values = c("FALSE" = "grey70", "TRUE" = "#0072B2"),
            labels = c("FALSE" = "Background", "TRUE" = "High evidence"),
            name = NULL
        ) +

        ggplot2::labs(
        title = title,
        subtitle = subtitle %||%
            sprintf("ROPE δ = %.2f, strong evidence ≥ %.2f", delta, prob_cut),
            x = "Posterior mean log2 fold-change",
            y = expression(-log[10](1 - P(abs(beta) > delta)))
        ) +

        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "top",
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold"),
            axis.title = ggplot2::element_text(face = "bold")
        )

    return(
        list(
            df = df, 
            plot = plot 
        )
    )
}


#' Volcano plot for limma DE results
#'
#'
#' @param spe A SpatialExperiment (or SummarizedExperiment) with `logcounts`.
#' @param norm_assay Assay name used for limma (default "logcounts").
#' @param group_col Column in colData(spe) with binary 0/1 group labels.
#' @param gene_type_col Column in rowData(spe) with gene type labels; REG used.
#' @param fdr_level FDR threshold for calling significant genes (default 0.05).
#' @param use_adj_p If TRUE, volcano y-axis uses -log10(adj.P.Val); else P.Value.
#' @param label_top Number of top genes to label (default 20).
#' @param label_by How to choose labeled genes: "fdr_then_abslogfc" or "signal".
#' @param title Plot title.
#' @param subtitle Optional subtitle.
#'
#' @return list(df = label_df, plot = ggplot_object, all = full_df)
#' @export
de_volcano_limma <- function(
    spe,
    norm_assay = "logcounts",
    group_col = "sim_group",
    gene_type_col = "gene_type",
    fdr_level = 0.05,
    use_adj_p = TRUE,
    label_top = 20,
    label_by = c("fdr_then_abslogfc", "signal"),
    title = "limma DE volcano",
    subtitle = NULL
) {
    label_by <- match.arg(label_by)

    for (pkg in c("SummarizedExperiment", "limma", "ggplot2",
                  "dplyr", "tibble", "ggrepel")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Missing package: ", pkg, call. = FALSE)
        }
    }

    validate_spe_contract(spe, assay = norm_assay)

    rd <- SummarizedExperiment::rowData(spe)
    cd <- SummarizedExperiment::colData(spe)

    if (!group_col %in% colnames(cd)) {
        stop("Missing group column: ", group_col, call. = FALSE)
    }

    y <- SummarizedExperiment::assay(spe, norm_assay)
    if (!is.matrix(y)) y <- as.matrix(y)

    if (any(!is.finite(y))) {
        stop("Assay contains non-finite values.", call. = FALSE)
    }

    group_vec <- as.factor(cd[[group_col]])

    if (!all(levels(group_vec) %in% c("0", "1"))) {
        stop("`group_col` must be binary with levels '0' and '1'.", call. = FALSE)
    }

    if (anyNA(group_vec)) {
        keep <- !is.na(group_vec)
        warning("Dropping ", sum(!keep), " samples with NA group.")
        spe <- spe[, keep, drop = FALSE]
        y   <- SummarizedExperiment::assay(spe, norm_assay)
        cd  <- SummarizedExperiment::colData(spe)
        group_vec <- as.factor(cd[[group_col]])
    }

    design <- stats::model.matrix(~ group_vec)

    fit <- limma::lmFit(y, design)
    fit <- limma::eBayes(fit)

    # --- SAFETY CHECKS ---

    n_genes_input <- nrow(y)
    n_genes_fit   <- nrow(fit$coefficients)

    if (n_genes_input != n_genes_fit) {
        stop("Gene count mismatch between input and limma fit.",
             call. = FALSE)
    }

    if (!identical(rownames(y), rownames(fit$coefficients))) {
        stop("Row order mismatch between assay and limma fit.",
             call. = FALSE)
    }

    na_coef <- sum(is.na(fit$coefficients[, 2]))
    if (na_coef > 0) {
        warning(na_coef, " genes have NA coefficient estimates.")
    }

    # --- Construct full results table (no gene dropping) ---

    df <- tibble::tibble(
        gene     = rownames(y),
        logFC    = fit$coefficients[, 2],
        p_value  = fit$p.value[, 2],
        adj_p    = p.adjust(fit$p.value[, 2], method = "BH")
    )

    # Attach gene_type if available
    if (!is.null(gene_type_col) && gene_type_col %in% colnames(rd)) {
        df$gene_type <- as.character(rd[[gene_type_col]])
    } else {
        df$gene_type <- "REG"
    }

    # Volcano calculations
    y_p <- if (use_adj_p) df$adj_p else df$p_value
    y_p <- pmax(y_p, .Machine$double.xmin)

    df <- dplyr::mutate(
        df,
        neglogp = -log10(y_p),
        is_sig  = .data$adj_p <= fdr_level
    )

    # --- Label selection ---

    label_df <- df

    if (label_by == "fdr_then_abslogfc") {
        label_df <- label_df |>
            dplyr::filter(.data$is_sig) |>
            dplyr::arrange(.data$adj_p,
                           dplyr::desc(abs(.data$logFC)))
    } else {
        label_df <- label_df |>
            dplyr::mutate(signal = .data$neglogp * abs(.data$logFC)) |>
            dplyr::arrange(dplyr::desc(.data$signal))
    }

    label_df <- dplyr::slice_head(label_df, n = label_top)

    # --- Plot ---

    hline <- -log10(fdr_level)
    ylab_txt <- if (use_adj_p) "-log10(FDR)" else "-log10(p-value)"

    plot <- ggplot2::ggplot(df,
        ggplot2::aes(x = logFC, y = neglogp)) +
        ggplot2::geom_point(
            ggplot2::aes(color = is_sig),
            alpha = 0.7,
            size = 1.8
        ) +
        ggplot2::geom_hline(
            yintercept = hline,
            linetype = "dotted",
            linewidth = 0.6,
            color = "grey40"
        ) +
        ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(label = gene),
            size = 3,
            max.overlaps = Inf,
            box.padding = 0.4,
            min.segment.length = 0
        ) +
        ggplot2::scale_color_manual(
            values = c("FALSE" = "grey70",
                       "TRUE" = "#D55E00"),
            name = NULL
        ) +
        ggplot2::labs(
            title = title,
            subtitle = subtitle,
            x = "limma log2 fold-change",
            y = ylab_txt
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            legend.position = "top",
            plot.title = ggplot2::element_text(face = "bold"),
            axis.title = ggplot2::element_text(face = "bold")
        )

    list(
        df   = label_df,
        plot = plot,
        all  = df
    )
}