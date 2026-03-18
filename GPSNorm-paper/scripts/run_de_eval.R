# ============================================================
# DE Evaluation
# ============================================================

library(dplyr)
library(purrr)
library(SummarizedExperiment)

devtools::load_all(".")
library(GPSNorm)

# ------------------------------------------------------------
# Global configuration
# ------------------------------------------------------------

SIM_BASE_DIR <- "simulations_4"
NORM_DIR <- file.path(SIM_BASE_DIR, "normalized")
RESULTS_DIR <- file.path(SIM_BASE_DIR, "results")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

list_norm_files <- function(norm_dir) {
    list.files(
        norm_dir,
        pattern = "_rep[0-9]+\\.rds$",
        recursive = TRUE,
        full.names = TRUE
    )
}

# ------------------------------------------------------------
# Main evaluation loop
# ------------------------------------------------------------

evaluate_all_de <- function(
    norm_dir,
    fdr_level = 0.05,
    top_k = NULL,
    verbose = TRUE
) {

    files <- list_norm_files(norm_dir)

    if (length(files) == 0L) {
        stop("No normalized files found.", call. = FALSE)
    }

    if (verbose) {
        message("Found ", length(files), " normalized datasets.")
    }

    summary_list <- list()
    gene_list <- list()

    for (f in files) {

        if (verbose) {
            message("Evaluating DE: ", basename(f))
        }

        obj <- readRDS(f)

        spe <- obj$spe_norm
        method <- obj$method

        # -----------------------------
        # limma-based DE
        # -----------------------------
        out_limma <- evaluate_de_performance_limma(
            spe = spe,
            fdr_level = fdr_level,
            top_k = top_k,
            return_table = TRUE       
        )

        tbl <- tibble::tibble(
            scenario   = obj$scenario,
            replicate  = obj$replicate,
            method     = method,
            engine     = "limma",
            mse        = out_limma$mse,
            bias       = out_limma$bias %||% NA_real_,
            spearman_r = out_limma$spearman_r %||% NA_real_,
            coverage   = out_limma$coverage %||% NA_real_,
            width      = out_limma$width %||% NA_real_,
            auc        = out_limma$auc,
            pr_auc     = out_limma$ap %||% NA_real_,
            tpr_at_fdr = out_limma$tpr_at_fdr %||% NA_real_,
            precision_k = out_limma$precision_k %||% NA_real_,
            recall_k   = out_limma$recall_k %||% NA_real_
        )

        summary_list[[length(summary_list)+1]] <- tbl

        gene_tbl <- out_limma$de_table |>
            mutate(
                scenario  = obj$scenario,
                replicate = obj$replicate,
                method    = method,
                engine    = "limma"
            )

        gene_list[[length(gene_list)+1]] <- gene_tbl

        # -----------------------------
        # Optional Bayesian DE (INLA)
        # -----------------------------
        if (!is.null(obj$fit)) {

            out_bayes <- evaluate_de_performance_bayes(
                spe = spe,
                fit_inla = obj$fit,
                return_table = TRUE
            )

            tbl_bayes <- tibble::tibble(
                scenario   = obj$scenario,
                replicate  = obj$replicate,
                method     = method,
                engine     = "bayes",
                mse        = out_bayes$mse,
                bias       = out_bayes$bias %||% NA_real_,
                spearman_r = out_bayes$spearman_r %||% NA_real_,
                coverage   = out_bayes$coverage %||% NA_real_,
                width      = out_bayes$width %||% NA_real_,
                auc        = out_bayes$auc,
                pr_auc     = out_bayes$ap %||% NA_real_,
                tpr_at_fdr = out_bayes$tpr_at_fdr %||% NA_real_,
                precision_k = NA_real_,
                recall_k   = NA_real_
            )

            summary_list[[length(summary_list)+1]] <- tbl_bayes

            gene_tbl_b <- out_bayes$de_table |>
                mutate(
                    scenario  = obj$scenario,
                    replicate = obj$replicate,
                    method    = method,
                    engine    = "bayes"
                )

            gene_list[[length(gene_list)+1]] <- gene_tbl_b
        }
    }

    summary_tbl <- bind_rows(summary_list)
    gene_tbl <- bind_rows(gene_list)

    if (verbose) {
        message("Summary rows: ", nrow(summary_tbl))
        message("Gene-level rows: ", nrow(gene_tbl))
    }

    list(
        summary = summary_tbl,
        gene_level = gene_tbl
    )

}

# ------------------------------------------------------------
# Run evaluation
# ------------------------------------------------------------

de_results <- evaluate_all_de(
    norm_dir = NORM_DIR,
    fdr_level = 0.05,
    top_k = NULL,
    verbose = TRUE
)

saveRDS(
    de_results$summary,
    file.path(RESULTS_DIR, "de_results.rds")
)

write.csv(
    de_results$summary,
    file.path(RESULTS_DIR, "de_results.csv"),
    row.names = FALSE
)

saveRDS(
    de_results$gene_level,
    file.path(RESULTS_DIR, "gene_level_all.rds")
)