# ============================================================
# ARI / Domain Preservation Evaluation
# ============================================================

library(SpatialExperiment)
library(dplyr)
library(reticulate)

devtools::load_all(".")
library(GPSNorm)

# ------------------------------------------------------------
# Python environment (required for some methods)
# ------------------------------------------------------------

use_condaenv("gpsnorm-ari", required = TRUE)

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

evaluate_all_ari <- function(
    norm_dir,
    norm_assay = "logcounts",
    methods = "SNN",
    strict = TRUE,
    verbose = TRUE
) {

    files <- list_norm_files(norm_dir)

    if (length(files) == 0L) {
        stop("No normalized files found.", call. = FALSE)
    }

    if (verbose) {
        message("Found ", length(files), " normalized datasets.")
    }

    res <- lapply(files, function(f) {

        if (verbose) {
            message("Evaluating ARI: ", basename(f))
        }

        obj <- readRDS(f)

        spe <- obj$spe_norm

        out <- evaluate_domain_preservation(
            spe,
            norm_assay = norm_assay,
            methods = methods,
            strict = strict
        )

        ari_tbl <- as.data.frame(t(out$ari))

        ari_tbl$scenario  <- obj$scenario
        ari_tbl$replicate <- obj$replicate
        ari_tbl$method    <- obj$method

        ari_tbl
    })

    bind_rows(res)
}

# ------------------------------------------------------------
# Run evaluation
# ------------------------------------------------------------

ari_results <- evaluate_all_ari(
    norm_dir = NORM_DIR,
    norm_assay = "logcounts",
    methods = "SNN",
    strict = TRUE,
    verbose = TRUE
)

saveRDS(
    ari_results,
    file.path(RESULTS_DIR, "ari_results.rds")
)

write.csv(
    ari_results,
    file.path(RESULTS_DIR, "ari_results.csv"),
    row.names = FALSE
)

