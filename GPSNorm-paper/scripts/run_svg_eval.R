# ============================================================
# SVG Evaluation
# ============================================================

library(SpatialExperiment)
library(dplyr)
library(reticulate)
library(GPSNorm)

# ------------------------------------------------------------
# Thread safety
# ------------------------------------------------------------
Sys.setenv(
    OMP_NUM_THREADS = 1,
    OPENBLAS_NUM_THREADS = 1,
    MKL_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1,
    NUMEXPR_NUM_THREADS = 1
)

# ------------------------------------------------------------
# Activate isolated SVG environment
# ------------------------------------------------------------
use_condaenv("gpsnorm-svg", required = TRUE)

# ------------------------------------------------------------
# Global configuration
# ------------------------------------------------------------
SIM_BASE_DIR <- "simulations_4"
NORM_DIR    <- file.path(SIM_BASE_DIR, "normalized")
RESULTS_DIR <- file.path(SIM_BASE_DIR, "results")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

RESULTS_CSV <- file.path(RESULTS_DIR, "svg_results_checkpoint.csv")

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
# Sequential SVG evaluation with CSV checkpointing
# ------------------------------------------------------------
evaluate_all_svg_sequential_csv <- function(
    norm_dir,
    results_csv,
    norm_assay = "logcounts",
    truth_col = "is_svg",
    k_frac = 0.10,
    methods = c("SPARKX", "nnSVG", "Seurat_Moran", "MERINGUE"),
    compute_iso_mse = TRUE,
    n_cores = 1L,
    strict = FALSE,
    verbose = TRUE
) {

    files <- list_norm_files(norm_dir)

    if (length(files) == 0L) {
        stop("No normalized files found.", call. = FALSE)
    }

    if (verbose) {
        message("Found ", length(files), " normalized datasets.")
    }

    # ----------------------------------------------------------
    # Load existing results (restart safeguard)
    # ----------------------------------------------------------
    if (file.exists(results_csv)) {
        results_existing <- read.csv(results_csv, stringsAsFactors = FALSE)

        completed_keys <- with(
            results_existing,
            paste(scenario, replicate, method, sep = "::")
        )

        if (verbose) {
            message("Loaded existing CSV with ", nrow(results_existing), " records.")
        }
    } else {
        completed_keys <- character()

        header <- tibble(
            scenario   = character(),
            replicate  = integer(),
            method     = character(),
            auc        = numeric(),
            pr_auc     = numeric(),
            recall_k   = numeric(),
            spearman   = numeric(),
            mse        = numeric(),
            iso_mse    = numeric()
        )

        write.table(
            header,
            file = results_csv,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE
        )

        if (verbose) {
            message("Created results CSV: ", results_csv)
        }
    }

    # ----------------------------------------------------------
    # Main loop
    # ----------------------------------------------------------
    for (f in files) {

        obj <- readRDS(f)
        key <- paste(obj$scenario, obj$replicate, obj$method, sep = "::")

        if (key %in% completed_keys) {
            if (verbose) {
                message("⏭️  Skipping completed: ", basename(f))
            }
            next
        }

        if (verbose) {
            message("▶️  Starting SVG evaluation: ", basename(f))
        }

        spe <- obj$spe_norm

        out <- tryCatch(
            GPSNorm::evaluate_svg_performance(
                spe = spe,
                norm_assay = norm_assay,
                truth_col = truth_col,
                k_frac = k_frac,
                compute_iso_mse = compute_iso_mse,
                methods = methods,
                n_cores = n_cores,
                verbose = FALSE,
                strict = strict
            ),
            error = function(e) {
                message(
                    "❌ SVG evaluation failed for ",
                    basename(f), ": ",
                    e$message
                )
                NULL
            }
        )

        if (is.null(out)) {
            next
        }

        row <- tibble(
            scenario   = obj$scenario,
            replicate  = obj$replicate,
            method     = obj$method,
            auc        = out$auc,
            pr_auc     = out$pr_auc,
            recall_k   = out$recall_k,
            spearman   = out$spearman,
            mse        = out$mse,
            iso_mse    = out$iso_mse
        )

        write.table(
            row,
            file = results_csv,
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE
        )

        completed_keys <- c(completed_keys, key)

        if (verbose) {
            message("✅ Completed and written: ", basename(f))
        }
    }

    invisible(NULL)
}

    # ------------------------------------------------------------
    # Run evaluation
    # ------------------------------------------------------------
evaluate_all_svg_sequential_csv(
    norm_dir = NORM_DIR,
    results_csv = RESULTS_CSV,
    norm_assay = "logcounts",
    truth_col = "is_svg",
    k_frac = 0.10,
    methods = "nnSVG",
    compute_iso_mse = TRUE,
    n_cores = 1L,
    strict = FALSE,
    verbose = TRUE
)