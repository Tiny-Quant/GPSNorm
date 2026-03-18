# ============================================================
# AOI Cluster Consistency Table
# ============================================================

library(SpaNorm)
library(SpatialExperiment)
library(dplyr)
library(readr)

devtools::load_all(".")
library(GPSNorm)

# ------------------------------------------------------------
# Global configuration
# ------------------------------------------------------------

SIM_BASE_DIR  <- "simulations_4"
RESULTS_DIR   <- file.path(SIM_BASE_DIR, "results")
CACHE_DIR     <- file.path(SIM_BASE_DIR, "cache")

RESULTS_CSV   <- file.path(RESULTS_DIR, "aoi_cluster_results_top.csv")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR,   recursive = TRUE, showWarnings = FALSE)

base_seed <- 1000
trials    <- 25

# Top-% settings (computed on common gene universe)
TOP_PROP  <- 0.10
TOP_K_MIN <- 25

# ------------------------------------------------------------
# Load base data
# ------------------------------------------------------------

data(HumanDLPFC)
rownames(HumanDLPFC) <- rowData(HumanDLPFC)$gene_name

# ------------------------------------------------------------
# Helper metrics
# ------------------------------------------------------------

smape <- function(y_true, y_pred) {
    mean(
        abs(y_pred - y_true) /
            ((abs(y_true) + abs(y_pred)) / 2)
    ) * 100
}

relative_rmse <- function(y_true, y_pred) {
    sqrt(mean((y_pred - y_true)^2)) / sd(y_true)
}

jaccard_index <- function(a, b) {
    length(intersect(a, b)) / length(union(a, b))
}

topk_overlap <- function(rank_ref, rank_cmp, k) {
    length(intersect(rank_ref[1:k], rank_cmp[1:k])) / k
}

# Rank Biased Overlap (top-weighted, handles full rankings)
rbo_score <- function(rank1, rank2, p = 0.9) {
    k <- min(length(rank1), length(rank2))
    cumulative <- 0
    for (d in 1:k) {
        overlap_d <- length(intersect(rank1[1:d], rank2[1:d]))
        cumulative <- cumulative + (overlap_d / d) * p^(d - 1)
    }
    (1 - p) * cumulative
}

# Robust beta extraction from de_volcano_bayes()
get_de_beta <- function(spe, fit) {
    out <- de_volcano_bayes(spe, fit)$df
    if (!("beta" %in% names(out))) stop("de_volcano_bayes() df missing 'beta' column.")
    gene_key <- NULL
    if ("gene" %in% names(out)) {
        gene_key <- out$gene
    } else if (!is.null(rownames(out))) {
        gene_key <- rownames(out)
    } else {
        stop("de_volcano_bayes() df has no 'gene' column and no rownames; cannot align.")
    }
    setNames(out$beta, gene_key)
}

# Robust "true DE" gene IDs from sim$truth_gene
get_true_de_genes <- function(sim) {
    tg <- sim$truth_gene
    if (is.null(tg)) return(character(0))
    tg_df <- as.data.frame(tg)
    if (!("is_de" %in% names(tg_df))) return(character(0))

    if ("gene" %in% names(tg_df)) {
        tg_df |> filter(is_de == 1) |> pull(gene) |> unique()
    } else {
        # common pattern: gene IDs are rownames
        rn <- rownames(tg_df)
        if (is.null(rn)) return(character(0))
        rn[tg_df$is_de == 1] |> unique()
    }
}

# ------------------------------------------------------------
# Resume logic (skip trials already in CSV)
# ------------------------------------------------------------

completed_trials <- integer(0)
if (file.exists(RESULTS_CSV)) {
    tmp <- tryCatch(
        read_csv(RESULTS_CSV, show_col_types = FALSE),
        error = function(e) NULL
    )
    if (!is.null(tmp) && ("trial" %in% names(tmp))) {
        completed_trials <- sort(unique(tmp$trial))
    }
}

message("Completed trials detected: ",
        if (length(completed_trials) == 0) "<none>" else paste(completed_trials, collapse = ", "))

# ------------------------------------------------------------
# Main evaluation loop
# ------------------------------------------------------------

for (i in 1:trials) {

    if (i %in% completed_trials) {
        message("Skipping trial ", i, " (already complete)")
        next
    }

    message("===================================================")
    message("Running trial ", i)

    sim_file   <- file.path(CACHE_DIR, paste0("sim_", i, ".rds"))
    full_file  <- file.path(CACHE_DIR, paste0("full_fit_", i, ".rds"))
    half_file  <- file.path(CACHE_DIR, paste0("half_fit_", i, ".rds"))
    quart_file <- file.path(CACHE_DIR, paste0("quart_fit_", i, ".rds"))

    # --------------------------------------------------------
    # Simulation (cached)
    # --------------------------------------------------------
    if (file.exists(sim_file)) {
        message("Loading cached simulation: ", basename(sim_file))
        sim <- readRDS(sim_file)
    } else {
        message("Simulating data...")
        sim <- GPSNorm::simulate_spe_benchmark(
            spe = HumanDLPFC,
            domain_col = "AnnotatedCluster",
            n_slide = 3L,
            spot_prop = 0.1,
            gene_prop = 0.1,
            zi_prob = 0.00,
            tech_spatial_sd = 1.6,
            de_logfc_log_mean = 0.1,
            dom_sd = 1.5,
            svg_sd = 1.5,
            slide_lib_sd = 1.2,
            seed = base_seed + i
        )
        saveRDS(sim, sim_file)
    }

    # --------------------------------------------------------
    # Model fits (cached)
    # --------------------------------------------------------
    run_model <- function(file, ...) {
        if (file.exists(file)) {
            message("Loading cached fit: ", basename(file))
            readRDS(file)
        } else {
            message("Fitting model: ", basename(file))
            fit <- GPSNorm::normalize_spe_inla_pac(...)
            saveRDS(fit, file)
            fit
        }
    }

    full_fit <- run_model(
        full_file,
        sim$spe,
        model_fun = GPSNorm::fit_inla_baseline_model,
        model_args = list(verbose = TRUE)
    )

    half_fit <- run_model(
        half_file,
        sim$spe,
        model_fun = GPSNorm::fit_inla_baseline_model,
        model_args = list(verbose = TRUE),
        cluster_aoi = TRUE,
        cluster_args = list(k = 4, target_per_slide = 200L, verbose = TRUE)
    )

    quart_fit <- run_model(
        quart_file,
        sim$spe,
        model_fun = GPSNorm::fit_inla_baseline_model,
        model_args = list(verbose = TRUE),
        cluster_aoi = TRUE,
        cluster_args = list(k = 8, target_per_slide = 100L, verbose = TRUE)
    )

    # --------------------------------------------------------
    # OLD metrics: AOIs + Runtime
    # --------------------------------------------------------
    aoi_full  <- nrow(full_fit$prep_full$aoi)
    aoi_half  <- nrow(half_fit$prep_fit$aoi)
    aoi_quart <- nrow(quart_fit$prep_fit$aoi)

    runtime_full  <- as.numeric(full_fit$fit$fit$cpu.used["Total"])
    runtime_half  <- as.numeric(half_fit$fit$fit$cpu.used["Total"])
    runtime_quart <- as.numeric(quart_fit$fit$fit$cpu.used["Total"])

    # --------------------------------------------------------
    # DE extraction (REG filtering happens inside de_volcano_bayes)
    # --------------------------------------------------------
    full_de_raw  <- get_de_beta(full_fit$spe,  full_fit$fit$fit)
    half_de_raw  <- get_de_beta(half_fit$spe,  half_fit$fit$fit)
    quart_de_raw <- get_de_beta(quart_fit$spe, quart_fit$fit$fit)

    # --------------------------------------------------------
    # Gene overlap diagnostics (% mismatch computed PRE-alignment)
    # --------------------------------------------------------
    genes_full  <- names(full_de_raw)
    genes_half  <- names(half_de_raw)
    genes_quart <- names(quart_de_raw)

    common_genes <- Reduce(intersect, list(genes_full, genes_half, genes_quart))
    n_common <- length(common_genes)

    pct_mismatch_full  <- if (length(genes_full)  == 0) NA_real_ else 100 * (1 - n_common / length(genes_full))
    pct_mismatch_half  <- if (length(genes_half)  == 0) NA_real_ else 100 * (1 - n_common / length(genes_half))
    pct_mismatch_quart <- if (length(genes_quart) == 0) NA_real_ else 100 * (1 - n_common / length(genes_quart))

    message(sprintf(
        "Gene mismatch (%% not shared vs common): full=%.2f%% | half=%.2f%% | quart=%.2f%% (common=%d)",
        pct_mismatch_full, pct_mismatch_half, pct_mismatch_quart, n_common
    ))

    if (n_common == 0) {
        warning("No overlapping DE genes across resolutions for trial ", i, ". Writing row with NAs for DE metrics.")
    }

    # Align DE vectors to common genes (required for all DE comparison metrics)
    full_de  <- full_de_raw[common_genes]
    half_de  <- half_de_raw[common_genes]
    quart_de <- quart_de_raw[common_genes]

    # --------------------------------------------------------
    # OLD metrics: abs RMSE / rel RMSE / SMAPE / Spearman (DE)
    # --------------------------------------------------------
    rmse_abs_half  <- if (n_common == 0) NA_real_ else sqrt(mean((full_de - half_de)^2))
    rmse_abs_quart <- if (n_common == 0) NA_real_ else sqrt(mean((full_de - quart_de)^2))

    rmse_rel_half  <- if (n_common == 0) NA_real_ else relative_rmse(full_de, half_de)
    rmse_rel_quart <- if (n_common == 0) NA_real_ else relative_rmse(full_de, quart_de)

    smape_half     <- if (n_common == 0) NA_real_ else smape(full_de, half_de)
    smape_quart    <- if (n_common == 0) NA_real_ else smape(full_de, quart_de)

    spearman_half  <- if (n_common == 0) NA_real_ else cor(full_de, half_de, method = "spearman")
    spearman_quart <- if (n_common == 0) NA_real_ else cor(full_de, quart_de, method = "spearman")

    # --------------------------------------------------------
    # Top-% metrics (computed on common universe)
    # --------------------------------------------------------
    TOP_K <- if (n_common == 0) NA_integer_ else max(TOP_K_MIN, ceiling(TOP_PROP * n_common))

    true_de <- get_true_de_genes(sim)

    half_top <- data.frame(
        TopK_overlap  = NA_real_,
        Jaccard       = NA_real_,
        RBO           = NA_real_,
        TopK_Spearman = NA_real_,
        Recall_at_K   = NA_real_
    )

    quart_top <- half_top

    if (n_common > 0) {

        # Rankings by absolute effect size
        rank_full  <- names(sort(abs(full_de),  decreasing = TRUE))
        rank_half  <- names(sort(abs(half_de),  decreasing = TRUE))
        rank_quart <- names(sort(abs(quart_de), decreasing = TRUE))

        # Clamp TOP_K to length (should already be <= n_common)
        TOP_K_use <- min(TOP_K, n_common)

        top_ref   <- rank_full[1:TOP_K_use]
        top_half  <- rank_half[1:TOP_K_use]
        top_quart <- rank_quart[1:TOP_K_use]

        half_top <- data.frame(
            TopK_overlap  = topk_overlap(rank_full, rank_half, TOP_K_use),
            Jaccard       = jaccard_index(top_ref, top_half),
            RBO           = rbo_score(rank_full, rank_half),
            TopK_Spearman = cor(full_de[top_ref], half_de[top_ref], method = "spearman"),
            Recall_at_K   = if (length(true_de) == 0) NA_real_ else length(intersect(top_half, true_de)) / TOP_K_use
        )

        quart_top <- data.frame(
            TopK_overlap  = topk_overlap(rank_full, rank_quart, TOP_K_use),
            Jaccard       = jaccard_index(top_ref, top_quart),
            RBO           = rbo_score(rank_full, rank_quart),
            TopK_Spearman = cor(full_de[top_ref], quart_de[top_ref], method = "spearman"),
            Recall_at_K   = if (length(true_de) == 0) NA_real_ else length(intersect(top_quart, true_de)) / TOP_K_use
        )

        message("Top-% comparison uses TOP_K = ", TOP_K_use, " (TOP_PROP=", TOP_PROP, ", common=", n_common, ")")
    }

    # --------------------------------------------------------
    # OLD metrics: Spatial field recovery (as in original script)
    # --------------------------------------------------------
    # Center spatial effects
    est_full <- scale(unique(full_fit$fit$latent$spatial_eff), center = TRUE, scale = FALSE)

    est_half <- data.frame(
        spatial_eff = half_fit$fit$latent_full$spatial_eff,
        AOI_ID = half_fit$prep_full$df$AOI_ID
    ) |> distinct(AOI_ID, .keep_all = TRUE) |>
        pull(spatial_eff) |>
        scale(center = TRUE, scale = FALSE)

    est_quart <- data.frame(
        spatial_eff = quart_fit$fit$latent_full$spatial_eff,
        AOI_ID = quart_fit$prep_full$df$AOI_ID
    ) |> distinct(AOI_ID, .keep_all = TRUE) |>
        pull(spatial_eff) |>
        scale(center = TRUE, scale = FALSE)

    sp_recover <- sim$truth_spot |>
        as.data.frame() |>
        mutate(
            est_full  = as.numeric(est_full),
            est_half  = as.numeric(est_half),
            est_quart = as.numeric(est_quart)
        ) |>
        group_by(Slide_ID) |>
        mutate(
            logL_ctr = as.numeric(scale(log_L_all, center = TRUE, scale = FALSE)),
            correct_full  = logL_ctr + est_full,
            correct_half  = logL_ctr + est_half,
            correct_quart = logL_ctr + est_quart
        ) |>
        ungroup()

    # tech_exposure_ctr should exist in truth_spot (per your original script)
    sp_cor_full  <- cor(sp_recover$correct_full,  sp_recover$tech_exposure_ctr)
    sp_cor_half  <- cor(sp_recover$correct_half,  sp_recover$tech_exposure_ctr)
    sp_cor_quart <- cor(sp_recover$correct_quart, sp_recover$tech_exposure_ctr)

    sp_rmse_rel_full  <- relative_rmse(sp_recover$correct_full,  sp_recover$tech_exposure_ctr)
    sp_rmse_rel_half  <- relative_rmse(sp_recover$correct_half,  sp_recover$tech_exposure_ctr)
    sp_rmse_rel_quart <- relative_rmse(sp_recover$correct_quart, sp_recover$tech_exposure_ctr)

    sp_smape_full  <- smape(sp_recover$correct_full,  sp_recover$tech_exposure_ctr)
    sp_smape_half  <- smape(sp_recover$correct_half,  sp_recover$tech_exposure_ctr)
    sp_smape_quart <- smape(sp_recover$correct_quart, sp_recover$tech_exposure_ctr)

    # --------------------------------------------------------
    # Collect results (3 rows per trial: 1x, 0.5x, 0.25x)
    # --------------------------------------------------------
    result_table <- bind_rows(

        data.frame(
            Resolution = "1x",
            trial = i,
            AOIs = aoi_full,
            Runtime = runtime_full,

            RMSE_abs = NA_real_,
            RMSE_rel = NA_real_,
            SMAPE    = NA_real_,
            Spearman = NA_real_,

            Top_prop = NA_real_,
            Top_k = NA_integer_,
            TopK_overlap = NA_real_,
            Jaccard = NA_real_,
            RBO = NA_real_,
            TopK_Spearman = NA_real_,
            Recall_at_K = NA_real_,

            pct_mismatch = pct_mismatch_full,

            sp_cor = sp_cor_full,
            sp_RMSE_rel = sp_rmse_rel_full,
            sp_SMAPE = sp_smape_full
        ),

        data.frame(
            Resolution = "0.5x",
            trial = i,
            AOIs = aoi_half,
            Runtime = runtime_half,

            RMSE_abs = rmse_abs_half,
            RMSE_rel = rmse_rel_half,
            SMAPE    = smape_half,
            Spearman = spearman_half,

            Top_prop = TOP_PROP,
            Top_k = TOP_K,
            pct_mismatch = pct_mismatch_half,

            sp_cor = sp_cor_half,
            sp_RMSE_rel = sp_rmse_rel_half,
            sp_SMAPE = sp_smape_half
        ) |>
            bind_cols(half_top),

        data.frame(
            Resolution = "0.25x",
            trial = i,
            AOIs = aoi_quart,
            Runtime = runtime_quart,

            RMSE_abs = rmse_abs_quart,
            RMSE_rel = rmse_rel_quart,
            SMAPE    = smape_quart,
            Spearman = spearman_quart,

            Top_prop = TOP_PROP,
            Top_k = TOP_K,
            pct_mismatch = pct_mismatch_quart,

            sp_cor = sp_cor_quart,
            sp_RMSE_rel = sp_rmse_rel_quart,
            sp_SMAPE = sp_smape_quart
        ) |>
            bind_cols(quart_top)
    )

    message("Writing trial results...")

    write.table(
        result_table,
        file = RESULTS_CSV,
        sep = ",",
        row.names = FALSE,
        col.names = !file.exists(RESULTS_CSV),
        append = file.exists(RESULTS_CSV)
    )

    message("Trial ", i, " complete and saved.")
}