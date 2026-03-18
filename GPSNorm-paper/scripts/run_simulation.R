# ============================================================
# Simulation + Normalization
# ============================================================

library(dplyr)
library(purrr)
library(SpaNorm)
library(SpatialExperiment)
library(SummarizedExperiment)

devtools::load_all(".")
library(GPSNorm)

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

sim_paths <- function(base_dir) {
    list(
        base       = base_dir,
        normalized = file.path(base_dir, "normalized")
    )
}

# ------------------------------------------------------------
# File / replicate helpers
# ------------------------------------------------------------

rep_file_path <- function(base_dir, method, scenario_id, rep_id) {
    file.path(
        base_dir,
        method,
        sprintf("%s_rep%d.rds", scenario_id, rep_id)
    )
}

rep_exists <- function(base_dir, method, scenario_id, rep_id) {
    file.exists(
        rep_file_path(base_dir, method, scenario_id, rep_id)
    )
}

# ------------------------------------------------------------
# Find first incomplete replicate (global barrier)
# ------------------------------------------------------------

first_incomplete_rep <- function(
    scenarios,
    methods,
    base_dir,
    verbose = TRUE
) {

    r <- 1L

    repeat {

        missing <- list()

        for (sid in scenarios$scenario_id) {
            for (m in methods) {

                if (!rep_exists(base_dir, m, sid, r)) {
                    missing[[length(missing) + 1]] <-
                        paste(sid, m, sep = "::")
                }
            }
        }

        if (length(missing) > 0L) {

            if (verbose) {
                message(
                    "STATUS: replicate ", r,
                    " is INCOMPLETE (",
                    length(missing),
                    " missing outputs)"
                )
                message(
                    "  Missing examples: ",
                    paste(head(missing, 5L), collapse = ", "),
                    if (length(missing) > 5L) " ..."
                )
            }

            return(r)
        }

        if (verbose) {
            message("STATUS: replicate ", r, " is complete")
        }

        r <- r + 1L
    }
}


# ------------------------------------------------------------
# Global configuration
# ------------------------------------------------------------

SIM_BASE_DIR <- "simulations_4"
paths <- sim_paths(SIM_BASE_DIR)

EXPECTED_METHODS <- c(
    "none",
    "q3",
    "spanorm",
    "scran",
    "sct",
    "giotto",
    "ls",
    "inla",
    #"inla_mix"
    "inla_tuned",
    "inla_heavy"
    #"inla_bin"
)

inla_methods <- list(

    inla = list(
        model_fun  = GPSNorm::fit_inla_baseline_model,
        model_args = list(k = 4, verbose = FALSE)
    ),

    inla_mix = list(
        model_fun  = GPSNorm::fit_inla_mix,
        model_args = list(k = 4, verbose = TRUE)
    ),

    inla_tuned = list(
        model_fun  = GPSNorm::fit_inla_tuned_pc,
        model_args = list(k = 4, verbose = TRUE)
    ),

    inla_heavy = list(
        model_fun  = GPSNorm::fit_inla_heavy_tail,
        model_args = list(k = 4, verbose = TRUE)
    )

    # inla_bin = list(
    #     model_fun  = GPSNorm::fit_inla_bin,
    #     model_args = list(k = 4, verbose = TRUE)
    # )

)

# ------------------------------------------------------------
# Load reference dataset
# ------------------------------------------------------------

data(HumanDLPFC, package = "SpaNorm")
rownames(HumanDLPFC) <- rowData(HumanDLPFC)$gene_name

# ------------------------------------------------------------
# Simulation scenarios
# ------------------------------------------------------------

scenarios <- tibble::tribble(
    ~scenario_id, ~tech_spatial_sd, ~slide_lib_sd, ~dom_sd,
    ~svg_sd, ~de_svg_overlap, ~de_dom_overlap,
    ~group_confounded, ~group_conf_rho, ~group_conf_mode,

    "A_overlap0",        0.8, 0.6, 1.0, 1.0, 0.0, 0.0,
                         FALSE, 0.0, "logit",

    "B_overlap30",       0.8, 0.6, 1.0, 1.0, 0.3, 0.3,
                         FALSE, 0.0, "logit",

    "C_overlap70",       0.8, 0.6, 1.0, 1.0, 0.7, 0.7,
                         FALSE, 0.0, "logit",

    "D_slideTech",       0.6, 1.2, 1.0, 1.0, 0.3, 0.3,
                         FALSE, 0.0, "logit",

    "E_aoiSlideTech",    1.6, 1.2, 1.0, 1.0, 0.3, 0.3,
                         FALSE, 0.0, "logit",

    "F_hard",            1.6, 1.2, 1.5, 1.5, 0.7, 0.7,
                         FALSE, 0.0, "logit",

    # --------------------------------------------------------
    # Scenario G: spatially confounded group
    # --------------------------------------------------------

    "G_conf_strong",     1.2, 1.0, 1.0, 1.0, 0.3, 0.3,
                         TRUE, 2.0, "logit"
)

# ------------------------------------------------------------
# Simulation wrapper
# ------------------------------------------------------------

run_one_sim <- function(scen, seed) {

    GPSNorm::simulate_spe_benchmark(
        spe = HumanDLPFC,
        domain_col = "AnnotatedCluster",
        n_slide = 3L,
        spot_prop = 0.05,
        gene_prop = 0.05,
        zi_prob = 0.00,
        tech_spatial_sd = scen$tech_spatial_sd,
        slide_lib_sd    = scen$slide_lib_sd,
        de_logfc_log_mean = 0.1,
        dom_sd = scen$dom_sd,
        svg_sd = scen$svg_sd,
        de_svg_overlap = scen$de_svg_overlap,
        de_dom_overlap = scen$de_dom_overlap,
        seed = seed
    )
}

# ------------------------------------------------------------
# Normalization + saving (skip-safe)
# ------------------------------------------------------------

normalize_and_save <- function(
    sim,
    scenario_id,
    rep_id,
    seed,
    paths,
    verbose = TRUE
) {

    spe_raw <- sim$spe

    save_one <- function(method, spe_norm, fit = NULL) {

        out_dir <- file.path(paths$normalized, method)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        out_path <- rep_file_path(
            paths$normalized,
            method,
            scenario_id,
            rep_id
        )

        if (file.exists(out_path)) {
            if (verbose) {
                message(
                    "  SKIP  ",
                    scenario_id,
                    " rep ",
                    rep_id,
                    " [",
                    method,
                    "] already exists"
                )
            }
            return(invisible(FALSE))
        }

        obj <- list(
            scenario   = scenario_id,
            replicate  = rep_id,
            seed       = seed,
            method     = method,
            spe_norm   = spe_norm,
            truth_gene = sim$truth_gene,
            truth_spot = sim$truth_spot,
            fit        = fit
        )

        saveRDS(obj, out_path)

        if (verbose) {
            message(
                "  WRITE ",
                scenario_id,
                " rep ",
                rep_id,
                " [",
                method,
                "]"
            )
        }

        invisible(TRUE)
    }

    spe_none <- spe_raw
    logcounts(spe_none) <- log2(counts(spe_raw) + 1)
    save_one("none", spe_none)

    save_one("q3", apply_q3_logcounts(spe_none))
    save_one("spanorm", SpaNorm::SpaNorm(spe_raw, sample.p = 0.5))
    save_one("scran", normalise_scran(spe_raw))
    save_one("sct", normalise_sct(spe_raw))
    save_one("giotto", normalise_Giotto(spe_raw))
    save_one("ls", normalise_ls(spe_raw))

    should_run <- function(method) {

        out_path <- rep_file_path(
            paths$normalized,
            method,
            scenario_id,
            rep_id
        )

        !file.exists(out_path)
    }

    for (method_name in names(inla_methods)) {

        if (!should_run(method_name)) {

            if (verbose) {
                message(
                    "  SKIP  ",
                    scenario_id,
                    " rep ",
                    rep_id,
                    " [",
                    method_name,
                    "] exists"
                )
            }

            next
        }

        if (verbose) {
            message(
                "  RUN   ",
                scenario_id,
                " rep ",
                rep_id,
                " [",
                method_name,
                "]"
            )
        }

        spec <- inla_methods[[method_name]]

        spe_inla <- GPSNorm::normalize_spe_inla_pac(
            spe_raw,
            model_fun  = spec$model_fun,
            model_args = spec$model_args
        )

        save_one(
            method_name,
            spe_inla$spe,
            fit = spe_inla$fit$fit
        )
    }

    invisible(TRUE)
}

# ------------------------------------------------------------
# Main execution loop (transactional catch-up)
# ------------------------------------------------------------

resume_simulation <- function(
    scenarios,
    paths,
    n_new_reps = 1L,
    base_seed = 1000,
    verbose = TRUE
) {

    for (iter in seq_len(n_new_reps)) {

        r <- first_incomplete_rep(
            scenarios = scenarios,
            methods   = EXPECTED_METHODS,
            base_dir  = paths$normalized,
            verbose   = verbose
        )

        message(
            "\n=== STARTING REPLICATE ",
            r,
            " ==="
        )

        for (i in seq_len(nrow(scenarios))) {

            sid <- scenarios$scenario_id[i]
            seed <- base_seed * i + r

            missing <- FALSE

            for (m in EXPECTED_METHODS) {
                if (!rep_exists(
                        paths$normalized,
                        m,
                        sid,
                        r
                    )) {
                    missing <- TRUE
                    break
                }
            }

            if (!missing) {
                message(
                    "SKIP scenario ",
                    sid,
                    " rep ",
                    r,
                    " (already complete)"
                )
                next
            }

            message(
                "SIMULATE ",
                sid,
                " rep ",
                r,
                " (seed=",
                seed,
                ")"
            )

            sim <- run_one_sim(scenarios[i, ], seed)

            normalize_and_save(
                sim,
                scenario_id = sid,
                rep_id = r,
                seed = seed,
                paths = paths,
                verbose = verbose
            )
        }

        message(
            "=== FINISHED REPLICATE ",
            r,
            " ===\n"
        )
    }

    invisible(TRUE)
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

resume_simulation(
    scenarios = scenarios,
    paths = paths,
    n_new_reps = 10,
    base_seed = 20260217,
    verbose = TRUE
)
