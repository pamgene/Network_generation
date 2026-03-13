# Golden Score Helper Functions
# 
# This file contains all functions related to golden score calculations.
# It supports both full analysis (with sensitivity data) and kinase-only analysis,
# with both sequential and parallel execution options.
#
# Main functions:
# - make_golden_score_full(): Full analysis with UKA and sensitivity data
# - make_golden_score_kinase(): Kinase-only analysis 
# - make_golden_score(): Unified interface (detects mode automatically)
#
# Author: Extracted from original helper.R and helper_only_kinase_paralell.R files

library(tidyverse)
library(data.table)
library(PCSF)

# =============================================================================
# MAIN GOLDEN SCORE FUNCTIONS
# =============================================================================

#' Unified Golden Score Function
#' 
#' This function automatically detects whether to run full (UKA + sensitivity) 
#' or kinase-only analysis based on the provided parameters.
#' 
#' @param uka UKA data
#' @param sens Optional sensitivity data (if NULL, runs kinase-only mode)
#' @param parallel Logical, run in parallel mode (default: FALSE)
#' @param ... Other parameters passed to specific functions
make_golden_score <- function(uka, sens = NULL, parallel = FALSE, ...) {
  if (is.null(sens)) {
    cat("Running kinase-only golden score analysis...\n")
    make_golden_score_kinase(uka = uka, parallel = parallel, ...)
  } else {
    cat("Running full golden score analysis (UKA + sensitivity)...\n")
    make_golden_score_full(uka = uka, sens = sens, parallel = parallel, ...)
  }
}

#' Full Golden Score Analysis (UKA + Sensitivity data)
#'
#' This function performs the full golden score analysis including both
#' kinase activity (UKA) and sensitivity data.
#'
#' @param uka UKA data frame
#' @param sens Sensitivity data frame  
#' @param del_cells Optional vector of cells to exclude
#' @param control Control condition
#' @param spec_cutoff Specificity score cutoff
#' @param zscore Logical, use z-score transformation
#' @param best_drug_per_target Optional drug selection data
#' @param respath Results path
#' @param score_overlap Logical, calculate overlap scores
#' @param uka_fam Kinase family data
#' @param score_network Logical, calculate network scores
#' @param perc_cutoffs Vector of percentile cutoffs
#' @param nperms_overlap Number of permutations for overlap scoring
#' @param nperms_network Number of permutations for network scoring
#' @param rank_uka_abs Logical, rank UKA by absolute values
#' @param balance Logical, balance dataset sizes
#' @param ppi_network PPI network data
#' @param relative_to What to make network scores relative to
#' @param b Parameter for network scoring
#' @param cs Logical, use cs scoring
#' @param parallel Logical, run in parallel
make_golden_score_full <- function(uka, sens, del_cells = NULL, 
                              control, spec_cutoff, zscore = FALSE, best_drug_per_target = NULL, 
                              respath, score_overlap = TRUE, uka_fam, score_network = TRUE, 
                              perc_cutoffs, nperms_overlap = 500, nperms_network = 50,
                              rank_uka_abs = TRUE, balance = FALSE, ppi_network,
                              relative_to = "n_nodes", b, cs, parallel = FALSE) {

  # Parse input data, filter for comparisons vs control
  sens_parsed <- clean_sens_to_kinograte(sens, control = control, zscore = zscore,
                                           best_drug_per_target = best_drug_per_target)
  uka_parsed <- clean_uka_to_kinograte_full(uka, spec_cutoff = spec_cutoff, control = "RL", cs = cs)
  
  # Find common cell lines in uka and sens
  common_cells <- intersect(unique(uka_parsed$cell_line), unique(sens_parsed$cell_line))
  if (!is.null(del_cells)){
    common_cells <- common_cells[!common_cells %in% del_cells]   
  }
  
  # Create empty results df
  results <- make_golden_results_df(common_cells, perc_cutoffs)
  
  # Initialize metrics df
  all_metrics_df <- data.frame()
  all_obs_metrics_df <- data.frame()
  

  # Main processing loop
  loop_start <- Sys.time()
  logs <- character()
  temp_files <- c()
  
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    # PARALLEL EXECUTION
    for (perc_cutoff in perc_cutoffs) {
      # Parallelize over common_cells using foreach
      cell_results <- foreach::foreach(cell = common_cells, 
                                       .packages = c("dplyr", "readr", "tidyr", "purrr", "rlang", "PCSF", "igraph")) %dopar% {
        process_cell_full(cell, perc_cutoff, uka_parsed, sens_parsed, uka_fam, 
                         score_overlap, score_network, nperms_overlap, nperms_network,
                         respath, spec_cutoff, rank_uka_abs, balance, ppi_network, 
                         relative_to, b)
      }
      
      # Process parallel results
      for (result in cell_results) {
        results <- update_golden_results_row_full(results, result$cell, result$perc_cutoff, result$vals)
        temp_files <- c(temp_files, result$temp_files)
        if (!is.null(result$metrics_df)) all_metrics_df <- bind_rows(all_metrics_df, result$metrics_df)
        if (!is.null(result$obs_metrics)) all_obs_metrics_df <- bind_rows(all_obs_metrics_df, result$obs_metrics)
      }
    }
  } else {
    # SEQUENTIAL EXECUTION
    for (perc_cutoff in perc_cutoffs) {
      for (cell in common_cells) {
        result <- process_cell_full(cell, perc_cutoff, uka_parsed, sens_parsed, uka_fam, 
                                   score_overlap, score_network, nperms_overlap, nperms_network,
                                   respath, spec_cutoff, rank_uka_abs, balance, ppi_network, 
                                   relative_to, b)
        
        results <- update_golden_results_row_full(results, result$cell, result$perc_cutoff, result$vals)
        temp_files <- c(temp_files, result$temp_files)
        if (!is.null(result$metrics_df)) all_metrics_df <- bind_rows(all_metrics_df, result$metrics_df)
        if (!is.null(result$obs_metrics)) all_obs_metrics_df <- bind_rows(all_obs_metrics_df, result$obs_metrics)
      }
    }
  }
  
  return(finalize_golden_score_results(results, all_metrics_df, all_obs_metrics_df, 
                                     temp_files, respath, loop_start, logs))
}


#' Kinase-Only Golden Score Analysis
#'
#' This function performs golden score analysis using only kinase activity (UKA) data,
#' without requiring sensitivity data.
#'
#' @param uka UKA data frame
#' @param spec_cutoff Specificity score cutoff
#' @param respath Results path
#' @param uka_fam Kinase family data (unused in kinase-only mode)
#' @param perc_cutoffs Vector of percentile cutoffs
#' @param nperms_network Number of permutations for network scoring
#' @param rank_uka_abs Logical, rank UKA by absolute values
#' @param ppi_network PPI network data
#' @param b Parameter for network scoring
#' @param cs Logical, use cs scoring
#' @param parallel Logical, run in parallel
make_golden_score_kinase <- function(uka, spec_cutoff, respath, uka_fam, perc_cutoffs, 
                                   nperms_network = 50, rank_uka_abs = TRUE, 
                                   ppi_network, b, cs, parallel = TRUE) {
  
  # Parse UKA
  uka_parsed <- clean_uka_to_kinograte_kinase(uka, cs = cs)
  conditions <- unique(uka_parsed$Sample)

  cat("Dataset has", length(conditions), "conditions for processing\n")

  # Initialize results df, metrics dfs
  results <- initialize_or_read_df(respath = respath, "results")
  all_metrics_df <- initialize_or_read_df(respath = respath, "metrics_permutations")
  all_obs_metrics_df <- initialize_or_read_df(respath = respath, "metrics_observed")

  loop_start <- Sys.time()
  logs <- character()
  temp_files <- c()

  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    # PARALLEL EXECUTION
    for (perc_cutoff in perc_cutoffs) {
      cat("Processing", length(conditions), "conditions in parallel for perc_cutoff", perc_cutoff, "\n")

      condition_results <- foreach::foreach(
        condition = conditions,
        .packages = c("dplyr", "readr", "tidyr", "purrr", "rlang", "PCSF", "igraph"),
        .combine = "list",
        .multicombine = TRUE
      ) %dopar% {
        process_condition_kinase(condition, perc_cutoff, uka_parsed, spec_cutoff, respath,
                               rank_uka_abs, ppi_network, nperms_network, b, cs)
      }

      # Process parallel results
      for (result in condition_results) {
        if (!is.null(result$skip) && result$skip) {
          cat("Skipped condition:", result$condition, "\n")
          next
        }
        if (!is.null(result$error)) {
          cat("Error:", result$error, "\n")
          next
        }
        if (!is.null(result$vals)) {
          results <- update_golden_results_row_kinase(results, result$vals)
        }
        if (!is.null(result$metrics_df)) {
          all_metrics_df <- bind_rows(all_metrics_df, result$metrics_df)
        }
        if (!is.null(result$obs_metrics)) {
          all_obs_metrics_df <- bind_rows(all_obs_metrics_df, result$obs_metrics)
        }
        if (!is.null(result$temp_files)) {
          temp_files <- c(temp_files, result$temp_files)
        }
      }
    }
  } else {
    # SEQUENTIAL EXECUTION
    for (perc_cutoff in perc_cutoffs) {
      for (condition in conditions) {
        # check if data is already in results
        if (perc_cutoff %in% results$perc_cutoff & condition %in% results$condition) {
          print(paste0("Results for perc cutoff: ", perc_cutoff, " and condition: ", condition, " are already calculated. Skipping to next condition."))
          next
        }

        result <- process_condition_kinase(condition, perc_cutoff, uka_parsed, spec_cutoff, respath,
                                         rank_uka_abs, ppi_network, nperms_network, b, cs)
        
        if (!is.null(result$vals)) {
          results <- update_golden_results_row_kinase(results, result$vals)
        }
        if (!is.null(result$metrics_df)) {
          all_metrics_df <- bind_rows(all_metrics_df, result$metrics_df)
        }
        if (!is.null(result$obs_metrics)) {
          all_obs_metrics_df <- bind_rows(all_obs_metrics_df, result$obs_metrics)
        }
        if (!is.null(result$temp_files)) {
          temp_files <- c(temp_files, result$temp_files)
        }
      }
    }
  }
  
  return(finalize_golden_score_results(results, all_metrics_df, all_obs_metrics_df, 
                                     temp_files, respath, loop_start, logs))
}

# =============================================================================
# PROCESSING HELPER FUNCTIONS
# =============================================================================

#' Process single cell for full analysis
process_cell_full <- function(cell, perc_cutoff, uka_parsed, sens_parsed, uka_fam, 
                             score_overlap, score_network, nperms_overlap, nperms_network,
                             respath, spec_cutoff, rank_uka_abs, balance, ppi_network, 
                             relative_to, b) {
  
  # Take top hits for uka and sens for observed score
  uka_filt <- uka_parsed %>% filter(cell_line == cell) %>% uka_top(rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
  sens_filt <- sens_parsed %>% filter(cell_line == cell) %>% sens_top(perc_cutoff, balance)
  # For random sampling, filter for cell
  uka_cell_all <- uka_parsed %>% filter(cell_line == cell)
  
  # Prepare values to update in results csv
  vals <- list(
    n_targets = nrow(sens_filt),
    n_kins = nrow(uka_filt),
    max_sens_value = max(sens_filt$LogFC),
    cell = cell,
    perc_cutoff = perc_cutoff
  )

  temp_files <- c()
  metrics_df <- NULL
  obs_metrics <- NULL
  
  # --- Golden score 1: Overlap ---
  if (score_overlap){
    overlap_res <- compute_golden_overlap_scores(
      uka_filt = uka_filt,
      sens_filt = sens_filt,
      uka_cell_all = uka_cell_all,
      uka_fam = uka_fam,
      cell = cell,
      perc_cutoff = perc_cutoff,
      nperms_overlap = nperms_overlap,
      respath = respath,
      spec_cutoff = spec_cutoff,
      rank_uka_abs = rank_uka_abs
    )
    vals <- c(vals, overlap_res$vals)
    temp_files <- c(temp_files, overlap_res$temp_files)
  } else {
    vals$obs_overlap = NA
    vals$score_sig_overlap = NA
    vals$obs_overlap_fam = NA
    vals$score_sig_overlap_fam = NA
  }

  # --- Golden score 2: Network ---
  if (score_network) {
    network_res <- compute_golden_network_scores_full(
      uka_filt = uka_filt,
      sens_filt = sens_filt,
      uka_cell_all = uka_cell_all,
      cell = cell,
      perc_cutoff = perc_cutoff,
      spec_cutoff = spec_cutoff,
      respath = respath,
      rank_uka_abs = rank_uka_abs, 
      ppi_network = ppi_network,
      nperms_network = nperms_network,
      relative_to = relative_to,
      b = b
    )
    vals <- c(vals, network_res$vals)
    temp_files <- c(temp_files, network_res$temp_files)
    metrics_df <- network_res$metrics_df
    obs_metrics <- network_res$obs_metrics
  } else {
    vals$obs_network = NA
    vals$obs_network_core = NA
    vals$score_sig_network = NA
    vals$score_sig_core_network = NA
  }

  list(cell = cell, perc_cutoff = perc_cutoff, vals = vals, temp_files = temp_files, 
       metrics_df = metrics_df, obs_metrics = obs_metrics)
}

#' Process single condition for kinase-only analysis
process_condition_kinase <- function(condition, perc_cutoff, uka_parsed, spec_cutoff, respath,
                                    rank_uka_abs, ppi_network, nperms_network, b, cs) {
  
  # check if data is already in results (for parallel execution)
  existing_results <- tryCatch({
    if (file.exists(paste0(respath, "/results.csv"))) {
      existing <- read_csv(paste0(respath, "/results.csv"), show_col_types = FALSE)
      if (perc_cutoff %in% existing$perc_cutoff && condition %in% existing$condition) {
        return(list(skip = TRUE, condition = condition))
      }
    }
    FALSE
  }, error = function(e) FALSE)

  if (is.list(existing_results) && existing_results$skip) {
    return(existing_results)
  }

  # Process this condition
  tryCatch({
    # Take top hits for uka for observed score
    uka_filt <- uka_parsed %>%
      filter(Sample == condition) %>%
      uka_top_kinase(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff, cs = cs)

    # For random sampling, filter for condition
    uka_cell_all <- uka_parsed %>% filter(Sample == condition)

    # Prepare values to update in results csv
    vals <- list(
      n_kins = nrow(uka_filt),
      condition = condition,
      spec_cutoff = spec_cutoff,
      perc_cutoff = perc_cutoff
    )

    network_res <- compute_golden_network_scores_kinase(
      uka_filt = uka_filt,
      uka_cell_all = uka_cell_all,
      condition = condition,
      perc_cutoff = perc_cutoff,
      spec_cutoff = spec_cutoff,
      respath = respath,
      rank_uka_abs = rank_uka_abs,
      ppi_network = ppi_network,
      nperms_network = nperms_network,
      b = b,
      cs = cs
    )

    if (!is.null(network_res)) {
      vals <- c(vals, network_res$vals)
      return(list(
        vals = vals,
        metrics_df = network_res$metrics_df,
        obs_metrics = network_res$obs_metrics,
        temp_files = network_res$temp_files,
        condition = condition
      ))
    } else {
      return(list(vals = vals, condition = condition, skip = FALSE))
    }
  }, error = function(e) {
    return(list(error = paste("Error processing condition", condition, ":", e$message), condition = condition))
  })
}

# =============================================================================
# RESULT PROCESSING FUNCTIONS  
# =============================================================================

#' Create empty results data frame for full analysis
make_golden_results_df <- function(common_cells, perc_cutoffs) {
  n_cells <- length(common_cells)
  n_perc_cutoffs <- length(perc_cutoffs)
  data.frame(
    perc_cutoff = rep(perc_cutoffs, each = n_cells),
    cell = rep(common_cells, n_perc_cutoffs),
    n_targets = rep(NA, n_cells * n_perc_cutoffs),
    n_kins = rep(NA, n_cells * n_perc_cutoffs),
    max_sens_value = rep(NA, n_cells * n_perc_cutoffs),
    obs_overlap = rep(NA, n_cells * n_perc_cutoffs),
    score_sig_overlap = rep(NA, n_cells * n_perc_cutoffs),
    obs_overlap_fam = rep(NA, n_cells * n_perc_cutoffs),
    score_sig_overlap_fam = rep(NA, n_cells * n_perc_cutoffs),
    obs_network = rep(NA, n_cells * n_perc_cutoffs),
    score_sig_network = rep(NA, n_cells * n_perc_cutoffs),
    obs_network_core = rep(NA, n_cells * n_perc_cutoffs),
    score_sig_core_network = rep(NA, n_cells * n_perc_cutoffs),
    obs_network_core_inv = rep(NA, n_cells * n_perc_cutoffs),
    score_sig_core_network_inv = rep(NA, n_cells * n_perc_cutoffs)
  )
}

#' Update results row for full analysis
update_golden_results_row_full <- function(results, cell, perc_cutoff, vals) {
  idx <- which(results$perc_cutoff == perc_cutoff & results$cell == cell)
  for (nm in names(vals) ) {
    results[idx, nm] <- vals[[nm]]
  }
  results
}

#' Update results row for kinase-only analysis
update_golden_results_row_kinase <- function(results, vals) {
  vals <- as.data.frame(vals)
  results %>% bind_rows(vals)
}

#' Initialize or read existing dataframe
initialize_or_read_df <- function(respath, dfname) {
  df_filename <- paste0(respath, "/", dfname, ".csv")
  if (file.exists(df_filename)) {
    print(paste0(dfname, " df already exists. New data will be added to the existing one."))
    read_csv(df_filename, show_col_types = FALSE)
  } else {
    data.frame()
  }
}

#' Finalize golden score results 
finalize_golden_score_results <- function(results, all_metrics_df, all_obs_metrics_df, 
                                        temp_files, respath, loop_start, logs) {
  
  # write results & metrics df & plot metrics
  write_csv(results, paste0(respath, "/results.csv"))

  write_csv(all_metrics_df, file.path(respath, "metrics_permutations.csv"))
  write_csv(all_obs_metrics_df, file.path(respath, "metrics_observed.csv"))
  plot_network_metric_hist_facet(all_metrics_df, all_obs_metrics_df, respath)
  
  # measure time
  loop_end <- Sys.time()
  elapsed <- loop_end - loop_start
  print(paste0("Run time: ", elapsed))

  log_time <- file.path(respath, "log_time.txt")
  writeLines(as.character(elapsed), log_time)
  
  logs <- c(logs, paste0("Run time: ", elapsed))
  logs <- c(logs, capture.output(print(results)))

  # Delete all temp files
  for (f in temp_files) {
    if (file.exists(f)) file.remove(f)
  }
  
  return(list(results = results, logs = logs))
}

# =============================================================================
# CORE CALCULATION FUNCTIONS
# =============================================================================

#' Compute overlap and family-wise overlap scores
compute_golden_overlap_scores <- function(uka_filt, sens_filt, uka_cell_all, uka_fam, cell, 
                                        perc_cutoff, nperms_overlap, respath, spec_cutoff, rank_uka_abs) {
  vals <- list()
  # observed score
  score_obs <- overlap_uka_sens(uka_filt$name, sens_filt$name)
  vals$obs_overlap <- round(score_obs, 2)
  # significance score - overlap
  overlap_sig <- pOverlapScore(
    score_obs = score_obs,
    uka_cell_all = uka_cell_all,
    sens_filt = sens_filt,
    cell = cell,
    nsample = nrow(uka_filt),
    nPerms = nperms_overlap,
    spec_cutoff = spec_cutoff,
    rank_uka_abs = rank_uka_abs,
    perc_cutoff = perc_cutoff
  )
  vals$score_sig_overlap <- overlap_sig

  # Write temp file for score_sig_overlap
  temp_overlap_file <- file.path(respath, paste0("temp_overlap_", cell, "_", perc_cutoff, ".txt"))
  writeLines(as.character(overlap_sig), temp_overlap_file)
  temp_files <- c(temp_overlap_file)

  # - Calculate family-wise overlap significance -
  uka_filt_fam <- uka_fam %>% filter(Kinase_Name %in% uka_filt$name) %>% pull(Kinase_family) %>% unique()
  sens_filt_fam <- uka_fam %>% filter(Kinase_Name %in% sens_filt$name) %>% pull(Kinase_family) %>% unique()
  score_obs_fam <- overlap_uka_sens(uka_filt_fam, sens_filt_fam)
  vals$obs_overlap_fam <- round(score_obs_fam, 2)

  # significance score - family
  overlap_sig_fam <- pOverlapScore(
    score_obs = score_obs_fam,
    uka_cell_all = uka_cell_all,
    sens_filt = sens_filt,
    cell = cell,
    nsample = nrow(uka_filt),
    nPerms = nperms_overlap,
    spec_cutoff = spec_cutoff,
    rank_uka_abs = rank_uka_abs,
    perc_cutoff = perc_cutoff,
    uka_fam = uka_fam,
    family_mode = TRUE
  )
  vals$score_sig_overlap_fam <- overlap_sig_fam

  # Write temp file for score_sig_overlap_fam
  temp_overlap_fam_file <- file.path(respath, paste0("temp_overlap_fam_", cell, "_", perc_cutoff, ".txt"))
  writeLines(as.character(overlap_sig_fam), temp_overlap_fam_file)
  temp_files <- c(temp_files, temp_overlap_fam_file)

  list(vals = vals, temp_files = temp_files)
}

#' Compute network scores for full analysis
compute_golden_network_scores_full <- function(uka_filt, sens_filt, uka_cell_all, cell, 
                                             perc_cutoff, spec_cutoff, respath, rank_uka_abs, 
                                             ppi_network, nperms_network, relative_to, b) {
  vals <- list()
  
  # observed score
  network_res_obs <- make_network_and_stats(uka = uka_filt, sens = sens_filt, perc_cutoff = perc_cutoff, 
                                           spec_cutoff = spec_cutoff, res.path = respath, 
                                           condition = cell, write = FALSE, ppi_network = ppi_network,
                                           relative_to = relative_to, b = b)
  if (!is.null(network_res_obs)) {
    score_obs <- network_res_obs$rel_med_path_KT
    score_obs_core <- network_res_obs$rel_med_path_KT_core
    score_obs_core_inv <- network_res_obs$rel_med_path_KT_core_inv

    obs_metrics <- as.data.frame(as.list(network_res_obs), stringsAsFactors = FALSE) %>% 
      mutate(cell = cell, perc_cutoff = perc_cutoff)
    obs_metrics_df_long <- obs_metrics %>% 
      pivot_longer(cols = !c(cell, perc_cutoff), names_to = "metric", values_to = "value")
    
    vals$obs_network <- round(score_obs, 3)
    vals$obs_network_core <- round(score_obs_core, 3)
    vals$obs_network_core_inv <- round(score_obs_core_inv, 3)

    # significance score - network
    perm_res <- pNetworkScore(
      score_obs = score_obs,
      score_obs_core = score_obs_core, 
      score_obs_core_inv = score_obs_core_inv,
      uka_cell_all = uka_cell_all,
      sens_filt = sens_filt,
      cell = cell,
      nsample = nrow(uka_filt),
      nPerms = nperms_network,
      perc_cutoff = perc_cutoff,
      spec_cutoff = spec_cutoff,
      respath = respath,
      condition = cell,
      rank_uka_abs = rank_uka_abs,
      ppi_network = ppi_network,
      relative_to = relative_to,
      b = b
    )

    score_sig_network <- perm_res$network_sig
    score_sig_core_network <- perm_res$network_core_sig
    score_sig_core_network_inv <- perm_res$network_core_inv_sig
    vals$score_sig_network <- score_sig_network
    vals$score_sig_core_network <- score_sig_core_network
    vals$score_sig_core_network_inv <- score_sig_core_network_inv

    # Write temp file for score_sig_network
    temp_network_file <- file.path(respath, paste0("temp_network_", cell, "_", perc_cutoff, ".txt"))
    writeLines(c(as.character(score_sig_network),
                 as.character(score_sig_core_network),
                 as.character(score_sig_core_network_inv)),
               temp_network_file)
    temp_files <- c(temp_network_file)

    list(vals = vals, temp_files = temp_files, metrics_df = perm_res$metrics_df, obs_metrics = obs_metrics_df_long)
  } else {
    vals$obs_network = NA
    vals$obs_network_core = NA
    vals$score_sig_network = NA
    vals$score_sig_core_network = NA
    vals$score_sig_core_network_inv = NA
    list(vals = vals, temp_files = c(), metrics_df = NULL, obs_metrics = NULL)
  }
}

#' Compute network scores for kinase-only analysis
compute_golden_network_scores_kinase <- function(uka_filt, uka_cell_all, condition,
                                                perc_cutoff, spec_cutoff, respath, rank_uka_abs,
                                                ppi_network, nperms_network, b, cs) {
  vals <- list()
  
  # observed score (for kinase-only, no sens data)
  network_res_obs <- make_network_and_stats_kinase(uka = uka_filt, perc_cutoff = perc_cutoff, 
                                                  spec_cutoff = spec_cutoff, res.path = respath, 
                                                  condition = condition, write = FALSE, 
                                                  ppi_network = ppi_network, b = b)
  if (!is.null(network_res_obs)) {
    score_obs <- network_res_obs$rel_med_path_KT
    score_obs_inv <- network_res_obs$rel_med_path_KT_core_inv
    
    obs_metrics <- as.data.frame(as.list(network_res_obs), stringsAsFactors = FALSE) %>% 
      mutate(condition = condition, perc_cutoff = perc_cutoff)
    obs_metrics_df_long <- obs_metrics %>% 
      pivot_longer(cols = !c(condition, perc_cutoff), names_to = "metric", values_to = "value")
    
    vals$obs_network <- round(score_obs, 3)
    vals$obs_network_inv <- round(score_obs_inv, 3)

    # significance score - network (kinase-only permutation)
    perm_res <- pNetworkScore_kinase(
      score_obs = score_obs,
      score_obs_inv = score_obs_inv,
      uka_cell_all = uka_cell_all,
      condition = condition,
      nsample = nrow(uka_filt),
      nPerms = nperms_network,
      perc_cutoff = perc_cutoff,
      spec_cutoff = spec_cutoff,
      respath = respath,
      rank_uka_abs = rank_uka_abs,
      ppi_network = ppi_network,
      b = b
    )

    score_sig_network <- perm_res$network_sig
    score_sig_network_inv <- perm_res$network_inv_sig
    vals$score_sig_network <- score_sig_network
    vals$score_sig_network_inv <- score_sig_network_inv

    # Write temp file for score_sig_network
    temp_network_file <- file.path(respath, paste0("temp_network_", condition, "_", perc_cutoff, ".txt"))
    writeLines(c(as.character(score_sig_network),
                 as.character(score_sig_network_inv)),
               temp_network_file)
    temp_files <- c(temp_network_file)

    list(vals = vals, temp_files = temp_files, metrics_df = perm_res$metrics_df, obs_metrics = obs_metrics_df_long)
  } else {
    vals$obs_network = NA
    vals$obs_network_inv = NA
    vals$score_sig_network = NA
    vals$score_sig_network_inv = NA
    list(vals = vals, temp_files = c(), metrics_df = NULL, obs_metrics = NULL)
  }
}

# =============================================================================
# PERMUTATION AND SCORING FUNCTIONS
# =============================================================================

#' Generalized permutation scoring function for overlap
pOverlapScore <- function(score_obs, uka_cell_all, sens_filt, cell, nsample, nPerms = 1000, 
                         spec_cutoff, rank_uka_abs, perc_cutoff, uka_fam = NULL, family_mode = FALSE) {
  score_perm <- rep(NA, nPerms)
  if (!family_mode) {
    for (i in 1:nPerms){
      # do random permutation of kinases and take the top percentile
      uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
      uka_filt_random <- uka_cell_all %>% uka_top(rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      score_perm[i] <- overlap_uka_sens(uka_filt_random$name, sens_filt$name)
    }
  } else {
    for (i in 1:nPerms){
      uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
      uka_filt_random <- uka_cell_all %>% uka_top(rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      uka_fam_rand <- uka_fam %>% filter(Kinase_Name %in% uka_filt_random$name) %>% pull(Kinase_family) %>% unique()
      sens_fam <- uka_fam %>% filter(Kinase_Name %in% sens_filt$name) %>% pull(Kinase_family) %>% unique()
      score_perm[i] <- overlap_uka_sens(uka_fam_rand, sens_fam)
    }
  }

  overlap_sig = max(sum(score_perm >= score_obs) / nPerms, 1 / nPerms)
  print(paste0("Cell: ", cell, ", overlap score", ifelse(family_mode, " family: ", ": "), overlap_sig))
  return(overlap_sig)
}

#' Network permutation scoring for full analysis
pNetworkScore <- function(score_obs, score_obs_core, score_obs_core_inv, uka_cell_all, sens_filt, cell, nsample, nPerms = 50, 
                         perc_cutoff, spec_cutoff, respath, condition,
                         rank_uka_abs, ppi_network, relative_to, b) {  
  
  score_perm <- rep(NA, nPerms)
  score_perm_core <- rep(NA, nPerms)
  score_perm_core_inv <- rep(NA, nPerms)
  # Replace for loop with lapply
  perm_results <- lapply(1:nPerms, function(i) {
    uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
    uka_filt_random <- uka_cell_all %>% uka_top(rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
    print(paste0("Cell: ", cell, ", Network permutation: ", i, "..."))
    res <- make_network_and_stats(uka = uka_filt_random, sens = sens_filt, perc_cutoff = perc_cutoff, 
                                  spec_cutoff = spec_cutoff, res.path = respath, condition = cell, write = FALSE, 
                                  ppi_network = ppi_network, relative_to = relative_to, b = b)
    if (!is.null(res)) {
      metrics <- as.data.frame(as.list(res), stringsAsFactors = FALSE) %>% mutate(nperm = i)
      score <- res$rel_med_path_KT
      score_core <- res$rel_med_path_KT_core
      score_core_inv <- res$rel_med_path_KT_core_inv
    } else {
      metrics <- data.frame(
        density = NA, clustering = NA, modularity = NA, assortativity = NA,
        rel_med_path_KT = NA, rel_med_path_KT_core = NA, 
        rel_med_path_all_nodes = NA, rel_med_path_all_nodes_core = NA, 
        rel_med_path_KT_vs_all_nodes = NA, 
        rel_med_path_KT_vs_all_nodes_core = NA,
        rel_med_path_KT_core_inv = NA,
        rel_med_path_all_core_inv = NA,
        rel_med_path_KT_vs_all_nodes_core_inv = NA,
        nperm = NA
      )
      score <- 1000
      score_core <- 1000
      score_core_inv <- 1000
    }
    list(score = score, score_core = score_core, score_core_inv = score_core_inv, metrics = metrics)
  })

  score_perm <- sapply(perm_results, function(x) x$score)
  score_perm_core <- sapply(perm_results, function(x) x$score_core)
  score_perm_core_inv <- sapply(perm_results, function(x) x$score_core_inv)
  metrics_list <- lapply(perm_results, function(x) x$metrics)

  network_sig <- max(sum(score_perm <= score_obs) / nPerms, 1 / nPerms)
  network_core_sig <- max(sum(score_perm_core <= score_obs_core) / nPerms, 1 / nPerms)
  network_core_inv_sig <- max(sum(score_perm_core_inv >= score_obs_core_inv) / nPerms, 1 / nPerms)
  print(paste0("Cell: ", cell, ", network score: ", network_sig))
  print(paste0("Cell: ", cell, ", network core score: ", network_core_sig))
  print(paste0("Cell: ", cell, ", network core inv score: ", network_core_inv_sig))
  metrics_df <- do.call(bind_rows, metrics_list)
  metrics_df_long <- metrics_df %>% 
    pivot_longer(cols = !'nperm', names_to = "metric", values_to = "value") %>% 
    mutate(cell = cell, perc_cutoff = perc_cutoff) %>%
    select(cell, perc_cutoff, nperm, metric, value)
  return(list(
    network_sig = network_sig, 
    network_core_sig = network_core_sig, 
    network_core_inv_sig = network_core_inv_sig,
    metrics_df = metrics_df_long
  ))
}

#' Network permutation scoring for kinase-only analysis
pNetworkScore_kinase <- function(score_obs, score_obs_inv, uka_cell_all, condition, nsample, nPerms = 50, 
                                perc_cutoff, spec_cutoff, respath,
                                rank_uka_abs, ppi_network, b) {  
  
  score_perm <- rep(NA, nPerms)
  score_perm_inv <- rep(NA, nPerms)
  
  # Replace for loop with lapply
  perm_results <- lapply(1:nPerms, function(i) {
    uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
    uka_filt_random <- uka_cell_all %>% uka_top_kinase(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
    print(paste0("Condition: ", condition, ", Network permutation: ", i, "..."))
    res <- make_network_and_stats_kinase(uka = uka_filt_random, perc_cutoff = perc_cutoff, 
                                       spec_cutoff = spec_cutoff, res.path = respath, condition = condition, 
                                       write = FALSE, ppi_network = ppi_network, b = b)
    if (!is.null(res)) {
      metrics <- as.data.frame(as.list(res), stringsAsFactors = FALSE) %>% mutate(nperm = i)
      score <- res$rel_med_path_KT
      score_inv <- res$rel_med_path_KT_core_inv
    } else {
      metrics <- data.frame(
        density = NA, clustering = NA, modularity = NA, assortativity = NA,
        rel_med_path_KT = NA, rel_med_path_KT_core_inv = NA,
        nperm = NA
      )
      score <- 1000
      score_inv <- 1000
    }
    list(score = score, score_inv = score_inv, metrics = metrics)
  })

  score_perm <- sapply(perm_results, function(x) x$score)
  score_perm_inv <- sapply(perm_results, function(x) x$score_inv)
  metrics_list <- lapply(perm_results, function(x) x$metrics)

  network_sig <- max(sum(score_perm <= score_obs) / nPerms, 1 / nPerms)
  network_inv_sig <- max(sum(score_perm_inv >= score_obs_inv) / nPerms, 1 / nPerms)
  print(paste0("Condition: ", condition, ", network score: ", network_sig))
  print(paste0("Condition: ", condition, ", network inv score: ", network_inv_sig))
  metrics_df <- do.call(bind_rows, metrics_list)
  metrics_df_long <- metrics_df %>% 
    pivot_longer(cols = !'nperm', names_to = "metric", values_to = "value") %>% 
    mutate(condition = condition, perc_cutoff = perc_cutoff) %>%
    select(condition, perc_cutoff, nperm, metric, value)
  return(list(
    network_sig = network_sig, 
    network_inv_sig = network_inv_sig,
    metrics_df = metrics_df_long
  ))
}

# =============================================================================
# DATA PROCESSING AND UTILITY FUNCTIONS
# =============================================================================

#' Clean UKA data for full analysis (with sensitivity data)
clean_uka_to_kinograte_full <- function(uka, control, spec_cutoff = 0, del_cell = NULL, cs = FALSE){
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  stat_col <- if (cs) "Kinase Statistic" else "Median Kinase Statistic"
  
  # --- Handle both "X vs control" and "control vs X" in contrast, flip sign if control is left ---
  uka_clean <- uka %>%
    clean_tercen_columns() %>%
    filter(.data[[finalscore_col]] > spec_cutoff) %>%
    filter(str_detect(contrast, paste0(" vs ", control)) | str_detect(contrast, paste0(control, " vs "))) %>%
    mutate(contrast = gsub(contrast, pattern = "-", replacement = "")) %>%
    mutate(
      control_left = str_detect(contrast, paste0("^", control, " vs ")),
      cell_line = ifelse(control_left,
                         str_replace(contrast, paste0(control, " vs "), ""),
                         str_replace(contrast, paste0(" vs ", control), "")),
      MKS = ifelse(control_left, -.data[[stat_col]], .data[[stat_col]])
    ) %>%
    dplyr::select('cell_line', 'Kinase Name', all_of(stat_col), all_of(finalscore_col)) %>%
    dplyr::rename("uniprotname" = "Kinase Name", 'LogFC' = !!stat_col, 
                  'fscore' = !!finalscore_col) %>%
    mutate(cell_line = gsub(cell_line, pattern = "-", replacement = "")) %>%
    distinct()
  
  return(uka_clean)
}

#' Clean UKA data for kinase-only analysis
clean_uka_to_kinograte_kinase <- function(uka, cs = FALSE) {
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  stat_col <- if (cs) "Kinase Statistic" else "Median Kinase Statistic"

  uka_clean <- uka %>%
    clean_tercen_columns() %>%
    dplyr::select("Sgroup_contrast", "Kinase Name", all_of(stat_col), all_of(finalscore_col)) %>%
    dplyr::rename(
      "uniprotname" = "Kinase Name", "LogFC" = !!stat_col,
      "fscore" = !!finalscore_col
    ) %>%
    distinct()

  return(uka_clean)
}

#' Clean sensitivity data for kinograte analysis
clean_sens_to_kinograte <- function(sens, control, zscore = FALSE, del_cell = NULL, best_drug_per_target){
  # if control is specified, it creates fold change.
  # if control is NULL, and zscore is True, creates z-score
  sens_filt <- sens %>% 
    dplyr::rename('cell_line' = 'CELL_LINE_NAME',
                  "uniprotname" = "TARGET_1")

  # if best drug per target is given, keep only that drug per target (based on data availability and sd)
  if (!is.null(best_drug_per_target)){
    sens_filt <- sens_filt %>%
      filter(DRUG_ID %in% best_drug_per_target$DRUG_ID) 
  }
  # If control is given, make LogFC
  if (!is.null(control) & zscore == FALSE){
    control_df <- sens_filt %>% filter(cell_line == control)
    sens_clean <- sens_filt %>%
      left_join(control_df, by = c("uniprotname"), suffix = c("", ".control")) %>%
      mutate(
        LogFC = LN_IC50 - LN_IC50.control
      ) %>%
      filter(cell_line != control, 
             !is.na(LogFC)) %>%
      dplyr::select(cell_line, uniprotname, LogFC)
  } else if (is.null(control) & zscore == FALSE){
    # use median sensitivity data based on LogFC
    sens_clean <- sens_filt
  } else if (is.null(control) & zscore == TRUE){
    sens_clean <- sens_filt %>%
      pivot_wider(id_cols = c("cell_line"), names_from = uniprotname, values_from = LN_IC50) %>%
      column_to_rownames("cell_line") %>%
      as.matrix() %>% scale() %>%
      as.data.frame() %>% rownames_to_column("cell_line") %>%
      pivot_longer(cols = !cell_line, names_to = "uniprotname", values_to = "LogFC")
  } 
  
  if(!is.null(del_cell)){
    sens_clean <- sens_clean %>% 
      filter(!cell_line %in% del_cell)
  }
  
  return(sens_clean)
}

#' Clean tercen column names
clean_tercen_columns <- function(df) {
  cols <- colnames(df)
  split <- str_split(cols, pattern = "\\.")
  cols <- sapply(split, tail, 1)
  colnames(df) <- cols
  return(df)
}

#' Calculate overlap between UKA and sensitivity data
overlap_uka_sens <- function(uka, sens){
  overlap <- length(intersect(uka, sens))
  # return relative overlap
  overlap/max(length(sens), 1)
}

#' Get top UKA hits for full analysis
uka_top <- function(uka, rank_uka_abs = TRUE, perc_cutoff){
  if (rank_uka_abs) {
    uka_rank <- percentile_rank_fast(uka, uniprotname, LogFC)
  } else {
    uka_rank <- percentile_rank_pg_noabs(uka, symbol = uniprotname, metric = LogFC, rank_lowest_highest = FALSE)
  }
  uka_rank %>% dplyr::filter(score >= perc_cutoff) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(type = 'Kinase') %>% 
      dplyr::select(name, prize = score, type, LogFC)
}

#' Get top UKA hits for kinase-only analysis  
uka_top_kinase <- function(uka, spec_cutoff, rank_uka_abs = TRUE, perc_cutoff, cs){
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  
  uka_filt <- uka %>% filter(.data[[finalscore_col]] > spec_cutoff)
  
  if (rank_uka_abs) {
    uka_rank <- percentile_rank_fast(uka_filt, uniprotname, LogFC)
  } else {
    uka_rank <- percentile_rank_pg_noabs(uka_filt, symbol = uniprotname, metric = LogFC, rank_lowest_highest = FALSE)
  }
  uka_rank %>% dplyr::filter(score >= perc_cutoff) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(type = 'Kinase') %>% 
      dplyr::select(name, prize = score, type, LogFC)
}

#' Get top sensitivity hits
sens_top <- function(sens, perc_cutoff, balance = TRUE){
  if (balance){
    perc_cutoff <- perc_cutoff - 0.2
  } 
  sens_rank <- percentile_rank_pg_noabs(sens, uniprotname, LogFC, rank_lowest_highest = TRUE)
  sens_rank %>% dplyr::filter(score >= perc_cutoff) %>% dplyr::ungroup() %>% 
    dplyr::mutate(type = "Sensitivity") %>% 
    dplyr::select(name, prize = score, type, LogFC)
}

#' Fast percentile rank calculation
percentile_rank_fast <- function(df, symbol, metric, desc = FALSE) {
  # Convert to data.table by reference
  dt <- as.data.table(df)
  
  # Get column names as strings
  sym_col <- deparse(substitute(symbol))
  met_col <- deparse(substitute(metric))
  
  # Filter NA, rename, and get unique in one operation
  dt <- dt[!is.na(get(met_col))]
  setnames(dt, sym_col, "name", skip_absent = TRUE)
  dt <- unique(dt, by = "name")
  
  # Calculate percentile rank
  abs_metric <- abs(dt[[met_col]])
  
  if (desc) {
    dt[, score := frank(-abs_metric, ties.method = "average") / .N]
  } else {
    dt[, score := frank(abs_metric, ties.method = "average") / .N]
  }
  
  # Sort and return
  setorder(dt, -score)
  dt[]
}

#' Percentile rank calculation without abs (from kinograte_PG_core.R)
percentile_rank_pg_noabs <- function(df, symbol, metric, rank_lowest_highest = FALSE) {
  # Rename the symbol column to "name"
  df <- df %>% 
    rename(name = !!enquo(symbol)) %>%
    filter(!is.na(!!enquo(metric))) %>%
    distinct(name, .keep_all = TRUE)
  
  if (rank_lowest_highest) {
    df <- dplyr::mutate(df, score = dplyr::percent_rank(dplyr::desc({{ metric }})))
  } else {
    df <- dplyr::mutate(df, score = dplyr::percent_rank({{ metric }}))
  }
  
  dplyr::arrange(df, desc(score)) 
}

# =============================================================================
# VISUALIZATION AND ANALYSIS FUNCTIONS
# =============================================================================

#' Plot network metric histograms
plot_network_metric_hist_facet <- function(all_metrics_df, all_obs_df, res.path) {
  dir.create(file.path(res.path, "Network_metrics_plots"), showWarnings = FALSE, recursive = TRUE)
  make_plot <- function(metric_i, cutoff_i) {
    
    df1_sub <- all_metrics_df %>%
      filter(metric == metric_i, perc_cutoff == cutoff_i,
             is.finite(value),
             !is.na(value))
    df2_sub <- all_obs_df %>%
      filter(metric == metric_i, perc_cutoff == cutoff_i,
             is.finite(value),
             !is.na(value))
    
    # Check if either dataframe is empty after filtering
    if (nrow(df1_sub) == 0 || nrow(df2_sub) == 0) {
      # Create empty plot with message
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No valid data available", size = 6) +
        theme_void() +
        labs(title = paste0("Metric: ", metric_i, " | Cutoff: ", cutoff_i))
      
      # Save empty plot
      fname <- paste0(res.path, "/Network_metrics_plots/", metric_i, "_cutoff", cutoff_i, ".png")
      ggsave(fname, p, width = 7, height = 5, dpi = 300)
      
      return(p)
    }
    
    # Determine grouping column
    grouping_col <- if ("cell" %in% colnames(df1_sub)) "cell" else "condition"
    
    p <- ggplot(df1_sub, aes(x = value)) +
      geom_histogram(bins = 30, fill = "grey70", color = "black") +
      geom_vline(
        data = df2_sub,
        aes(xintercept = value), color = "red", linetype = "dashed") +
      facet_wrap(as.formula(paste("~", grouping_col)), scales = "free_y") +
      theme_bw() +
      theme(
        legend.position = "top",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold")
      ) +
      labs(
        title = paste0("Metric: ", metric_i, " | Cutoff: ", cutoff_i),
        x = "Permutation values",
        y = "Count"
      )
    
    # Save plot
    fname <- paste0(res.path, "/Network_metrics_plots/", metric_i, "_cutoff", cutoff_i, ".png")
    ggsave(fname, p, width = 7, height = 5, dpi = 300)
    
    return(p)
  }
  
  combos <- all_metrics_df %>%
    distinct(metric, perc_cutoff)
  plots <- pmap(
    list(combos$metric, combos$perc_cutoff),
    make_plot
  )
}

# =============================================================================
# PARAMETER MANAGEMENT FUNCTIONS
# =============================================================================

#' Capture and organize function parameters
capture_params <- function(...) {
  dots <- enquos(...)
  param_names <- names(dots)
  # Handle unnamed arguments
  if (is.null(param_names)) {
    param_names <- rep("", length(dots))
  }
  
  # Replace blank names with original variable names
  symbol_names <- map_chr(dots, rlang::as_label)
  param_names <- ifelse(param_names == "" | is.na(param_names), symbol_names, param_names)
  
  # Evaluate objects
  values <- map(dots, rlang::eval_tidy)
  params <- set_names(values, param_names)
  
  # Attach original variable names
  attr(params, "symbol_names") <- set_names(symbol_names, param_names)
  
  # Create parameter folder based on analysis type
  if (!is.null(params$sens)) {
    # Full analysis folder structure
    uka_name <- symbol_names[["uka"]]
    sens_name <- symbol_names[["sens"]]
    ppi_network_name <- symbol_names[["ppi_network"]]
    param_folder <- file.path(params$respath, uka_name, 
          paste0(
            sens_name, "_thrs_", paste0(params$perc_cutoffs, collapse = "-"),
            "_", ppi_network_name,
            '_nperms', params$nperms_network,
            "_blnc", as.integer(params$balance), "_rel2", params$relative_to,
            "_ukaabs", as.integer(params$rank_uka_abs), "_b", params$b
    ))
  } else {
    # Kinase-only folder structure
    uka_name <- symbol_names[["uka"]]
    ppi_network_name <- symbol_names[["ppi_network"]]
    param_folder <- file.path(
      params$respath, 
      paste0(
        "thrs_", paste0(params$perc_cutoffs, collapse = "-"),
        "_spec", params$spec_cutoff,
        "_", ppi_network_name,
        "_nperms", params$nperms_network,
        "_ukaabs", as.integer(params$rank_uka_abs), "_b", params$b,
        ifelse(!is.null(params$art_nodes), "_art", "")
      )
    )
  }
  
  print(paste0("param folder: ", param_folder))
  if (!dir.exists(param_folder)) dir.create(param_folder, recursive = TRUE)
  params$respath <- param_folder
  save_params(params, respath = param_folder)
  return(params)
}

#' Save parameters to file
save_params <- function(params, respath){
  symbol_names <- attr(params, "symbol_names")
  
  df <- tibble(
    parameter = names(symbol_names),
    value = unname(symbol_names)
  )
  df <- df %>%
    mutate(value = gsub('\\"', '"', value)) %>%   # replace escaped quotes 
    mutate(value = gsub('c\\(|\\)', '', value)) %>% # remove c( and )
    mutate(value = gsub('^"|"$', '', value)) %>%       # remove leading/trailing quotes 
    mutate(value = str_remove_all(value, '\\\\')) %>%  # remove all backslashes
    mutate(value = str_remove_all(value, '"')) %>%
    filter(!parameter %in% c('uka_fam','respath'))
  
  write_csv(df, paste0(respath, "/params.csv"))
  print(df)
}

# =============================================================================
# RESULT ANALYSIS FUNCTIONS
# =============================================================================

#' Reconstruct partial results from temp files in a folder
reconstruct_golden_score_from_temp <- function(respath, cells, perc_cutoffs) {
  df <- expand.grid(cell = cells, perc_cutoff = perc_cutoffs, stringsAsFactors = FALSE)
  df$score_sig_overlap <- NA_real_
  df$score_sig_overlap_fam <- NA_real_
  df$score_sig_network <- NA_real_
  for (i in seq_len(nrow(df))) {
    cell <- df$cell[i]
    perc <- df$perc_cutoff[i]
    # Overlap
    f_overlap <- file.path(respath, paste0("temp_overlap_", cell, "_", perc, ".txt"))
    if (file.exists(f_overlap)) {
      val <- suppressWarnings(as.numeric(readLines(f_overlap, warn = FALSE)))
      if (!is.na(val)) df$score_sig_overlap[i] <- val
    }
    # Overlap fam
    f_overlap_fam <- file.path(respath, paste0("temp_overlap_fam_", cell, "_", perc, ".txt"))
    if (file.exists(f_overlap_fam)) {
      val <- suppressWarnings(as.numeric(readLines(f_overlap_fam, warn = FALSE)))
      if (!is.na(val)) df$score_sig_overlap_fam[i] <- val
    }
    # Network
    f_network <- file.path(respath, paste0("temp_network_", cell, "_", perc, ".txt"))
    if (file.exists(f_network)) {
      val <- suppressWarnings(as.numeric(readLines(f_network, warn = FALSE)))
      if (!is.na(val)) df$score_sig_network[i] <- val
    }
  }
  df
}

#' Find missing cell/perc_cutoff combinations (where any temp file is missing)
find_missing_golden_score_combinations <- function(respath, cells, perc_cutoffs) {
  df <- expand.grid(cell = cells, perc_cutoff = perc_cutoffs, stringsAsFactors = FALSE)
  missing_idx <- c()
  for (i in seq_len(nrow(df))) {
    cell <- df$cell[i]
    perc <- df$perc_cutoff[i]
    f_overlap <- file.path(respath, paste0("temp_overlap_", cell, "_", perc, ".txt"))
    f_overlap_fam <- file.path(respath, paste0("temp_overlap_fam_", cell, "_", perc, ".txt"))
    f_network <- file.path(respath, paste0("temp_network_", cell, "_", perc, ".txt"))
    if (!(file.exists(f_overlap) && file.exists(f_overlap_fam) && file.exists(f_network))) {
      missing_idx <- c(missing_idx, i)
    }
  }
  df[missing_idx, , drop = FALSE]
}

#' Create visualization plots for golden score results
golden_results_to_barplot <- function(results_list, result_path = NULL, perccutoff = 0.7) {
  
  if (!is.null(result_path) && !dir.exists(paste0(result_path, "/Score_plots"))) {
    dir.create(paste0(result_path, "/Score_plots"), recursive = TRUE)
  }
  
  # If file paths are given, read them
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = FALSE, id = "path") else x
  }) %>%
    bind_rows()
  
  # Detect analysis type and create appropriate plots
  if ("condition" %in% colnames(results_df)) {
    # Kinase-only analysis plotting
    golden_results_to_barplot_kinase(results_df, result_path)
  } else {
    # Full analysis plotting  
    golden_results_to_barplot_full(results_df, perccutoff)
  }
}

#' Barplot for kinase-only results
golden_results_to_barplot_kinase <- function(results_df, result_path) {
  results_df <- results_df %>% select(-`path...2`) %>% rename("path" = "path...1")
  
  results_df_clean <- results_df %>%
    mutate(path = sub(".*Onlykin_fscore1.3\\/", "", path)) %>%
    mutate(path = sub("\\/results_notna.csv.*", "", path)) %>%
    separate_wider_delim(cols = "path", delim = "/", names = c("path", "params")) %>%
    filter(params == "thrs_0.01_ppi_networkv12_nperms30_ukaabs1_b1")

  p <- ggplot(results_df_clean, aes(x = condition, y = score_sig_network)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    labs(
      x = "Condition", y = "score (p-value)",
      title = unique(results_df_clean$params)
    ) +
    facet_wrap(~path, ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  if (!is.null(result_path)) {
    ggsave(paste0(result_path, "/Score_plots/Barplot_all_scores", unique(results_df_clean$params), ".png"),
           p, w = 18, h = 20, unit = "cm")
  }
  
  return(p)
}

#' Barplot for full analysis results
golden_results_to_barplot_full <- function(results_df, perccutoff) {
  # Golden score columns
  my_cols <- c("score_sig_overlap", "score_sig_overlap_fam", "score_sig_network", 
               "score_sig_core_network", "score_sig_core_network_inv")
  
  results_df_long <- results_df %>%
    filter(perc_cutoff == perccutoff) %>% 
    pivot_longer(cols = all_of(my_cols), names_to = "Metric", values_to = "value")
  
  results_df_network <- results_df_long %>% 
    filter(Metric %in% c("score_sig_network", "score_sig_core_network", "score_sig_core_network_inv"))
  
  p <- ggplot(results_df_network, aes(x = cell, y = value, fill = Metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    labs(x = "Cell", y = "score (p-value)",
         title = paste0("Percentile cutoff: ", perccutoff)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p)
}

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

# Note: The following functions (make_network_and_stats, make_network_and_stats_kinase) 
# are dependencies that should exist in other R files (likely kinograte_PG.R or similar)
# They are referenced but not implemented here as they are part of the network analysis core
