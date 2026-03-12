plot_network_metric_hist_facet <- function(all_metrics_df, all_obs_df, res.path) {
  dir.create(file.path(res.path, "Network_metrics_plots"), showWarnings = FALSE, recursive = TRUE)
  make_plot <- function(metric_i, cutoff_i) {
    df1_sub <- all_metrics_df %>%
      filter(
        metric == metric_i, perc_cutoff == cutoff_i,
        is.finite(value),
        !is.na(value)
      )
    df2_sub <- all_obs_df %>%
      filter(
        metric == metric_i, perc_cutoff == cutoff_i,
        is.finite(value),
        !is.na(value)
      )
    Sys.getenv()
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

    # browser()
    p <- ggplot(df1_sub, aes(x = value)) +
      geom_histogram(bins = 30, fill = "grey70", color = "black") +
      geom_vline(
        data = df2_sub,
        aes(xintercept = value), color = "red", linetype = "dashed"
      ) +
      facet_wrap(~condition, scales = "free_y") +
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



pNetworkScore <- function(score_obs, score_obs_inv, uka_cell_all, condition, nsample, nPerms = DFT.PERMS,
                          perc_cutoff, spec_cutoff, respath,
                          rank_uka_abs, ppi_network, b, cs) {
  score_perm <- rep(NA, nPerms)
  score_perm_inv <- rep(NA, nPerms)
  metrics_list <- vector("list", nPerms)
  for (i in 1:nPerms) {
    uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
    uka_filt_random <- uka_cell_all %>% uka_top(
      spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs,
      perc_cutoff = perc_cutoff, cs = cs
    )
    print(paste0("Comparison: ", condition, ", Network permutation: ", i, "..."))
    res <- make_network_and_stats(
      uka = uka_filt_random, perc_cutoff = perc_cutoff,
      spec_cutoff = spec_cutoff, res.path = respath, condition = condition,
      write = FALSE, ppi_network = ppi_network, b = b
    )

    if (!is.null(res)) {
      metrics_list[[i]] <- as.data.frame(as.list(res), stringsAsFactors = FALSE) %>% mutate(nperm = i)
      score_perm[i] <- res$rel_med_path
      score_perm_inv[i] <- res$rel_med_path_inv
    } else {
      # Fill with NA for all expected columns
      metrics_list[[i]] <- data.frame(
        density = NA, clustering = NA,
        rel_med_path = NA, rel_med_path_inv = NA,
        nperm = NA
      )
      score_perm[i] <- 1000
      score_perm_inv[i] <- 1000
    }
  }
  nPerms_notna <- length(score_perm[!is.na(score_perm)])
  network_sig_score <- max(sum(score_perm <= score_obs, na.rm = T) / nPerms_notna, 1 / nPerms_notna)
  network_inv_sig_score <- max(sum(score_perm_inv >= score_obs_inv, na.rm = T) / nPerms_notna, 1 / nPerms_notna)
  print(paste0("Condition: ", condition, ", network score: ", network_sig_score))
  print(paste0("Condition: ", condition, ", network inv score: ", network_inv_sig_score))
  metrics_df <- do.call(bind_rows, metrics_list)
  metrics_df_long <- metrics_df %>%
    pivot_longer(cols = !"nperm", names_to = "metric", values_to = "value") %>%
    mutate(condition = condition, perc_cutoff = perc_cutoff) %>%
    select(condition, perc_cutoff, nperm, metric, value)
  return(list(
    network_sig = network_sig_score,
    network_inv_sig = network_inv_sig_score,
    metrics_df = metrics_df_long
  ))
}


uka_top <- function(uka, spec_cutoff, rank_uka_abs = T, perc_cutoff, cs) {
  # filter final score
  uka <- uka %>% filter(fscore >= spec_cutoff)
  if (rank_uka_abs) {
    uka_rank <- percentile_rank_fast(uka, uniprotname, LogFC) # kinograte::percentile_rank of abs LogFC of kinases
  } else {
    uka_rank <- percentile_rank_pg_noabs(uka, symbol = uniprotname, metric = LogFC, rank_lowest_highest = F)
  }
  uka_rank %>%
    dplyr::filter(score >= perc_cutoff) %>% # & score >= perc_cutoff
    dplyr::ungroup() %>%
    dplyr::mutate(type = "Kinase") %>%
    dplyr::select(name, prize = score, type, LogFC)
}

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



initialize_or_read_df <- function(respath, dfname) {
  df_filename <- paste0(respath, "/", dfname, ".csv")
  if (file.exists(df_filename)) {
    print(paste0(dfname, " df already exists. New data will be added to the existing one."))
    read_csv(df_filename, show_col_types = FALSE)
  } else {
    data.frame()
  }
}
# Helper to update a row in the results data frame
update_golden_results_row <- function(results, vals) {
  vals <- as.data.frame(vals)
  results %>% bind_rows(vals)
}


# Helper to compute network scores and manage temp files


compute_golden_network_scores <- function(uka_filt, uka_cell_all, condition,
                                          perc_cutoff, spec_cutoff, respath, rank_uka_abs,
                                          ppi_network, nperms_network, b, cs) {
  vals <- list()
  # observed network score and metrics
  obs_res <- make_network_and_stats(
    uka = uka_filt,
    perc_cutoff = perc_cutoff, spec_cutoff = spec_cutoff,
    res.path = respath, condition = condition, write = F,
    ppi_network = ppi_network, b = b
  )
  obs_network <- obs_res$rel_med_path
  obs_network_inv <- obs_res$rel_med_path_inv
  obs_metrics_df_long <- as.data.frame(obs_res) %>%
    pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
    mutate(condition = condition, perc_cutoff = perc_cutoff) %>%
    select(condition, perc_cutoff, metric, value)
  vals$obs_network <- obs_network
  vals$obs_network_inv <- obs_network_inv

  # If observed result is NULL or NA, return NAs for permutation results
  if (all(is.na(obs_res))) {
    print(paste0("Network path of: ", condition, ", percentile cutoff: ", perc_cutoff, " has returned NULL"))
    return(list(
      vals = c(vals, score_sig_network = NA, score_sig_network_inv = NA),
      temp_files = NA,
      metrics_df = data.frame(condition = NA, perc_cutoff = NA, nperm = NA, metric = NA, value = NA),
      obs_metrics = obs_metrics_df_long
    ))
  }

  # collect metrics from permutations
  perm_res <- pNetworkScore(
    score_obs = obs_network,
    score_obs_inv = obs_network_inv,
    uka_cell_all = uka_cell_all,
    condition = condition,
    nsample = nrow(uka_filt),
    nPerms = nperms_network,
    perc_cutoff = perc_cutoff,
    spec_cutoff = spec_cutoff,
    respath = respath,
    rank_uka_abs = rank_uka_abs,
    ppi_network = ppi_network,
    b = b,
    cs = cs
  )
  score_sig_network <- perm_res$network_sig
  score_sig_network_inv <- perm_res$network_inv_sig
  vals$score_sig_network <- score_sig_network
  vals$score_sig_network_inv <- score_sig_network_inv

  # Write temp file for score_sig_network
  temp_network_file <- file.path(respath, paste0("temp_network_", condition, "_", perc_cutoff, ".txt"))
  writeLines(
    c(
      as.character(score_sig_network),
      as.character(score_sig_network_inv)
    ),
    temp_network_file
  )
  temp_files <- c(temp_network_file)

  list(vals = vals, temp_files = temp_files, metrics_df = perm_res$metrics_df, obs_metrics = obs_metrics_df_long)
}

make_golden_score <- function(uka,
                              spec_cutoff, respath,
                              uka_fam, perc_cutoffs, nperms_network = 50,
                              rank_uka_abs = T, ppi_network = ppi_networkv12, b, cs,
                              parallel_conditions = TRUE) {
  # Parse UKA
  uka_parsed <- clean_uka_to_kinograte(uka, cs = cs)
  conditions <- unique(uka_parsed$Sample)

  cat("Dataset has", length(conditions), "conditions for parallel processing\n")

  # Initialize results df, metrics dfs
  results <- initialize_or_read_df(respath = respath, "results")
  all_metrics_df <- initialize_or_read_df(respath = respath, "metrics_permutations")
  all_obs_metrics_df <- initialize_or_read_df(respath = respath, "metrics_observed")

  loop_start <- Sys.time()
  logs <- character() # Collect log messages
  temp_files <- c()

  for (perc_cutoff in perc_cutoffs) {
    # PARALLEL VERSION: Process conditions in parallel instead of sequential
    if (parallel_conditions && requireNamespace("foreach", quietly = TRUE)) {
      cat("Processing", length(conditions), "conditions in parallel for perc_cutoff", perc_cutoff, "\n")

      condition_results <- foreach(
        condition = conditions,
        .packages = c("dplyr", "readr", "tidyr", "purrr", "rlang", "PCSF", "igraph"),
        .combine = "list",
        .multicombine = TRUE
      ) %dopar% {
        # check if data is already in results (using exists conditions from main thread)
        existing_results <- tryCatch(
          {
            if (file.exists(paste0(respath, "/results.csv"))) {
              existing <- read_csv(paste0(respath, "/results.csv"), show_col_types = FALSE)
              if (perc_cutoff %in% existing$perc_cutoff && condition %in% existing$condition) {
                return(list(skip = TRUE, condition = condition))
              }
            }
            FALSE
          },
          error = function(e) FALSE
        )

        if (is.list(existing_results) && existing_results$skip) {
          return(existing_results)
        }

        # Process this condition
        tryCatch(
          {
            # Take top hits for uka for observed score
            uka_filt <- uka_parsed %>%
              filter(Sample == condition) %>%
              uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff, cs = cs)

            # For random sampling, filter for cell
            uka_cell_all <- uka_parsed %>% filter(Sample == condition)

            # Prepare values to update in results csv
            vals <- list(
              n_kins = nrow(uka_filt),
              condition = condition,
              spec_cutoff = spec_cutoff,
              perc_cutoff = perc_cutoff
            )

            network_res <- compute_golden_network_scores(
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
          },
          error = function(e) {
            return(list(error = paste("Error processing condition", condition, ":", e$message), condition = condition))
          }
        )
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
          results <- update_golden_results_row(results, result$vals)
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
    } else {
      # SEQUENTIAL VERSION: Original sequential loop over conditions
      for (condition in conditions) {
        # check if data is already in results
        if (perc_cutoff %in% results$perc_cutoff & condition %in% results$condition) {
          print(paste0("Results for perc cutoff: ", perc_cutoff, " and condition: ", condition, " are already calculated. Skipping to next condition."))
          next
        }

        # Take top hits for uka for observed score
        uka_filt <- uka_parsed %>%
          filter(Sample == condition) %>%
          uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff, cs = cs)
        # For random sampling, filter for cell
        uka_cell_all <- uka_parsed %>% filter(Sample == condition)

        # Prepare values to update in results csv
        vals <- list(
          n_kins = nrow(uka_filt),
          condition = condition,
          spec_cutoff = spec_cutoff,
          perc_cutoff = perc_cutoff
        )

        network_res <- compute_golden_network_scores(
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
          temp_files <- c(temp_files, network_res$temp_files)
          # Accumulate metrics
          all_metrics_df <- bind_rows(all_metrics_df, network_res$metrics_df)
          all_obs_metrics_df <- bind_rows(all_obs_metrics_df, network_res$obs_metrics)
        }
        # Update results df even if network_res is NULL (so condition/perc_cutoff are not skipped)
        results <- update_golden_results_row(results, vals)
      }
    }
  }

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
  writeLines(
    c(
      as.character(elapsed)
    ),
    log_time
  )

  logs <- c(logs, paste0("Run time: ", elapsed))

  logs <- c(logs, capture.output(print(results)))


  # Delete all temp files
  for (f in temp_files) {
    if (file.exists(f)) file.remove(f)
  }
  return(list(results = results, logs = logs))
}


# data cleaning functions
clean_uka_to_kinograte <- function(uka, cs = F) {
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

clean_tercen_columns <- function(df) {
  cols <- colnames(df)
  split <- str_split(cols, pattern = "\\.")
  cols <- sapply(split, tail, 1)
  colnames(df) <- cols
  return(df)
}

# Parameter capturing functions
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

  # Create parameter folder using the X variable name from symbol_names
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
  print(paste0("param folder: ", param_folder))
  if (!dir.exists(param_folder)) dir.create(param_folder, recursive = TRUE)
  params$respath <- param_folder
  save_params(params, respath = param_folder)
  return(params)
}

save_params <- function(params, respath) {
  symbol_names <- attr(params, "symbol_names")

  df <- tibble(
    parameter = names(symbol_names),
    value = unname(symbol_names)
  )
  df <- df %>%
    mutate(value = gsub('\\"', '"', value)) %>% # replace escaped quotes
    mutate(value = gsub("c\\(|\\)", "", value)) %>% # remove c( and )
    mutate(value = gsub('^"|"$', "", value)) %>% # remove leading/trailing quotes
    mutate(value = str_remove_all(value, "\\\\")) %>% # remove all backslashes
    mutate(value = str_remove_all(value, '"')) %>%
    filter(!parameter %in% c("uka_fam", "respath"))


  write_csv(df, paste0(respath, "/params.csv"))
  print(df)
}

# Visualization function
golden_results_to_barplot <- function(results_list, result_path) {
  
  if(!dir.exists(paste0(result_path, "/Score_plots"))){
    dir.create(paste0(result_path, "/Score_plots"), recursive = TRUE)
  }
  
  # If file paths are given, read them
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = F, id = "path") else x
  }) %>%
    bind_rows()
  results_df <- results_df %>% select(-`path...2`) %>% rename("path" = "path...1")
  
  browser()
  results_df_clean <- results_df %>%
    mutate(path = sub(".*Onlykin_fscore1.3\\/", "", path)) %>%
    mutate(path = sub("\\/results_notna.csv.*", "", path)) %>%
    separate_wider_delim(cols = "path", delim = "/", names = c("path", "params")) %>%
    filter(params == "thrs_0.01_ppi_networkv12_nperms30_ukaabs1_b1")

  
  p <- ggplot(results_df_clean, aes(x = condition, y = network_sig_score)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    labs(
      x = "Cell", y = "score (p-value)",
      title = unique(results_df_clean$params)
    ) +
    facet_wrap(~path, ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p
  
  
  ggsave(paste0(result_path, "/Score_plots/Barplot_all_scores", unique(results_df_clean$params), ".png"),
         p,
         w = 18, h = 20, unit = "cm"
  )
  
  
  # two metrics
  # results_df_long <- results_df_clean %>%
  #   pivot_longer(
  #     cols = c("score_sig_network", "score_sig_network_inv"),
  #     names_to = "Metric", values_to = "value"
  #   )
  # p <- ggplot(results_df_long, aes(x = condition, y = value, fill = Metric)) +
  #   geom_bar(stat = "identity", position = position_dodge()) +
  #   geom_hline(yintercept = 0.05, linetype = "dashed") +
  #   labs(
  #     x = "Cell", y = "score (p-value)",
  #     title = paste0("Percentile cutoff: ", perccutoff)
  #   ) +
  #   facet_wrap(~path, ncol = 2) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # p
  # 
  # 
  # ggsave("results/only_kinase_PTKSTK/Score_plots/Barplot_all_scores.png",
  #   p,
  #   w = 18, h = 20, unit = "cm"
  # )

  summary <- results_df %>%
    group_by(path) %>%
    summarize(
      n_sig_comparison = sum(score_sig_network < 0.07 | score_sig_network_inv < 0.07, na.rm = TRUE)
    )
  summary1 <- results_df %>%
    group_by(path, condition) %>%
    summarize(
      min_score = min(score_sig_network, score_sig_network_inv, na.rm = TRUE)
    )

  p1 <- ggplot(summary1, aes(x = condition, y = min_score)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    labs(x = "Cell", y = "network or network inv score (p-value)") +
    facet_wrap(~path, ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p1

  ggsave("results/only_kinase_PTKSTK/Score_plots/Barplot_min_score.png",
    p1,
    w = 16, h = 18, unit = "cm"
  )
}



golden_results_to_lineplot <- function(results_list) {
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = F, id = "path") else x
  }) %>%
    bind_rows() %>%
    mutate(path = sub(".*Score_results.*\\/(.*)\\/sens.*", "\\1", path)) %>%
    pivot_longer(
      cols = c(
        "score_sig_overlap", "score_sig_overlap_fam", "score_sig_network",
        "score_sig_core_network", "score_sig_core_network_inv"
      ),
      names_to = "Metric", values_to = "value"
    )
  # across cutoffs

  make_lineplot <- function(df, metric) {
    df_filt <- df %>% filter(Metric == metric)
    p <- ggplot(df_filt, aes(x = perc_cutoff, y = value, color = cell, linetype = Metric)) +
      geom_line() +
      facet_wrap(~path, ncol = 1) +
      labs(
        x = "Percentile cutoff",
        y = metric
      )
    p

    ggsave(paste0("results/Score_plots/Lineplot_", metric, ".png"), p, w = 12, h = 14, unit = "cm")
  }

  make_lineplot(results_df, "score_sig_network")
  make_lineplot(results_df, "score_sig_overlap")
  make_lineplot(results_df, "score_sig_overlap_fam")
  make_lineplot(results_df, "score_sig_core_network")
  make_lineplot(results_df, "score_sig_core_network_int")
}


create_combined_heatmap <- function(heatmap_data_list, w_combined, h, save_folder){
  # Find overlapping and union pathways
  all_pathways <- lapply(heatmap_data_list, function(x) rownames(x$matrix))
  common_pathways <- Reduce(intersect, all_pathways)
  union_pathways <- Reduce(union, all_pathways)
  
  # Calculate global color limits across all comparisons
  all_logfc_vals <- unlist(lapply(heatmap_data_list, function(x) {
    matrix_subset <- x$matrix[rownames(x$matrix) %in% common_pathways, , drop = FALSE]
    as.numeric(matrix_subset)
  }))
  all_logfc_vals <- all_logfc_vals[!is.na(all_logfc_vals) & is.finite(all_logfc_vals)]
  
  if (length(all_logfc_vals) > 0) {
    p25_neg <- quantile(all_logfc_vals[all_logfc_vals < 0], probs = 0.25, na.rm = TRUE)
    p75_pos <- quantile(all_logfc_vals[all_logfc_vals > 0], probs = 0.75, na.rm = TRUE)
    limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
    
    if (!is.finite(limit_val) || limit_val == 0) {
      limit_val <- max(abs(all_logfc_vals), na.rm = TRUE)
    }
    
    global_color_limits <- c(-limit_val, limit_val)
  } else {
    global_color_limits <- c(-1, 1)
  }
  
  col_fun <- colorRamp2(c(global_color_limits[1], 0, global_color_limits[2]), 
                        c("#0000CC", "#BDBDBD", "#C62828"))
  
  
  
  if (length(common_pathways) > 0) {
    
    # Create individual heatmaps for each comparison (common pathways only)
    individual_heatmaps <- list()
    
    for (comp_name in names(heatmap_data_list)) {
      # Get matrix subset for common pathways only
      comp_matrix <- heatmap_data_list[[comp_name]]$matrix
      comp_matrix_subset <- comp_matrix[rownames(comp_matrix) %in% common_pathways, , drop = FALSE]

      comp_matrix_subset <- comp_matrix_subset[
        !apply(comp_matrix_subset == 0, 1, all),
        !apply(comp_matrix_subset == 0, 2, all),
        drop = FALSE
      ]
      
      # Ensure same pathway order across all comparisons
      comp_matrix_subset <- comp_matrix_subset[common_pathways, , drop = FALSE]
      comp_matrix_subset[is.na(comp_matrix_subset)] <- 0
      
      # Create heatmap for this comparison
      ht_individual <- Heatmap(
        comp_matrix_subset,
        name = "LogFC",
        col = col_fun,
        cluster_rows = FALSE,  # Keep pathway order consistent
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = comp_name,
        column_title_gp = gpar(fontsize = 8, fontface = "bold"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(8, "cm"),
        rect_gp = gpar(col = "black", lwd = 0.5),
        show_heatmap_legend = FALSE  # Only show legend on first heatmap
      )
      
      individual_heatmaps[[comp_name]] <- ht_individual
    }
    
    # Show legend only on the first heatmap
    individual_heatmaps[[1]]@heatmap_param$show_heatmap_legend <- TRUE
    individual_heatmaps[[1]]@heatmap_param$heatmap_legend_param <- list(
      title = "LogFC", 
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
    
    # Combine heatmaps horizontally
    ht_combined <- Reduce(`+`, individual_heatmaps)
    
    # Calculate dimensions based on number of comparisons
    n_comparisons <- length(heatmap_data_list)
    n_pathways_combined <- length(common_pathways)
    
    # Adjust width based on number of comparisons
    if (!is.null(w_combined)){
      if (n_comparisons == 2) {
        w_combined <- 25
      } else if (n_comparisons == 3) {
        w_combined <- 35
      } else {
        w_combined <- 45
      }
    }
    
    
    
    # Save combined heatmap
    combined_filename <- paste0(save_folder, "/Kinase_Pathway_heatmap_overlapping_comparisons.png")
    png(combined_filename, 
        width = w_combined, 
        height = h, 
        units = "cm",
        res = 300)
    draw(ht_combined)
    dev.off()
    
    cat("Saved combined heatmap:", combined_filename, "\n")
    cat("Found", length(common_pathways), "overlapping pathways across", 
        length(heatmap_data_list), "comparisons\n")
    
  } else {
    cat("No overlapping pathways found across comparisons\n")
  }
}

# Function to create combined heatmap using union of pathways (for exactly 2 comparisons)
create_union_combined_heatmap <- function(heatmap_data_list, w_combined, h, save_folder) {
  # Only works with exactly 2 comparisons
  if (length(heatmap_data_list) != 2) {
    cat("Union combined heatmap only works with exactly 2 comparisons. Found", 
        length(heatmap_data_list), "comparisons.\n")
    return()
  }
  
  # Get pathway sets for each comparison
  comp_names <- names(heatmap_data_list)
  pathways_comp1 <- rownames(heatmap_data_list[[1]]$matrix)
  pathways_comp2 <- rownames(heatmap_data_list[[2]]$matrix)
  
  # Find common and unique pathways
  common_pathways <- intersect(pathways_comp1, pathways_comp2)
  unique_comp1 <- setdiff(pathways_comp1, pathways_comp2)
  unique_comp2 <- setdiff(pathways_comp2, pathways_comp1)
  
  # Create ordered pathway list: common, unique_comp1, unique_comp2
  ordered_pathways <- c(common_pathways, unique_comp1, unique_comp2)
  
  # Create pathway category factor for row splitting
  pathway_categories <- c(
    rep("1", length(common_pathways)),
    rep("2", length(unique_comp1)),
    rep("3", length(unique_comp2))
  )
  names(pathway_categories) <- ordered_pathways
  pathway_split <- factor(pathway_categories, levels = c("1", "2", "3"))
  
  cat("Found", length(common_pathways), "common pathways,", 
      length(unique_comp1), "unique to", comp_names[1], ",", 
      length(unique_comp2), "unique to", comp_names[2], "\n")
  
  # Calculate global color limits across all comparisons (including all pathways)
  all_logfc_vals <- unlist(lapply(heatmap_data_list, function(x) {
    as.numeric(x$matrix)
  }))
  all_logfc_vals <- all_logfc_vals[!is.na(all_logfc_vals) & is.finite(all_logfc_vals)]
  
  if (length(all_logfc_vals) > 0) {
    p25_neg <- quantile(all_logfc_vals[all_logfc_vals < 0], probs = 0.25, na.rm = TRUE)
    p75_pos <- quantile(all_logfc_vals[all_logfc_vals > 0], probs = 0.75, na.rm = TRUE)
    limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
    
    if (!is.finite(limit_val) || limit_val == 0) {
      limit_val <- max(abs(all_logfc_vals), na.rm = TRUE)
    }
    
    global_color_limits <- c(-limit_val, limit_val)
  } else {
    global_color_limits <- c(-1, 1)
  }
  
  col_fun <- colorRamp2(c(global_color_limits[1], 0, global_color_limits[2]), 
                        c("#0000CC", "#BDBDBD", "#C62828"))
  
  # Create matrices for each comparison with union of pathways
  individual_heatmaps <- list()
  
  for (i in 1:2) {
    comp_name <- comp_names[i]
    comp_matrix <- heatmap_data_list[[comp_name]]$matrix
    
    # Create a new matrix with all pathways, filled with zeros for missing pathways
    union_matrix <- matrix(0, nrow = length(ordered_pathways), ncol = ncol(comp_matrix),
                          dimnames = list(ordered_pathways, colnames(comp_matrix)))
    
    # Fill in the actual values for pathways that exist in this comparison
    existing_pathways <- intersect(ordered_pathways, rownames(comp_matrix))
    union_matrix[existing_pathways, ] <- comp_matrix[existing_pathways, , drop = FALSE]
    
    # Replace NA with 0
    union_matrix[is.na(union_matrix)] <- 0
    
    # Create heatmap for this comparison
    ht_individual <- Heatmap(
      union_matrix,
      name = "LogFC",
      col = col_fun,
      cluster_rows = FALSE,  # Keep pathway order consistent
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_split = pathway_split,  # Split rows by pathway category
      column_title = comp_name,
      column_title_gp = gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 6),
      row_names_max_width = unit(8, "cm"),
      rect_gp = gpar(col = "black", lwd = 0.5),
      show_heatmap_legend = FALSE  # Only show legend on first heatmap
    )
    
    individual_heatmaps[[comp_name]] <- ht_individual
  }
  
  # Show legend only on the first heatmap
  individual_heatmaps[[1]]@heatmap_param$show_heatmap_legend <- TRUE
  individual_heatmaps[[1]]@heatmap_param$heatmap_legend_param <- list(
    title = "LogFC", 
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
  
  # Combine heatmaps horizontally
  ht_combined <- Reduce(`+`, individual_heatmaps)
  
  # Set width for 2 comparisons
  if (is.null(w_combined)) {
    w_combined <- 35  # Default for union heatmap with 2 comparisons
  }
  
  # Save combined union heatmap
  union_filename <- paste0(save_folder, "/Kinase_Pathway_heatmap_union_comparisons.png")
  png(union_filename, 
      width = w_combined, 
      height = h, 
      units = "cm",
      res = 300)
  draw(ht_combined)
  dev.off()
  
  cat("Saved union combined heatmap:", union_filename, "\n")
}

# Function to create kinase-pathway heatmaps from network analysis results
plot_kinase_pathway_heatmaps <- function(result_folder, save_folder = NULL, w = 20, h = 14, w_combined = 25) {
  
  if (is.null(save_folder)) {
    save_folder <- file.path(result_folder, "kinase_pathway_heatmaps")
  }
  
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)
  }
  
  # Find all nodes files with pathways
  nodes_pw_files <- list.files(result_folder, 
                              pattern = "nodes_.*_with_pathways\\.csv$", 
                              full.names = TRUE, recursive = TRUE)
  
  if (length(nodes_pw_files) == 0) {
    stop("No nodes_*_with_pathways.csv files found in the specified folder")
  }
  
  # Extract comparison names from filenames
  extract_comparison_name <- function(filepath) {
    basename_file <- basename(filepath)
    # Remove nodes_ prefix and _with_pathways.csv suffix
    comparison <- gsub("^nodes_(.*)_with_pathways\\.csv$", "\\1", basename_file)
    # Remove _p## suffix (percentile cutoff)
    comparison <- gsub("_p[0-9]+$", "", comparison)
    return(comparison)
  }
  
  comparison_names <- sapply(nodes_pw_files, extract_comparison_name)
  names(nodes_pw_files) <- comparison_names
  
  # Process each comparison
  heatmap_data_list <- list()
  
  for (i in seq_along(nodes_pw_files)) {
    comparison <- names(nodes_pw_files)[i]
    pw_file <- nodes_pw_files[i]
    
    cat("Processing comparison:", comparison, "\n")
    
    # Find corresponding nodes file (without pathways)
    nodes_file_pattern <- gsub("_with_pathways\\.csv$", "\\.csv", basename(pw_file))
    comparison_p <- gsub("nodes_", "", nodes_file_pattern)
    nodes_file <- list.files(dirname(pw_file), 
                            pattern = paste0("^nodes_", comparison_p, "$"), 
                            full.names = TRUE)[1]
    
    if (is.na(nodes_file) || !file.exists(nodes_file)) {
      warning("Could not find corresponding nodes file for: ", comparison)
      next
    }

    # Read the files
    tryCatch({
      nodes_pw <- read_csv(pw_file, show_col_types = FALSE)
      nodes <- read_csv(nodes_file, show_col_types = FALSE)
      
      # Separate longer delim for pathways and join with nodes data
      nodes_pw_expanded <- nodes_pw %>%
        separate_longer_delim(pathway, delim = ",") %>%
        filter(!is.na(pathway), pathway != "") %>%
        left_join(nodes %>% select(Protein, type, LogFC), 
                 by = "Protein") %>%
        filter(type %in% c("Kinase", "Artificial")) %>%
        select(Protein, pathway, LogFC) %>%
        distinct()
      
      if (nrow(nodes_pw_expanded) == 0) {
        warning("No kinase data found for comparison: ", comparison)
        next
      }
      
      # Create matrix for heatmap
      heatmap_matrix <- nodes_pw_expanded %>%
        pivot_wider(names_from = Protein, values_from = LogFC, 
                   values_fn = mean, values_fill = NA) %>%
        column_to_rownames("pathway") %>%
        as.matrix()
      
      # Remove rows and columns with all NA
      heatmap_matrix <- heatmap_matrix[!apply(is.na(heatmap_matrix), 1, all), 
                                     !apply(is.na(heatmap_matrix), 2, all), drop = FALSE]
      heatmap_matrix[is.na(heatmap_matrix)] <- 0
      
      if (ncol(heatmap_matrix) == 0 || nrow(heatmap_matrix) == 0) {
        warning("Empty matrix for comparison: ", comparison)
        next
      }
      
      # Store for later combined analysis
      heatmap_data_list[[comparison]] <- list(
        matrix = heatmap_matrix,
        data = nodes_pw_expanded
      )
      
      # Apply the same color logic as visualize_network_pg
      logfc_vals <- as.numeric(heatmap_matrix)
      logfc_vals <- logfc_vals[!is.na(logfc_vals) & is.finite(logfc_vals)]
      
      if (length(logfc_vals) > 0) {
        p25_neg <- quantile(logfc_vals[logfc_vals < 0], probs = 0.25, na.rm = TRUE)
        p75_pos <- quantile(logfc_vals[logfc_vals > 0], probs = 0.75, na.rm = TRUE)
        limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
        
        if (!is.finite(limit_val) || limit_val == 0) {
          limit_val <- max(abs(logfc_vals), na.rm = TRUE)
        }
        
        color_limits <- c(-limit_val, limit_val)
      } else {
        color_limits <- c(-1, 1)
      }
      
      # Color function (blue -> gray -> red)
      col_fun <- colorRamp2(c(color_limits[1], 0, color_limits[2]), 
                           c("#0000CC", "#BDBDBD", "#C62828"))
      
      # Calculate reasonable dimensions
      n_kinases <- ncol(heatmap_matrix)
      n_pathways <- nrow(heatmap_matrix)
      
      # Base width and height with reasonable limits
      if (!is.null(w)){
        w <- max(20, max(8, n_kinases * 0.4))  # Max 1600px wide
      }
      if (!is.null(h)){
        h <- max(14, max(4, n_pathways * 0.25 + 2))  # Extra space for row names
      }
      
      
      # Create heatmap
      ht <- Heatmap(
        heatmap_matrix,
        name = "LogFC",
        col = col_fun,
        cluster_rows = TRUE,  # cluster pathways
        cluster_columns = TRUE,  # Cluster kinases
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = paste("Kinase-Pathway Heatmap:", comparison),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_names_gp = gpar(fontsize = 6),  # Smaller font for long pathway names
        column_names_gp = gpar(fontsize = 8),
        row_names_max_width = unit(8, "cm"),  # More space for row names
        rect_gp = gpar(col = "black", lwd = 0.5),
        heatmap_legend_param = list(title = "LogFC", 
                                   title_gp = gpar(fontsize = 10),
                                   labels_gp = gpar(fontsize = 8))
      )
      
      # Save individual heatmap
      filename <- paste0(save_folder, "/Kinase_Pathway_heatmap_", 
                        gsub("[^A-Za-z0-9_-]", "_", comparison), ".png")
      png(filename, 
          width = w, 
          height = h,
          units = "cm",
          res = 300)
      draw(ht)
      dev.off()
      
      cat("Saved heatmap:", filename, "\n")
      
    }, error = function(e) {
      warning("Error processing ", comparison, ": ", e$message)
    })
  }
  
  # Create combined heatmap for overlapping pathways
  if (length(heatmap_data_list) >= 2 && length(heatmap_data_list) <= 4) {
    
    create_combined_heatmap(heatmap_data_list = heatmap_data_list, w_combined = w_combined, h = h, 
                            save_folder = save_folder)
    
    # Create union combined heatmap for exactly 2 comparisons
    if (length(heatmap_data_list) == 2) {
      create_union_combined_heatmap(heatmap_data_list = heatmap_data_list, w_combined = w_combined, h = h, 
                                   save_folder = save_folder)
    }
    
  } else if (length(heatmap_data_list) > 4) {
    cat("Too many comparisons (", length(heatmap_data_list), ") for combined heatmap. Maximum 4 allowed. Skipping combined heatmap.\n")
  }
  
  cat("Heatmap generation completed. Files saved in:", save_folder, "\n")
  invisible(heatmap_data_list)
}
