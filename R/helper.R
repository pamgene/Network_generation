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
        aes(xintercept = value), color = "red", linetype = "dashed") +
      facet_wrap(~cell, scales = "free_y") +
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


# Generalized permutation scoring function
# Permute UKA, keep sens constant
pOverlapScore = function(score_obs, uka_cell_all, sens_filt, cell, nsample, nPerms = DFT.PERMS, 
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

pNetworkScore = function(score_obs, score_obs_core, score_obs_core_inv, uka_cell_all, sens_filt, cell, nsample, nPerms = DFT.PERMS, 
                         perc_cutoff, spec_cutoff, respath, condition,
                         rank_uka_abs, ppi_network, relative_to, b = b) {  
  
  score_perm <- rep(NA, nPerms)
  score_perm_core <- rep(NA, nPerms)
  score_perm_core_inv <- rep(NA, nPerms)
  # Replace for loop with lapply
  perm_results <- lapply(1:nPerms, function(i) {
    uka_cell_all$uniprotname <- sample(uka_cell_all$uniprotname)
    uka_filt_random <- uka_cell_all %>% uka_top(rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
    print(paste0("Cell: ", cell, ", Network permutation: ", i, "..."))
    res <- make_network_and_stats(uka = uka_filt_random, sens = sens_filt, perc_cutoff = perc_cutoff, 
                                  spec_cutoff = spec_cutoff, res.path = respath, condition = cell, write = FALSE, ppi_network = ppi_network,
                                  relative_to = relative_to, b = b)
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



uka_top <- function(uka, rank_uka_abs = T, perc_cutoff){
  
  if (rank_uka_abs) {
    uka_rank <- percentile_rank_fast(uka, uniprotname, LogFC) # kinograte::percentile_rank of abs LogFC of kinases
  } else {
    uka_rank <- percentile_rank_pg_noabs(uka, symbol = uniprotname, metric = LogFC, rank_lowest_highest = F)
  }
  uka_rank %>% dplyr::filter(score >= perc_cutoff) %>% # & score >= perc_cutoff
      dplyr::ungroup() %>% 
      dplyr::mutate(type = 'Kinase') %>% 
      dplyr::select(name, prize = score, type, LogFC)
}

sens_top <- function(sens, perc_cutoff, balance = T){
  if (balance){
    perc_cutoff <- perc_cutoff - 0.2
  } 
  sens_rank <- percentile_rank_pg_noabs(sens, uniprotname, LogFC, rank_lowest_highest = T)
  sens_rank %>% dplyr::filter(score >= perc_cutoff) %>% dplyr::ungroup() %>% 
    dplyr::mutate(type = "Sensitivity") %>% 
    dplyr::select(name, prize = score, type, LogFC)
}

overlap_uka_sens <- function(uka, sens){
  overlap <- length(intersect(uka, sens))
  # return relative overlap
  overlap/max(length(sens), 1)
}


# Helper to create empty results data frame for make_golden_score
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

# Helper to update a row in the results data frame
update_golden_results_row <- function(results, cell, perc_cutoff, vals) {
  idx <- which(results$perc_cutoff == perc_cutoff & results$cell == cell)
  for (nm in names(vals) ) {
    results[idx, nm] <- vals[[nm]]
  }
  results
}

# Helper to compute overlap and family-wise overlap scores and manage temp files
compute_golden_overlap_scores <- function(uka_filt, sens_filt, uka_cell_all, uka_fam, cell, perc_cutoff, nperms_overlap, respath, spec_cutoff, rank_uka_abs) {
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
    family_mode = T
  )
  vals$score_sig_overlap_fam <- overlap_sig_fam

  # Write temp file for score_sig_overlap_fam
  temp_overlap_fam_file <- file.path(respath, paste0("temp_overlap_fam_", cell, "_", perc_cutoff, ".txt"))
  writeLines(as.character(overlap_sig_fam), temp_overlap_fam_file)
  temp_files <- c(temp_files, temp_overlap_fam_file)

  list(vals = vals, temp_files = temp_files)
}

# Helper to compute network scores and manage temp files

# New version: returns metrics for this cell/perc_cutoff only; accumulation is handled outside
compute_golden_network_scores <- function(uka_filt, sens_filt, uka_cell_all, cell, 
                                          perc_cutoff, spec_cutoff, respath, rank_uka_abs, 
                                          ppi_network, nperms_network, relative_to, b) {
  vals <- list()
  # observed network score and metrics
  
  obs_res <- make_network_and_stats(
    uka = uka_filt, sens = sens_filt,
    perc_cutoff = perc_cutoff, spec_cutoff = spec_cutoff,
    res.path = respath, condition = cell, write = F,
    ppi_network = ppi_network, relative_to = relative_to, b = b
  )
  obs_network <- obs_res$rel_med_path_KT
  
  obs_network_core <- obs_res$rel_med_path_KT_core
  obs_network_core_inv <- obs_res$rel_med_path_KT_core_inv
  obs_metrics_df_long <- as.data.frame(obs_res) %>% 
    pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>% 
    mutate(cell = cell, perc_cutoff = perc_cutoff) %>%
    select(cell, perc_cutoff, metric, value)
  vals$obs_network <- obs_network
  vals$obs_network_core <- obs_network_core
  vals$obs_network_core_inv <- obs_network_core_inv

  # collect metrics from permutations
  perm_res <- pNetworkScore(
    score_obs = obs_network,
    score_obs_core = obs_network_core,
    score_obs_core_inv = obs_network_core_inv,
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
}

make_golden_score <- function(uka, sens, del_cells = NULL, 
                              control, spec_cutoff, zscore = F, best_drug_per_target = NULL, respath, score_overlap = T, 
                              uka_fam, score_network = T, perc_cutoffs, nperms_overlap = 500, nperms_network = 50,
                              rank_uka_abs = T, balance = F, ppi_network = ppi_networkv12,
                              relative_to = "n_nodes", b, cs) {

  # input data uka: uka output for all cells
  # input data sens: sensitivity for all cells
  # del_cells: cells that have very few FDR significant nodes, should be deleted
  # rank_lowest_highest: If TRUE, the lowest values get the highest rank (i.e., most sensitive/lowest LogFC is ranked highest).
  # balance: if T, perc_cutoff for sens is 0.2 smaller to balance dataset sizes out
  # Parse input data, filter for comparisons vs control
  sens_parsed <- clean_sens_to_kinograte(sens, control = control, zscore = zscore,
                                           best_drug_per_target = best_drug_per_target)
  
  uka_parsed <- clean_uka_to_kinograte(uka, spec_cutoff = spec_cutoff, control = "RL", cs = cs)
  # Find common cell lines in uka and sens
  common_cells <- intersect(unique(uka_parsed$cell_line), unique(sens_parsed$cell_line))
  if (!is.null(del_cells)){
    common_cells <- common_cells[!common_cells %in% del_cells]   
  }
  # Create empty results df
  results <- make_golden_results_df(common_cells, perc_cutoffs)
  
  # For plotting: accumulate metrics for each perc_cutoff and cell
  metric_names <- c("density","modularity", "assortativity", 
                    "rel_med_path_KT", "rel_med_path_KT_core",
                    "rel_med_path_all_nodes", 'rel_med_path_all_nodes_core', 
                    "rel_med_path_KT_vs_all_nodes", 
                    "rel_med_path_KT_vs_all_nodes_core",
                    "rel_med_path_KT_core_inv", "rel_med_path_all_core_inv", "rel_med_path_KT_vs_all_nodes_core_inv",
                    "avg_sens_clustering")
  # Initialize metrics df using helper
  all_metrics_df <- data.frame()
  all_obs_metrics_df <- data.frame()
  

  # Sequential loop over percentile cutoffs
  loop_start <- Sys.time()
  logs <- character() # Collect log messages
  for (perc_cutoff in perc_cutoffs) {
    # Sequential loop over cells for this percentile cutoff
    temp_files <- c()
    for (cell in common_cells) {
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
        network_res <- compute_golden_network_scores(
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
        
        # Accumulate metrics
        all_metrics_df <- bind_rows(all_metrics_df, network_res$metrics_df)
        all_obs_metrics_df <- bind_rows(all_obs_metrics_df, network_res$obs_metrics)
      } else {
        vals$obs_network = NA
        vals$obs_network_core = NA
        vals$score_sig_network = NA
        vals$score_sig_core_network = NA
      }

      # Update results df
      results <- update_golden_results_row(results, cell, perc_cutoff, vals)
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
  print(elapsed)

  logs <- c(logs, paste0("Run time: ", elapsed))

  logs <- c(logs, capture.output(print(results)))
  

  # Delete all temp files
  for (f in temp_files) {
    if (file.exists(f)) file.remove(f)
  }
  
  return(list(results = results, logs = logs))
}

make_golden_score_p <- function(uka, sens, del_cells = NULL, 
                              control, spec_cutoff, zscore = F, best_drug_per_target = NULL, respath, score_overlap = T, 
                              uka_fam, score_network = T, perc_cutoffs, nperms_overlap = 500, nperms_network = 50,
                              rank_uka_abs = T, balance = F, ppi_network = ppi_networkv12,
                              relative_to = "n_nodes", b, cs) {

  # input data uka: uka output for all cells
  # input data sens: sensitivity for all cells
  # del_cells: cells that have very few FDR significant nodes, should be deleted
  # rank_lowest_highest: If TRUE, the lowest values get the highest rank (i.e., most sensitive/lowest LogFC is ranked highest).
  # balance: if T, perc_cutoff for sens is 0.2 smaller to balance dataset sizes out
  # Parse input data, filter for comparisons vs control
  sens_parsed <- clean_sens_to_kinograte(sens, control = control, zscore = zscore,
                                           best_drug_per_target = best_drug_per_target)
  uka_parsed <- clean_uka_to_kinograte(uka, spec_cutoff = spec_cutoff, control = "RL", cs = cs)
  # Find common cell lines in uka and sens
  common_cells <- intersect(unique(uka_parsed$cell_line), unique(sens_parsed$cell_line))
  if (!is.null(del_cells)){
    common_cells <- common_cells[!common_cells %in% del_cells]   
  }
  # Create empty results df
  results <- make_golden_results_df(common_cells, perc_cutoffs)
  
  # For plotting: accumulate metrics for each perc_cutoff and cell
  metric_names <- c("density","modularity", "assortativity", 
                    "rel_med_path_KT", "rel_med_path_KT_core",
                    "rel_med_path_all_nodes", 'rel_med_path_all_nodes_core', 
                    "rel_med_path_KT_vs_all_nodes", 
                    "rel_med_path_KT_vs_all_nodes_core",
                    "rel_med_path_KT_core_inv", "rel_med_path_all_core_inv", "rel_med_path_KT_vs_all_nodes_core_inv",
                    "avg_sens_clustering")
  # Initialize metrics df using helper
  all_metrics_df <- data.frame()
  all_obs_metrics_df <- data.frame()
  

  # Sequential loop over percentile cutoffs
  loop_start <- Sys.time()
  logs <- character() # Collect log messages
    for (perc_cutoff in perc_cutoffs) {
      temp_files <- c()
      # Parallelize over common_cells using foreach
      if (requireNamespace("foreach", quietly = TRUE)) {
        cell_results <- foreach::foreach(cell = common_cells, 
                                         .packages = c("dplyr", "readr", "tidyr", "purrr", "rlang", "PCSF", "igraph")) %dopar% {
          uka_filt <- uka_parsed %>% filter(cell_line == cell) %>% uka_top(rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
          sens_filt <- sens_parsed %>% filter(cell_line == cell) %>% sens_top(perc_cutoff, balance)
          uka_cell_all <- uka_parsed %>% filter(cell_line == cell)
          vals <- list(
            n_targets = nrow(sens_filt),
            n_kins = nrow(uka_filt),
            max_sens_value = max(sens_filt$LogFC),
            cell = cell,
            perc_cutoff = perc_cutoff
          )
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
            temp_files <- c(overlap_res$temp_files)
          } else {
            vals$obs_overlap = NA
            vals$score_sig_overlap = NA
            vals$obs_overlap_fam = NA
            vals$score_sig_overlap_fam = NA
          }
          # --- Golden score 2: Network ---
          metrics_df <- NULL
          obs_metrics <- NULL
          if (score_network) {
            network_res <- compute_golden_network_scores(
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
          list(cell = cell, perc_cutoff = perc_cutoff, vals = vals, temp_files = temp_files, metrics_df = metrics_df, obs_metrics = obs_metrics)
        }
        # Process parallel results
        for (result in cell_results) {
          results <- update_golden_results_row(results, result$cell, result$perc_cutoff, result$vals)
          temp_files <- c(temp_files, result$temp_files)
          if (!is.null(result$metrics_df)) all_metrics_df <- bind_rows(all_metrics_df, result$metrics_df)
          if (!is.null(result$obs_metrics)) all_obs_metrics_df <- bind_rows(all_obs_metrics_df, result$obs_metrics)
        }
      } else {
        stop("The 'foreach' package is required for parallel execution but is not available.")
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
  print(elapsed)

  logs <- c(logs, paste0("Run time: ", elapsed))

  logs <- c(logs, capture.output(print(results)))
  

  # Delete all temp files
  for (f in temp_files) {
    if (file.exists(f)) file.remove(f)
  }
  
  return(list(results = results, logs = logs))
}



# data cleaning functions
clean_uka_to_kinograte <- function(uka, control, spec_cutoff = 0, del_cell = NULL, cs = F){

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

clean_sens_to_kinograte <- function(sens, control, zscore = F, del_cell = NULL, best_drug_per_target){
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
  if (!is.null(control) & zscore == F){
    control_df <- sens_filt %>% filter(cell_line == control)
    sens_clean <- sens_filt %>%
      left_join(control_df, by = c("uniprotname"), suffix = c("", ".control")) %>%
      mutate(
        LogFC = LN_IC50 - LN_IC50.control
      ) %>%
      filter(cell_line != control, 
             !is.na(LogFC)) %>%
      dplyr::select(cell_line, uniprotname, LogFC)
  } else if (is.null(control) & zscore == F){
    # use median sensitivity data based on LogFC
    sens_clean <- sens_filt
  } else if (is.null(control) & zscore == T){
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
  sens_name <- symbol_names[["sens"]]
  # best_drug_target_name <- symbol_names[["best_drug_per_target"]]
  ppi_network_name <- symbol_names[["ppi_network"]]
  param_folder <- file.path(params$respath, uka_name, 
        paste0(
          sens_name, "_thrs_", paste0(params$perc_cutoffs, collapse = "-"),
          "_", ppi_network_name,
          '_nperms', params$nperms_network,
          "_blnc", as.integer(params$balance), "_rel2", params$relative_to,
          "_ukaabs", as.integer(params$rank_uka_abs), "_b", params$b
  ))
  print(paste0("param folder: ", param_folder))
  if (!dir.exists(param_folder)) dir.create(param_folder, recursive = TRUE)
  params$respath <- param_folder
  save_params(params, respath = param_folder)
  return(params)
}

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

# Visualization function
golden_results_to_barplot_singlemetric <- function(results_list, perccutoff = 0.7) {
  # If file paths are given, read them
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = F, id = "path") else x
  }) %>% bind_rows() %>%
    mutate(path = sub('.*_results\\/(.*)\\/sens.*', '\\1', path)) %>%
    filter(perc_cutoff == perccutoff) %>% 
    pivot_longer(
      cols = c("score_sig_overlap", "score_sig_overlap_fam", "score_sig_network", 
               "score_sig_core_network", "score_sig_core_network_inv"), 
      names_to = "Metric", values_to = "value")
  
  browser()
  
  # filter (always different filter)
  results_df <- results_df %>% filter(str_detect(path, "df_simple"))
  
  
  make_barplot <- function(df_long, metric){
    df_long_filt <- df_long %>% filter(Metric == metric)
    p <- ggplot(df_long_filt, aes(x = cell, y = value, fill = path)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      labs(x = "Cell line",
           y = "significance scores") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~path, ncol = 1) +
      ggtitle(paste0("Percentile cutoff: ", perccutoff, ", metric: ", metric))
    p
    ggsave(paste0("results/Score_plots/Barplot_all_scores_p", metric, "_", perccutoff, ".png"), 
           p, w = 13, h = 5, unit = "cm")
    
  }
  browser()
  make_barplot(results_df, "score_sig_network")
  make_barplot(results_df, "score_sig_core_network")
  make_barplot(results_df, "score_sig_core_network_inv")
  
  
  
  plot(results_df$n_kins, results_df$score_sig_overlap, xlab= "Number of kinases", ylab = "Overlap score")
  plot(results_df$n_kins, results_df$score_sig_overlap_fam, xlab= "Number of kinases", ylab = "Overlap family score")
  plot(results_df$n_kins, results_df$score_sig_network, xlab= "Number of kinases", ylab = "Network score")
  plot(results_df$path, results_df$score_sig_overlap, xlab= "dbs", ylab = "Overlap score")
  plot(results_df$path, results_df$score_sig_overlap_fam, xlab= "dbs", ylab = "Overlap family score")
  plot(results_df$path, results_df$score_sig_network, xlab= "dbs", ylab = "Network score")
  
}


golden_results_to_barplot <- function(results_list, perccutoff = 0.7) {
  # If file paths are given, read them
  browser()
  # golden cols
  my_cols <- c("score_sig_overlap", "score_sig_overlap_fam", "score_sig_network", 
                   "score_sig_core_network", "score_sig_core_network_inv")
  # other project cols
  my_cols <- c('score_sig_network', 'score_sig_network_inv')
  
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = F, id = "path") else x
  }) %>% bind_rows()%>%
    pivot_longer(
      cols = all_of(my_cols), 
      names_to = "Metric", values_to = "value")
  
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = F, id = "path") else x
  }) %>% bind_rows() %>%
    mutate(path = sub(".*gds_uka_cs_ivPNKL_log10_ab\\/(.*)\\/results.*", '\\1', path)) %>% 
    # ".*_results\\/(.*)\\/sens.*"
    # filter(perc_cutoff == perccutoff) %>% 
    # other projects
    select(-obs_network, -obs_network_inv) %>%
    # mutate(condition1 = paste0(condition, " (# kins: ", n_kins, ")")) %>%
    pivot_longer(
        cols = all_of(my_cols), 
      names_to = "Metric", values_to = "value")
  browser()
  # filter (always different filter)
  results_df <- results_df %>% filter(str_detect(path, "uka_iv"))
 
  results_df_network <- results_df %>% filter(Metric %in% c("score_sig_network", 
                                                            "score_sig_core_network", "score_sig_core_network_inv")) %>%
    distinct(perc_cutoff, cell, Metric, value)
  p <- ggplot(results_df_network, aes(x = cell, y = value, fill = Metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    labs(x = "condition", y = "score (p-value)",
         title = "Top 30%") +
    # facet_wrap(~path, ncol = 1)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p
  browser()
  ggsave(paste0("results/Score_plots/Scores_barplot_scoredf_simple_ivPNKL_ab_recognition", ".png"), 
         p, w = 13, h = 8, unit = "cm")
  
}


golden_results_to_lineplot <- function(results_list){
  results_df <- lapply(results_list, function(x) {
    if (is.character(x)) readr::read_csv(x, show_col_types = F, id = "path") else x
  }) %>% bind_rows() %>%
    mutate(path = sub('.*Score_results.*\\/(.*)\\/sens.*', '\\1', path)) %>%
    pivot_longer(
      cols = c("score_sig_overlap", "score_sig_overlap_fam", "score_sig_network", 
               "score_sig_core_network", "score_sig_core_network_inv"), 
      names_to = "Metric", values_to = "value")
  # across cutoffs
  
  make_lineplot <- function(df, metric){
    df_filt <- df %>% filter(Metric == metric)
    p <- ggplot(df_filt, aes(x = perc_cutoff, y = value, color = cell, linetype = Metric)) +
      geom_line() + 
      facet_wrap(~path, ncol = 1)+
      labs(x = "Percentile cutoff",
           y = metric)
    p
    
    ggsave(paste0("results/Score_plots/Lineplot_", metric, ".png"), p, w = 12, h = 14, unit = "cm")
    
  }
  
  make_lineplot(results_df, "score_sig_network")
  make_lineplot(results_df, "score_sig_overlap")
  make_lineplot(results_df, "score_sig_overlap_fam")
  make_lineplot(results_df, "score_sig_core_network")
  make_lineplot(results_df, "score_sig_core_network_int")
  
}


# Reconstruct partial results from temp files in a folder
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

# Find missing cell/perc_cutoff combinations (where any temp file is missing)
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

