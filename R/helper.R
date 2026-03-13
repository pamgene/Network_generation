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


# Helper to compute overlap and family-wise overlap scores and manage temp files
# Helper to compute network scores and manage temp files

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
# data cleaning functions
clean_uka_to_kinograte1 <- function(uka, control, spec_cutoff = 0, del_cell = NULL, cs = F){

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
  sens_name <- if("sens" %in% names(symbol_names)) symbol_names[["sens"]] else NULL
  # best_drug_target_name <- symbol_names[["best_drug_per_target"]]
  ppi_network_name <- symbol_names[["ppi_network"]]
  
  # Build folder structure
  param_folder <- file.path(params$respath, 
        if(!is.null(sens_name)) uka_name else "",
        paste0(
          if(!is.null(sens_name)) paste0(sens_name, "_") else "",
          "thrs_", paste0(params$spec_cutoffs, collapse = "-"),
          if(is.null(sens_name)) paste0("_perc", params$perc_cutoff) else "",
          "_", ppi_network_name,
          if("nperms_network" %in% names(params)) paste0("_nperms", params$nperms_network) else "",
          if("relative_to" %in% names(params)) paste0("_rel2", params$relative_to) else "",
          "_ukaabs", as.integer(params$rank_uka_abs), "_b", params$b,
          if("cs" %in% names(params)) paste0("_cs", as.integer(params$cs)) else "",
          if(!is.null(params$art_nodes)) "_art" else ""
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



