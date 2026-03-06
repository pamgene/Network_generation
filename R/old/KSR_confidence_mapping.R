library(pgUpstream)
library(pgSupertest)
library(pgscales)
library(pgFCS)

# Parameters for KSR-UKA

make_K <- function(scores_df, min_score, min_score_p = F, db, rmDualSpecKinases = T, minSeqHom = 0.9, 
                          minKxsScore = 300, maxRank = 12, kinaseFamily){
  # Computes matrix K (KSR score between peptides (rows) and kinases (columns)), by processing UKA DB to normalized confidence scores.
  # Scores under the min_score are converted to 0 to select specific KSRs. 
  # min_score_p = True, means that min_score is a percentile on non-0 data.E.g. 0.2 p means that the bottom 20 % of confidence values are converted to 0.
  
  db_in = db %>% 
    mutate(ID = as.factor(ID)) %>%
    filter(PepProtein_SeqSimilarity >= minSeqHom, ) %>% 
    filter(Kinase_PKinase_PredictorVersion2Score >= minKxsScore) %>% 
    filter(Kinase_Rank <= maxRank) %>% 
    filter(family == kinaseFamily)
  
  if (rmDualSpecKinases){
    db_in <- if (kinaseFamily == "PTK"){
      db_in %>% filter(Kinase_group == "TYR")
    } else {
      db_in %>% filter(Kinase_group != "TYR")
    }
  }

  nKinases = db_in %>% 
    pull(Kinase_Name) %>% 
    length()
  
  scores_df <- scores_df %>% 
    mutate(cscore = ifelse(cscore < 1/nKinases, 1/nKinases, cscore))
  
  

  # get UKA db with normalized weights:
  # dbw = db_in %>%
  #   add_scores(scores_df) %>%
  #   cscore2w(na.wt = 1/nKinases) %>%
  #   arrange(-wn)
  
  # instead of add_scores
  addDbWeights = function(db, scores_df){
    db = db %>%
      left_join(scores_df, by = c("Database", "Kinase_Rank"))
    
    if (any(is.na(db$cscore))){
      stop("Missing or incorrect database weights!")
    }
    return(db)
  }
  
  combinedScores = function(db, minscore, outOfGroupScore = 0){
    dbc = db %>%
      mutate(Kinase_Name = Kinase_Name %>% as.factor,
             ID = ID %>% as.factor) %>%
      group_by(Kinase_Name, ID) %>%
      summarise(s = 1-prod(1-cscore)) #, outOfGroup = all(outOfGroup
    
    
    db.all = expand.grid(levels(droplevels(dbc$Kinase_Name)), levels(droplevels(dbc$ID)) )
    colnames(db.all) = c("Kinase_Name", "ID")
    db.all %>%
      left_join(dbc, by = c("Kinase_Name", "ID")) %>%
        mutate(s = case_when(
          s < minscore ~ minscore,
          s >= minscore  ~ s,
          is.na(s) ~ 0)) %>% # 0 where there is no KSR
        select(Kinase_Name, ID, s)

  }
  
  normalizeScores = function(db){
    sumdf = db %>% group_by(ID) %>% summarise(sumsc = sum(s))
    db %>%
      left_join(sumdf, by = "ID") %>%
      group_by(ID) %>%
      mutate(sc.nor = s/sumsc)
      # mutate(sc.nor = (s/sumsc)*100)
  }

  
  dbw = db_in %>%
    addDbWeights(scores_df) %>%
    arrange(ID, Kinase_Name, Database, Kinase_Rank) %>%
    distinct(ID, Kinase_Name, Database, .keep_all = T
             ) %>% # ID-Kinase-Db data can have multiple PepProteins in between and multiple Psite-residues. One datapoint / ID-Kinase-Db is kept with the lowest rank
    combinedScores(minscore = min_score) %>%
    normalizeScores()

  temp <- dbw %>% 
    reshape2::acast(ID ~ Kinase_Name, value.var = "sc.nor",fill = 0)
  # print(paste0("K dimensions before filtering min score: ", dim(temp)[1], ", ", dim(temp)[2]))
  
  hist(temp, 1000, xlim = c(0, 0.1), main = "normalized confidence scores")
  
  
  if (min_score_p){
    # if min_score_p is True, that is min_score is a percentile on non-0 data:
    temp1 <- dbw %>% filter(sc.nor != 0)
    min_score_q <- quantile(temp1$sc.nor, min_score, na.rm = TRUE)
    # print(paste0("Min score at ", min_score, " percentile of non-0 data: ", min_score_q))
  } else {
    min_score_q <- min_score
  }
  
  
  K = dbw %>% 
    filter(sc.nor > min_score_q) %>% 
    reshape2::acast(ID ~ Kinase_Name, value.var = "sc.nor",fill = 0)
  
  # print(paste0("K dimensions after filtering min score: ", dim(K)[1], ", ", dim(K)[2]))
  
  return(K)
  
}


make_X <- function(S, K){
  # Computes matrix X of cell lines (rows) and kinases (columns).
  # the signal for kinase1 is taken as the KSR weighted average of the peptide signals that are included for kinase1.
  # [sample x spot] * [spot x kinase] = [sample x kinase]
  # [VSN] * [KSR] = Kinase data
  
  
  ov = intersect(colnames(S), rownames(K))
  # print(paste0("Number of IDs in both H and K: ", length(ov)))
  
  S = S[, colnames(S) %in% ov]
  K = K[rownames(K) %in% ov,]
  
  X = S %*% K
  # print(paste0("Number of possible kinases: ", dim(X)[2]))
  X <- X[, colSums(X != 0) > 0]
  # print(paste0("Number of kinases for which at least 1 sample is non-0: ", dim(X)[2]))
  
  s0 = apply(X,2, sd)
  m0 = apply(X,2, mean)
  
  # plot(m0, s0) # X log normal?
  
  X = X[, s0>0] # remove kinases with sd = 0
  
  
  return(X)
}

X_to_uka <- function(X, control = NULL){
  if (is.null(control)) {
    as.data.frame(X) %>%
      rownames_to_column("sample") %>%
      pivot_longer(cols = !sample, names_to = "Kinase_Name", values_to = "norm_conf_score")
  } else {
    if (!(control %in% rownames(X))) {
      stop("Control cell not found in rownames of X")
    }
    X_control <- X[control, , drop = FALSE]
    # Data is already log-transformed: compute difference
    FC <- sweep(X, 2, as.numeric(X_control), FUN = "-")
    as.data.frame(FC) %>%
      rownames_to_column("sample") %>%
      filter(sample != control) %>%
      pivot_longer(cols = !sample, names_to = "Kinase_Name", values_to = "LogFC")
  }
}

binarize_matrix <- function(mat, quantile_threshold = 0.5) {
  # Validate input
  if (!is.matrix(mat)) {
    stop("Input must be a matrix")
  }
  if (quantile_threshold < 0 || quantile_threshold > 1) {
    stop("Quantile threshold must be between 0 and 1")
  }
  
  # Calculate the threshold value
  threshold_value <- quantile(mat, probs = quantile_threshold, na.rm = TRUE)
  
  # Binarize the matrix
  binarized_mat <- matrix(as.integer(mat >= threshold_value), 
                          nrow = nrow(mat), 
                          ncol = ncol(mat))
  
  # Preserve row and column names if they exist
  rownames(binarized_mat) <- rownames(mat)
  colnames(binarized_mat) <- colnames(mat)
  
  return(binarized_mat)
}


compare_uka_cs_comp <- function(X, uka_comp, control = "Pooled REF", fscore_cutoff = 2, 
                                del_cell = c('A3KAW', "DB", "HT", "MC116")){
  
  uka_cs <- X_to_uka(X, control = control) %>% filter(!sample %in% del_cell)
  print(paste0("UKA - Confidence score - n kins: ", nrow(uka_cs)))
  uka_comp_to_control <- uka_comp %>%
    clean_uka_to_kinograte(control = control, fscore_cutoff = fscore_cutoff, del_cell = del_cell)
  print(paste0("UKA - legacy - n kins: ", nrow(uka_comp_to_control)))
  check <- uka_cs %>% full_join(uka_comp_to_control, 
                                by = c("sample" = "cell_line", "Kinase_Name" = "uniprotname"), 
                                suffix = c(".cs", ".comp"))  %>% distinct()
  check_filt <- check %>% filter(!is.na(LogFC.cs)) %>% filter(!is.na(LogFC.comp) ) %>%
    filter(abs(LogFC.comp) > 0.1 & abs(LogFC.cs) > 0.1)
  print(paste0("n common kins in confidence-score based UKA and legacy: ", nrow(check_filt)))
  browser()
  plot(check_filt$LogFC.cs, check_filt$LogFC.comp, xlab = "UKA - Conf score - FC", ylab  = "UKA - legacy - FC")
}

evaluate_X_vs_uka <- function(X, uka_df, control = "Pooled REF", del_cell = c('A3KAW', "DB", "HT", "MC116")) {
  # Returns a numeric similarity score between X and uka_df (e.g., mean Pearson correlation across cell lines)
  # X: [sample x kinase] matrix
  # uka_df: data.frame with columns cell_line, uniprotname, LogFC
  # Returns: mean Pearson correlation of LogFC for overlapping kinases per cell line

  # Convert X to long format (LogFC vs control)
  X_long <- X_to_uka(X, control = control) %>%
    filter(!sample %in% del_cell) %>%
    rename(cell_line = sample, uniprotname = Kinase_Name, LogFC_X = LogFC)

  # Prepare uka_df
  uka_long <- uka_df %>%
    clean_uka_to_kinograte(control = control, del_cell = del_cell) %>%
    select(cell_line, uniprotname, LogFC) %>%
    rename(LogFC_uka = LogFC)

  # Merge on cell_line and uniprotname
  merged <- inner_join(X_long, uka_long, by = c("cell_line", "uniprotname"))

  # For each cell_line, compute Pearson correlation of LogFC_X and LogFC_uka
  corrs <- merged %>%
    group_by(cell_line) %>%
    summarise(cor = ifelse(length(unique(LogFC_X)) > 1 & length(unique(LogFC_uka)) > 1,
                           cor(LogFC_X, LogFC_uka, use = "pairwise.complete.obs"), NA_real_)) %>%
    pull(cor)

  # Return mean correlation (excluding NA)
  mean_corr <- mean(corrs, na.rm = TRUE)
  return(mean_corr)
}

# Helper: generate scores_df from a parameter vector
scores_df_from_vector <- function(param_vec, template_scores_df) {
  # param_vec: numeric vector of cscore values (same order as template_scores_df)
  # template_scores_df: data.frame with columns Database, Kinase_Rank
  stopifnot(length(param_vec) == nrow(template_scores_df))
  out <- template_scores_df
  out$cscore <- param_vec
  return(out)
}

# Wrapper: given a parameter vector, compute X and evaluate similarity to uka_02
objective_X_vs_uka <- function(param_vec, template_scores_df, min_score, min_score_p, db, kinaseFamily, H, uka_02, control = "Pooled REF") {
  scores_df <- scores_df_from_vector(param_vec, template_scores_df)
  K <- make_K(scores_df = scores_df, min_score = min_score, min_score_p = min_score_p, db = db, kinaseFamily = kinaseFamily)
  X <- make_X(S = H, K = K)
  score <- evaluate_X_vs_uka(X, uka_02, control = control)
  print(score)
  # For maximization, return score; for minimization, return -score
  return(score)
}

