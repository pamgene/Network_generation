library(kinograte)
library(visNetwork)
library(igraph)
library(httr)
library(PCSF)


extract_data_percentile <- function(data_list, res.path, perc_cutoff, spec_cutoff, rank_lowest_highest = FALSE){
  # data list includes a named list of dataframes, 
  # where names should be from the following: Kinase, Pep, MS, RNA
  # dataframe columns: uniprotname, LogFC + column for kinases: fscore 
  
  data_top_list <- list()
  
  for (i in 1:length(data_list)){
    data_name <- names(data_list[i])
    #ranking <- ranking_vect[i]
    # 1. rank data according to value column
    
    # Sensitivity data should not be ranked based on absolute LFC. 
    if(data_name != "Sensitivity"){
      data_rank <- percentile_rank(data_list[[i]], symbol = uniprotname, metric = LogFC)
    } else {
      data_rank <- percentile_rank_pg_noabs(data_list[[i]], uniprotname, LogFC, rank_lowest_highest = rank_lowest_highest)
    }
    
    # extract top hits
    data_top <- top_hits_pg(data_rank, perc_cutoff = perc_cutoff, spec_cutoff = spec_cutoff,
                            omic_type = data_name)
    data_top_list[[data_name]] <- data_top
  }
  
  combined_df <- bind_rows(data_top_list, .id = "type") %>% 
    dplyr::group_by(name) %>%
    summarize(
      type = paste0(type, collapse = "-"),
      prize = mean(prize),
      # direction = sign(LogFC),
      LogFC_all = paste0(round(LogFC, 2), collapse = ", "),
      LogFC = mean(LogFC)
    ) %>%
    distinct()
  
  # combined_df$direction[grep(",", combined_df$LogFC)] <- 0
  sumdf <- combined_df %>% group_by(type) %>% summarize(n_nodes = n())
  print(sumdf) 
  
  return(combined_df)
}



top_hits_pg <- function(df, perc_cutoff, omic_type, spec_cutoff){
  if (!omic_type == "Kinase"){
    df %>% dplyr::filter(score >= perc_cutoff) %>% dplyr::ungroup() %>% 
      dplyr::mutate(rank = 1:base::nrow(.), prize = score, 
                    type = omic_type) %>% 
      dplyr::select(name, prize, type, LogFC)
  } else {
    df %>% dplyr::filter(fscore >= spec_cutoff) %>% dplyr::ungroup() %>% 
      dplyr::mutate(rank = 1:base::nrow(.), prize = score, 
                    type = omic_type) %>% 
      dplyr::select(name, prize, type, LogFC)
  }
  
}

PCSF_rand_pg <- function(edges_df, terminals, n = 8, r = 0.1, 
                         w = 2, b = 1, mu = 0.0005, dummies = NULL) {
  # set.seed(12345)

  if (missing(edges_df))
    stop("Need to specify edge dataframe with columns: head, tail, cost")
  if (!all(c("head", "tail", "cost") %in% colnames(edges_df)))
    stop("edges_df must have columns: head, tail, cost")
  if (missing(terminals))
    stop("  Need to provide terminal nodes as a named numeric vector, 
    where node names must be same as in the interaction network.")
  if (is.null(names(terminals)))
    stop("  The terminal nodes must be provided as a named numeric vector, 
    where node names must be same as in the interaction network.")
  
  # Gather terminal information
  terminal_names <- names(terminals)
  terminal_values <- as.numeric(terminals)
  
  # === REPLICATE construct_interactome EXACTLY ===
  node_names <- unique(c(as.character(edges_df[, 1]), as.character(edges_df[, 2])))
  # Build temporary igraph to get simplified structure
  # This replicates: graph.data.frame() + simplify()
  temp_graph <- igraph::graph.data.frame(
    edges_df[, 1:2], 
    vertices = node_names, 
    directed = FALSE
  )
  igraph::E(temp_graph)$weight <- as.numeric(edges_df[, 3])
  temp_graph <- igraph::simplify(temp_graph)
  
  # calculate degrees from simplified graph
  node_degrees <- igraph::degree(temp_graph)
  # Ensure order matches node_names
  node_degrees <- node_degrees[node_names]
  # Get simplified edges for the algorithm
  edges_simplified <- igraph::ends(temp_graph, es = igraph::E(temp_graph))
  edge_weights_simplified <- igraph::E(temp_graph)$weight
  # To test equivalence:
  # write_csv(as.data.frame(edge_weights_simplified), "results/pcsf_rand_equivalent/pcsf_rand_pg_edgew_simplified.csv")
  
  # Create node prize vector
  node_prz <- vector(mode = "numeric", length = length(node_names))
  
  # Match terminals to nodes
  index <- match(terminal_names, node_names)
  percent <- signif((length(index) - sum(is.na(index))) / length(index) * 100, 4)
  
  if (percent < 5)
    stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
  
  cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
  
  terminal_names <- terminal_names[!is.na(index)]
  terminal_values <- terminal_values[!is.na(index)]
  index <- index[!is.na(index)]
  node_prz[index] <- terminal_values
  
  # Set dummies
  if (missing(dummies) || is.null(dummies) || is.na(dummies))
    dummies <- terminal_names
  
  cat("  Solving the PCSF by adding random noise to the edge costs...\n")
  
  # Calculate hub penalization
  hub_penalization <- -mu * node_degrees
  
  # Update node prizes
  node_prizes <- b * node_prz
  index <- which(node_prizes == 0)
  node_prizes[index] <- hub_penalization[index]
  
  # === PREPARE EDGES - use simplified edges ===
  from <- c(rep("DUMMY", length(dummies)), edges_simplified[, 1])
  to <- c(dummies, edges_simplified[, 2])
  base_costs <- edge_weights_simplified
  
  all_nodes <- NULL


  for (i in 1:n) {
    # set.seed(123)
    random_noise <- stats::runif(length(base_costs), 0, r)
    cat("First 3 random values:", head(random_noise, 3), "\n")
    
    # write_csv(as.data.frame(random_noise), "results/pcsf_rand_equivalent/pcsf_rand_pg_random_vals.csv")
    cost <- c(
      rep(w, length(dummies)), 
      base_costs + (base_costs*random_noise)
    )
    # write_csv(as.data.frame(cost), "results/pcsf_rand_equivalent/pcsf_rand_pg_cost.csv")
    
    
    # set.seed(123)
    output <- call_sr(from, to, cost, node_names, node_prizes)
    
    edge <- data.frame(
      as.character(output[[1]]),
      as.character(output[[2]]),
      stringsAsFactors = FALSE
    )
    colnames(edge) <- c("source", "target")
    
    edge <- edge[which(edge[, 1] != "DUMMY"), ]
    edge <- edge[which(edge[, 2] != "DUMMY"), ]
    
    assign(paste0("graph_", i), edge)
    all_nodes <- c(all_nodes, edge$source, edge$target)
  }
  
  # === CALCULATE STATISTICS ===
  node_frequency <- table(all_nodes)
  node_names_out <- names(node_frequency)
  node_prizes_out <- as.numeric(node_frequency)
  
  # === BUILD ADJACENCY MATRIX ===
  adj_matrix <- matrix(0, length(node_names_out), length(node_names_out))
  colnames(adj_matrix) <- node_names_out
  rownames(adj_matrix) <- node_names_out
  
  for (i in 1:n) {
    graph <- get(paste0("graph_", i))
    edges <- graph
    
    x <- match(edges[, 1], node_names_out)
    y <- match(edges[, 2], node_names_out)
    
    if (length(x) > 0 & length(y) > 0) {
      for (j in 1:length(x)) {
        if (x[j] >= y[j]) {
          k <- x[j]
          l <- y[j]
        } else {
          k <- y[j]
          l <- x[j]
        }
        adj_matrix[k, l] <- adj_matrix[k, l] + 1
      }
    }
  }
  
  # === CHECK AND RETURN ===
  # browser()
  # write_csv(as.data.frame(adj_matrix) %>% rownames_to_column("from"), "results/pcsf_rand_equivalent/adj_mat_new.csv")
  if (sum(adj_matrix) != 0) {
    
    # Convert to edge list
    edge_indices <- which(adj_matrix > 0, arr.ind = TRUE)
    result_edges <- data.frame(
      from = node_names_out[edge_indices[, 1]],
      to = node_names_out[edge_indices[, 2]],
      weight = adj_matrix[edge_indices],
      stringsAsFactors = FALSE
    )
    
    # Create nodes dataframe
    index <- match(node_names_out, node_names_out)
    result_nodes <- data.frame(
      name = node_names_out,
      prize = node_prizes_out[index],
      type = ifelse(node_names_out %in% terminal_names, "Terminal", "Steiner"),
      stringsAsFactors = FALSE
    )
    result <- list(
      edges = result_edges,
      nodes = result_nodes,
      n_runs = n
    )
    
    class(result) <- c("PCSF_df", "list")
    return(result)
    
  } else {
    stop("  Subnetwork can not be identified for a given parameter set.
    Provide a compatible b or mu value with your terminal prize list...\n\n")
  }
}

call_sr <- function(from, to, cost, node_names, node_prizes) {
  # write_csv(as.data.frame(from), "results/pcsf_rand_equivalent/call_sr_from_new.csv")
  # write_csv(as.data.frame(to), "results/pcsf_rand_equivalent/call_sr_to_new.csv")
  # write_csv(as.data.frame(node_names), "results/pcsf_rand_equivalent/call_sr_node_names_new.csv")
  # write_csv(as.data.frame(node_prizes), "results/pcsf_rand_equivalent/call_sr_node_prizes_new.csv")
  .Call('_PCSF_call_sr', PACKAGE = 'PCSF', from, to, cost, node_names, node_prizes)
}

kinograte_pg_pcsf <- function(df, ppi_network, maintitle, n = 8, w = 10, r = 0.1, b = 1.5, 
                         mu = 0.005, cluster = TRUE, seed = NULL, res.path, perc_cutoff, 
                         condition, write) {
  
  print('Using fast kinograte...')
  maintitle <- paste0(condition, " - Network with percentile rank cutoff = ", 
                      perc_cutoff, ", Number of nodes = ", nrow(df))
  
  terms <- df$prize
  names(terms) <- df$name
  
  if (!is.null(seed)) {
    set.seed(12345)
  }
  
  # Ensure ppi_network is a dataframe
  ppi_network <- base::as.data.frame(ppi_network)

  # === CALL THE NEW DATAFRAME-BASED PCSF ===
  
  # set.seed(123)
  subnet <- tryCatch({
    PCSF_rand_pg(
      edges_df = ppi_network,  # Assuming columns: head, tail, cost
      terminals = terms,
      n = n,
      w = w,
      r = r,
      b = b,
      mu = mu
    )
  }, error = function(e) {
    message("PCSF_rand_df error: ", e$message)
    return(NULL)
  })
  
  # write_csv(subnet$edges, "results/pcsf_rand_equivalent/pcsf_rand_pg_output.csv")
  
  if (is.null(subnet)) {
    return(NULL)
  }
  
  # === EXTRACT EDGES AND NODES FROM DATAFRAME OUTPUT ===
  edges <- subnet$edges  # Already a dataframe with: from, to, weight
  
  # Get unique nodes from edges
  nodes <- dplyr::tibble(
    Protein = c(edges$from %>% base::unique(), edges$to %>% base::unique())
  ) %>% 
    dplyr::distinct()
  
  # Join with original df to get prizes and types
  nodes <- dplyr::left_join(nodes, df, by = c("Protein" = "name")) %>% 
    dplyr::mutate(
      type = ifelse(is.na(type), "Hidden", type), 
      prize = ifelse(is.na(prize), 2, prize)
    )
  
  # === CALCULATE DEGREES ===
  if (!'degree' %in% colnames(nodes)) {
    # Calculate degree from edge list
    degree_from <- table(edges$from)
    degree_to <- table(edges$to)
    all_degrees <- c(degree_from, degree_to)
    degree_sum <- tapply(all_degrees, names(all_degrees), sum)
    
    node_degrees <- dplyr::tibble(
      Protein = names(degree_sum),
      degree = as.numeric(degree_sum)
    )
    
    nodes <- dplyr::left_join(nodes, node_degrees, by = "Protein")
  }
  
  # === CHECK FOR MISSING NODES ===
  subnet_nodes <- unique(c(edges$from, edges$to))
  missing_nodes <- setdiff(names(terms), subnet_nodes)
  
  if (length(missing_nodes) > 0 && write) {
    missing_df <- df %>% dplyr::filter(name %in% missing_nodes)
    readr::write_csv(
      missing_df, 
      paste0(res.path, "/missing_nodes_", condition, "_p", 
             sub(".*0\\.", "", perc_cutoff), ".csv")
    )
  } else {
    missing_df <- NULL
  }
  
  # === PREPARE OUTPUT NODES ===
  kinograte_nodes <- nodes %>% 
    as.data.frame() %>% 
    dplyr::select(-any_of(c('title')))
  
  nodes_to_write <- kinograte_nodes %>% select(-prize)
  
  if (write) {
    readr::write_csv(
      nodes_to_write, 
      paste0(res.path, "/nodes_", condition, "_p", 
             sub(".*0\\.", "", perc_cutoff), ".csv")
    )
  }
  
  # === ADD TITLE COLUMN ===
  if (!"pathway" %in% colnames(nodes)) {
    nodes$title <- if ("LogFC_all" %in% colnames(nodes)) {
      paste0("LFC = ", nodes$LogFC_all)
    } else {
      paste0("LFC = ", round(nodes$LogFC, 2))
    }
  } else {
    nodes$title <- nodes$pathway
  }
  
  # === CLUSTERING (if requested) ===
  if (cluster) {
    # Build igraph object from edge list for clustering
    subnet_igraph <- igraph::graph_from_data_frame(
      d = edges[, c("from", "to", "weight")],
      directed = FALSE
    )
    
    wc <- suppressWarnings(
      igraph::edge.betweenness.community(
        subnet_igraph, 
        weights = igraph::E(subnet_igraph)$weight, 
        directed = FALSE, 
        bridges = TRUE
      )
    )
    
    wc_df <- dplyr::tibble(
      id = wc$names, 
      cluster = wc$membership
    )
    
    return(list(
      network = subnet_igraph,  # Return igraph for compatibility
      nodes = nodes, 
      edges = edges, 
      missing = missing_df,
      wc_df = wc_df, 
      maintitle = maintitle
    ))
  } else {
    # Build igraph object for compatibility even without clustering
    subnet_igraph <- igraph::graph_from_data_frame(
      d = edges[, c("from", "to", "weight")],
      directed = FALSE
    )
    
    return(list(
      network = subnet_igraph, 
      nodes = nodes, 
      edges = edges, 
      missing_nodes = missing_df, 
      maintitle = maintitle
    ))
  }
}

kinograte_pg <- function (df, ppi_network, maintitle, n = 8, w = 10, r = 0.1, b = 1.5, 
                          mu = 0.005, cluster = TRUE, seed = NULL, res.path, perc_cutoff, condition, write) {
  
  maintitle = paste0(condition, " - Network with percentile rank cutoff = ", perc_cutoff, ", Number of nodes = ", nrow(df))
  
  terms <- df$prize
  names(terms) <- df$name
  ppi_net <- base::as.data.frame(ppi_network)
  ppi <- PCSF::construct_interactome(ppi_net)
  if (!is.null(seed)) {
    base::set.seed(12345)
  }
  # browser()
  # set.seed(123)
  subnet <- tryCatch({
    PCSF_rand(ppi, terms, n = n, w = w, r = r, b = b, mu = mu)
  }, error = function(e) {
    message("PCSF_rand error: ", e$message)
    return(NULL)
  })
  # equivalence_check
  # subnet_to_save <- as_data_frame(subnet)
  # write_csv(subnet_to_save, "results/pcsf_rand_equivalent/pcsf_rand_output.csv")
  # 
  if (is.null(subnet) || igraph::vcount(subnet) == 0 || igraph::ecount(subnet) == 0) {
    # If error, skip and return NULL so the calling for loop can continue
    message("PCSF_rand returned empty network.")
    return(NULL)
  }
  edges <- igraph::as_data_frame(subnet, what = c("edges", 
                                                  "vertices", "both"))

  nodes <- dplyr::tibble(Protein = c(edges$from %>% base::unique(), 
                                     edges$to %>% base::unique())) %>% dplyr::distinct()

  nodes <- dplyr::left_join(nodes, df, by = c(Protein = "name")) %>% 
    dplyr::mutate(type = ifelse(is.na(type), "Hidden", type), 
                  prize = ifelse(is.na(prize), 2, prize))

  if(!'degree' %in% colnames(nodes)){
    # this is the case when kinograte_pg is run "natively", not on selected nodes
    node_degrees <- dplyr::tibble(
      Protein = igraph::V(subnet)$name,
      degree = igraph::degree(subnet)
    )
    # Join degrees into the nodes tibble
    nodes <- dplyr::left_join(nodes, node_degrees, by = "Protein")
  } 

  # Print nodes in terms but not in subnet
  subnet_nodes <- igraph::V(subnet)$name
  missing_nodes <- setdiff(names(terms), subnet_nodes)
  if (length(missing_nodes) > 0 && write) {
    missing_df <- df %>% filter(name %in% missing_nodes)
    write_csv(missing_df, paste0(res.path, "/missing_nodes_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv"))
  } else {
    missing_df = NULL
  }

  kinograte_nodes <- nodes %>% as.data.frame() %>% 
    dplyr::select(-any_of(c('title')))
  
  nodes_to_write <- kinograte_nodes %>% select(-prize)
  
  if (write) {
    readr::write_csv(
      nodes_to_write, 
      paste0(res.path, "/nodes_", condition, "_p", 
             sub(".*0\\.", "", perc_cutoff), ".csv")
    )
  }
  

  if (!"pathway" %in% colnames(nodes)){
    # Use LogFC_all if present, otherwise fallback to LogFC
    nodes$title <- if ("LogFC_all" %in% colnames(nodes)) {
      paste0("LFC = ", nodes$LogFC_all)
    } else {
      paste0("LFC = ", round(nodes$LogFC, 2))
    }
  } else {
    nodes$title <- nodes$pathway
  }

  if (cluster) {
    wc <- suppressWarnings(igraph::edge.betweenness.community(subnet, 
                                                              weights = igraph::E(subnet)$weight, directed = FALSE, 
                                                              bridges = TRUE))
    wc_df <- dplyr::tibble(id = wc$names, cluster = wc$membership)
    
    return(list(network = subnet, nodes = nodes, edges = edges, missing = missing_df,
                wc_df = wc_df, maintitle = maintitle))
  }
  else {
    return(list(network = subnet, nodes = nodes, edges = edges, missing_nodes = missing_df, maintitle = maintitle))
  }
}


PCSF_rand <-
  function(ppi, terminals, n = 8, r = 0.1, w = 2, b = 1, mu = 0.0005,dummies){
    # set.seed(12345)
    # Checking function arguments
    if (missing(ppi))
      stop("Need to specify an interaction network \"ppi\".")
    if (class(ppi) != "igraph")
      stop("The interaction network \"ppi\" must be an igraph object.")
    if (missing(terminals))
      stop("  Need to provide terminal nodes as a named numeric vector, 
    where node names must be same as in the interaction network.")
    if(is.null(names(terminals)))
      stop("  The terminal nodes must be provided as a named numeric vector, 
    where node names must be same as in the interaction network.")
    
    
    # Gather the terminal genes to be analyzed, and their scores
    terminal_names = names(terminals)
    terminal_values = as.numeric(terminals)
    
    # Incorporate the node prizes
    node_names = V(ppi)$name
    node_prz = vector(mode = "numeric", length = length(node_names))
    index = match(terminal_names, node_names)
    percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
    if (percent < 5)
      stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
    cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
    terminal_names = terminal_names[!is.na(index)]
    terminal_values = terminal_values[!is.na(index)]
    index = index[!is.na(index)]
    node_prz[index] =  terminal_values
    
    ## Prepare input file for MST-PCSF implementation in C++
    if(missing(dummies)||is.null(dummies)||is.na(dummies))
      dummies = terminal_names #re-assign this to allow for input
    
    cat("  Solving the PCSF by adding random noise to the edge costs...\n")
    
    # Calculate the hub penalization scores
    node_degrees = igraph::degree(ppi)
    hub_penalization = - mu*node_degrees
    
    
    # Update the node prizes
    node_prizes = b*node_prz
    index = which(node_prizes==0)
    node_prizes[index] = hub_penalization[index]
    
    # Construct the list of edges 
    edges = ends(ppi,es = E(ppi))
    # test equivalence of pcsf_rand_pg with this
    # write_csv(as.data.frame(E(ppi)$weight), "results/pcsf_rand_equivalent/pcsf_rand_edgew.csv")
    from = c(rep("DUMMY", length(dummies)), edges[,1])
    to = c(dummies, edges[,2])
    

    all_nodes=NULL
    
    # Run the MST-PCSF algorithm for n times with random noise added to edge costs at each time
    for(i in 1:n){

      # Randomize the edge costs
      # set.seed(123)
      random_vals = stats::runif(length(E(ppi)), 0, r)
      # write_csv(as.data.frame(random_vals), "results/pcsf_rand_equivalent/pcsf_rand_random_vals.csv")
      cat("First 3 random values:", head(random_vals, 3), "\n")
      
      cost = c(rep(w, length(dummies)), E(ppi)$weight + E(ppi)$weight*random_vals)
      # write_csv(as.data.frame(cost), "results/pcsf_rand_equivalent/pcsf_rand_cost.csv")
      # Feed in the input files into MSt-PCSF algorithm
      
      # set.seed(123)
      output = call_sr(from,to,cost,node_names,node_prizes)
      
      # Construct an igraph object for current parameter set, and save it
      edge = data.frame(as.character(output[[1]]),as.character(output[[2]]))
      colnames(edge)= c("source", "target")
      edge = edge[which(edge[,1]!="DUMMY"), ]
      edge = edge[which(edge[,2]!="DUMMY"), ]
      
      graph = graph.data.frame(edge,directed=F)
      assign(paste0("graph_",i), get("graph"))
      all_nodes = c(all_nodes, V(graph)$name)
    }
    
    # Calculate graph statistics
    node_frequency = table(all_nodes)
    print(head(node_frequency, 10))
    
    node_names = names(node_frequency)
    node_prizes =  as.numeric(node_frequency)
    # Combine the graphs in order to get unionized graph
    adj_matrix = matrix(0, length(node_names), length(node_names))
    colnames(adj_matrix) = node_names
    rownames(adj_matrix) = node_names
    for(i in 1:n){
      assign("graph", get(paste0("graph_",i)))
      edges = ends(graph,es = E(graph))
      x = match(edges[,1],node_names)
      y = match(edges[,2],node_names)
      if(length(x)>0 & length(y)>0){
        for( j in 1:length(x)){
          if(x[j]>=y[j]){
            k = x[j]
            l = y[j]
          }else{
            k = y[j]
            l = x[j]
          }
          adj_matrix[k,l] = adj_matrix[k,l]+1
        }
      }
    }
    # write_csv(as.data.frame(adj_matrix) %>% rownames_to_column("from"), "results/pcsf_rand_equivalent/adj_mat_orig.csv")
    # Check the size of output subnetwork and print a warning if it is 0
    if(sum(adj_matrix) != 0){
      
      # Construct the igraph object from the union graph

      subnet = graph_from_adjacency_matrix(adj_matrix, weighted=TRUE, mode="undirected")
      index = match(V(subnet)$name, node_names)
      V(subnet)$prize = node_prizes[index]
      
      #   # Associate the type of nodes to shape
      #   V(subnet)$shape = "triangle"
      #   index = match(terminal_names, V(subnet)$name)
      #   index = index[!is.na(index)]
      #   V(subnet)$shape[index] = "circle"
      
      # Associate the type of nodes
      V(subnet)$type = "Steiner"
      index = match(terminal_names, V(subnet)$name)
      index = index[!is.na(index)]
      V(subnet)$type[index] = "Terminal"
      
      
      class(subnet) <- c("PCSF", "igraph")
      
      return (subnet)
      
    } else {
      
      stop("  Subnetwork can not be identified for a given parameter set.
    Provide a compatible b or mu value with your terminal prize list...\n\n")
    }
    
    
  }



percentile_rank_pg_noabs <- function(df, symbol, metric, rank_lowest_highest = FALSE) {
  # rank_lowest_highest: If T, the lowest values get the highest rank (i.e., most sensitive/lowest LogFC is ranked highest).
  df %>% 
    dplyr::filter(!is.na({{ metric }})) %>% 
    dplyr::rename(name = {{ symbol }}) %>% 
    dplyr::distinct(name, .keep_all = T) -> df
  if(rank_lowest_highest) {
    df <- dplyr::mutate(df, score = dplyr::percent_rank(dplyr::desc({{ metric }})))
  } else {
    df <- dplyr::mutate(df, score = dplyr::percent_rank({{ metric }}))
  }
  
  dplyr::arrange(df, desc(score)) 
}




