library(kinograte)
library(visNetwork)
library(igraph)
library(httr)
library(ggplot2)
library(PCSF)
source("R/kinograte_PG_core.R")

augment_by_threshold_steps <- function(pcsf_graph, full_graph, cost_threshold = 0.1, steps = 1) {
  # Get the core nodes from the PCSF result
  core_nodes <- V(pcsf_graph)$name
  print(paste0("number of input nodes: ", length(unique(core_nodes))))
  
  # Start with core nodes
  current_nodes <- core_nodes
  
  # Iteratively expand by threshold-filtered edges for 'steps' iterations
  for (step in 1:steps) {
    print(paste0("Step ", step, " of ", steps))
    
    # Find all edges connected to current nodes with cost <= threshold
    current_vids <- which(V(full_graph)$name %in% current_nodes)
    
    # Get all edges incident to current nodes
    incident_edges <- incident_edges(full_graph, current_vids, mode = "all")
    incident_edges <- unique(unlist(incident_edges))
    
    # Filter edges by cost threshold
    strong_edges <- incident_edges[E(full_graph)[incident_edges]$cost <= cost_threshold]
    
    # Get all vertices incident to these strong edges
    strong_verts <- ends(full_graph, strong_edges, names = TRUE)
    strong_verts <- unique(as.vector(strong_verts))
    
    # Update current_nodes for next iteration
    current_nodes <- unique(c(current_nodes, strong_verts))
    
    print(paste0("After step ", step, ", nodes count: ", length(current_nodes)))
  }
  
  # Induce the final subgraph
  augmented_graph <- induced_subgraph(full_graph, which(V(full_graph)$name %in% current_nodes))
  print(paste0("Number of nodes in final augmented graph: ", length(unique(V(augmented_graph)$name))))
  
  return(augmented_graph)
}

filter_for_core_network <- function(g) {
  # excludes kinases that are not on a path between two sensitivity or sensitivity-kinase nodes
  # Get node types (assuming stored as vertex attribute 'type')
  node_types <- V(g)$type
  
  
  # Get indices of different node types
  k_nodes <- which(node_types == "Kinase")
  t_nodes <- which(node_types == "Sensitivity") 
  d_nodes <- which(node_types == "Sensitivity-Kinase")
  td_nodes <- c(t_nodes, d_nodes)
  h_nodes <- which(node_types == "Hidden")
  
  # Find K nodes that are on paths between T/D nodes
  k_nodes_to_keep <- c()
  h_nodes_to_keep <- c()
  
  # Check all pairs of T/D nodes
  for(i in 1:(length(td_nodes)-1)) {
    for(j in (i+1):length(td_nodes)) {
      node1 <- td_nodes[i]
      node2 <- td_nodes[j]
      
      # Find all shortest paths between this pair
      paths <- all_shortest_paths(g, from = node1, to = node2, mode = "all")
      
      # Extract K nodes from these paths
      for(path in paths$res) {
        path_nodes <- as.numeric(path)
        k_in_path <- intersect(path_nodes, k_nodes)
        h_in_path <- intersect(path_nodes, h_nodes)
        k_nodes_to_keep <- c(k_nodes_to_keep, k_in_path)
        h_nodes_to_keep <- c(h_nodes_to_keep, h_in_path)
      }
    }
  }
  
  # Remove duplicates
  k_nodes_to_keep <- unique(k_nodes_to_keep)
  h_nodes_to_keep <- unique(h_nodes_to_keep)
  
  # Find K and H nodes to remove
  k_nodes_to_remove <- setdiff(k_nodes, k_nodes_to_keep)
  h_nodes_to_remove <- setdiff(h_nodes, h_nodes_to_keep)
  nodes_to_remove <- c(k_nodes_to_remove, h_nodes_to_remove)
  
  # Remove these nodes from the graph
  if(length(nodes_to_remove) > 0) {
    g_filtered <- delete_vertices(g, nodes_to_remove)
  } else {
    g_filtered <- g
  }
  
  return(list(
    graph = g_filtered,
    removed_kinase_nodes = V(g)$name[k_nodes_to_remove],
    removed_hidden_nodes = V(g)$name[h_nodes_to_remove],
    kept_kinase_nodes = V(g)$name[k_nodes_to_keep],
    kept_hidden_nodes = V(g)$name[h_nodes_to_keep]
  ))
}



# compute path statistics between kinases and targets
compute_path_stats_kin_tar <- function(nodes, edges, missing_nodes, a = 1, relative_to, res.path, condition, perc_cutoff) {
  # Check if edges has at least two columns
  if (is.null(edges) || !is.data.frame(edges) || ncol(edges) < 2) {
    message("make_network_stats: edges data frame does not have at least two columns, skipping.")
    return(list(
      rel_med_path_KT = NA,
      rel_med_path_KT_core = NA,
      rel_med_path_all_nodes = NA,
      rel_med_path_all_nodes_core = NA,
      rel_med_path_KT_vs_all_nodes = NA,
      rel_med_path_KT_vs_all_nodes_core = NA,
      rel_med_path_KT_core_inv = NA,
      rel_med_path_all_core_inv = NA,
      rel_med_path_KT_vs_all_nodes_core_inv = NA
    ))
  }
  g <- graph_from_data_frame(edges, directed = F, vertices = nodes)
  g_core <- filter_for_core_network(g)$graph

  save_network_plot <- function(g, maintitle, res.path, condition, perc_cutoff){
    nodes <- igraph::as_data_frame(g, what = "vertices") %>% dplyr::rename("Protein" = 'name')
    # Add degree column to nodes
    nodes$degree <- igraph::degree(g)
    edges <- igraph::as_data_frame(g, what = "edges")
    mynet <- visualize_network_pg(nodes = nodes, edges = edges,
                               maintitle = maintitle, cluster_df = NULL,
                               options_by = 'pathway', highlight_degree = 5)
    mynet
    # Ensure directory exists before saving
    if (!dir.exists(res.path)) {
      dir.create(res.path, recursive = TRUE, showWarnings = FALSE)
    }

    # Save with error handling to prevent HTML widget issues
    html_file <- paste0(res.path, "/", maintitle, "_p", perc_cutoff, ".html")
    tryCatch({
      visSave(mynet, html_file, selfcontained = TRUE, background = "white")
    }, error = function(e) {
      warning(paste("Failed to save HTML file:", html_file, "Error:", e$message))
      # Try saving without selfcontained as fallback
      tryCatch({
        visSave(mynet, html_file, selfcontained = FALSE, background = "white")
      }, error = function(e2) {
        warning(paste("Fallback save also failed:", e2$message))
      })
    })
  }

  # save_network_plot(g, maintitle = paste0("Full network - ", condition), res.path, condition, perc_cutoff)
  # save_network_plot(g_core, maintitle = paste0(condition, "_core network"), res.path, condition, perc_cutoff)


  # --- Metric 1: median path between kinase and target nodes ---
  make_shortestpath_kin_sens <- function(g, missing_nodes, relative_to, add_missing = T, inverse_score = FALSE, a = 1){
    # relative to = "n_nodes" (this includes missing_nodes nodes) or "n_sens_all" (sensitivity and double nodes)
    nodes_kin <- V(g)[V(g)$type == "Kinase"]
    nodes_sens <- V(g)[V(g)$type == "Sensitivity"]
    nodes_double <- V(g)[V(g)$type == "Sensitivity-Kinase"]
    nodes_sd <- V(g)[V(g)$type != "Kinase"]
    
    # Compute shortest paths from kinases to sensitivity and double nodes
    shortestpaths_K_TD <- distances(g, v = nodes_kin, to = nodes_sd, mode = "all") %>% 
      as.data.frame() %>% rownames_to_column("from") %>% 
      pivot_longer(cols = !"from", names_to = "to", values_to = "d")
    # count double nodes as kinase
    if (length(nodes_double) > 0){
      shortestpaths_T_D <- distances(g, v = nodes_sens, to = nodes_double, mode = "all") %>% 
        as.data.frame() %>% rownames_to_column("from") %>% 
        pivot_longer(cols = !"from", names_to = "to", values_to = "d")
      shortestpaths_kin_tar <- bind_rows(shortestpaths_K_TD, shortestpaths_T_D) %>%
        distinct() %>% pull(d) 
    } else {
      shortestpaths_kin_tar <- shortestpaths_K_TD %>% pull(d)
    }   
    # Add overlap nodes as shortest path = 0
    if (length(nodes_double) > 0){
      shortestpaths_kin_tar <- c(shortestpaths_kin_tar,  rep(0, length(nodes_double))) 
    }
    # Add missing nodes as maximum path length
    if (!is.null(missing_nodes) && nrow(missing_nodes) > 0 && add_missing) {
      n_nodes <- vcount(g) + nrow(missing_nodes)
      max_path_kin_tar <- max(shortestpaths_kin_tar, na.rm = TRUE)
      shortestpaths_kin_tar <- c(shortestpaths_kin_tar, rep(max_path_kin_tar, nrow(missing_nodes)))
    } else {
      n_nodes <- vcount(g)
    }
    #  inverse shortest path so that higher score means shorter
    if (inverse_score) {
      shortestpaths_kin_tar <- 1/(a + shortestpaths_kin_tar)
    }

    if (relative_to == "n_sens_all"){
      relative_to <- length(nodes_sd)
    } else if (relative_to == "n_nodes"){
      relative_to <- n_nodes
    } else {
      stop("parameter defining what the network score is relative to is missing")
    }
    
    return(median(shortestpaths_kin_tar, na.rm = T) / relative_to)
  }

  # --- Metric 2: median path between all nodes ---
  make_shortestpath_all <- function(g, missing_nodes, relative_to = "n_nodes", inverse_score = FALSE, add_missing = T, a = 1){
    nodes_sd <- V(g)[V(g)$type != "Kinase"]
    all_shortestpaths <- as.vector(distances(g, mode = "all"))
    all_shortestpaths[is.infinite(all_shortestpaths)] <- NA
    # Add missing_nodes nodes as maximum path length
    if (!is.null(missing_nodes) && nrow(missing_nodes) > 0) {
      n_nodes <- vcount(g) + nrow(missing_nodes)
      max_path_all <- max(all_shortestpaths, na.rm = TRUE)
      all_shortestpaths <- c(all_shortestpaths, rep(max_path_all, nrow(missing_nodes)))
    } else {
      n_nodes <- vcount(g)
    }
    if (relative_to == "n_sens_all"){
      relative_to <- length(nodes_sd)
    } else if (relative_to == "n_nodes"){
      relative_to <- n_nodes
    } else {
      stop("parameter defining what the network score is relative to is missing")
    }
    #  inverse shortest path so that higher score means shorter
    if (inverse_score) {
      all_shortestpaths <- 1/(a + all_shortestpaths)
    }

    return(median(all_shortestpaths, na.rm = T) / relative_to)
  }
  rel_med_path_KT = make_shortestpath_kin_sens(g = g, missing_nodes = missing_nodes, relative_to = relative_to)
  rel_med_path_KT_core = make_shortestpath_kin_sens(g = g_core, missing_nodes = missing_nodes, relative_to = relative_to, add_missing = F)
  rel_med_path_all_nodes = make_shortestpath_all(g = g, missing_nodes = missing_nodes, relative_to = relative_to)
  rel_med_path_all_nodes_core = make_shortestpath_all(g = g_core, missing_nodes = missing_nodes, relative_to = relative_to, add_missing = F)
  rel_med_path_KT_core_inv = make_shortestpath_kin_sens(g = g, missing_nodes = missing_nodes, relative_to = relative_to, add_missing = F, inverse_score = TRUE)
  rel_med_path_all_core_inv = make_shortestpath_all(g = g, missing_nodes = missing_nodes, relative_to = relative_to, add_missing = F, inverse_score = TRUE)

  # Compute relative median paths
  list(
    rel_med_path_KT = rel_med_path_KT,
    rel_med_path_KT_core = rel_med_path_KT_core,
    rel_med_path_all_nodes = round(rel_med_path_all_nodes, 3), # rounding is here not in the make_shortestpath_all bc the value is used for division
    rel_med_path_all_nodes_core = round(rel_med_path_all_nodes_core, 3),
    rel_med_path_KT_vs_all_nodes = round(rel_med_path_KT / rel_med_path_all_nodes, 3),
    rel_med_path_KT_vs_all_nodes_core = round(rel_med_path_KT_core / rel_med_path_all_nodes_core, 3),
    rel_med_path_KT_core_inv = rel_med_path_KT_core_inv,
    rel_med_path_all_core_inv = rel_med_path_all_core_inv,
    rel_med_path_KT_vs_all_nodes_core_inv = round(rel_med_path_KT_core_inv / rel_med_path_all_core_inv, 5)
  )
}

# Helper: average local clustering for Sensitivity and Sensitivity-Kinase nodes
compute_avg_sens_clustering <- function(g, nodes_df = NULL) {
  # nodes_df: data.frame with columns 'Protein' and 'type'
  # If nodes_df is not provided, try to extract from g (will fail if not present)
  if (is.null(nodes_df)) {
    stop("nodes_df must be provided for correct node type mapping.")
  }
  # Get node names for Sensitivity and Sensitivity-Kinase
  sens_names <- nodes_df %>%
    filter(type %in% c("Sensitivity", "Sensitivity-Kinase")) %>%
    pull(Protein)
  # Only keep those present in the graph
  sens_names <- sens_names[sens_names %in% V(g)$name]
  if (length(sens_names) == 0) return(NA)
  sens_nodes <- V(g)[name %in% sens_names]
  sens_clustering <- igraph::transitivity(g, type = "local", vids = sens_nodes)
  sens_clustering <- sens_clustering[!is.nan(sens_clustering)]
  if (length(sens_clustering) > 0) {
    mean(sens_clustering)
  } else {
    NA
  }
}

compute_network_metrics <- function(g, nodes, edges, missing_nodes, relative_to, res.path = NULL, condition = NULL, perc_cutoff = NULL) {
  path_stats <- compute_path_stats_kin_tar(nodes = nodes, edges = edges, missing_nodes = missing_nodes, relative_to = relative_to, 
                                          res.path = res.path, condition = condition, perc_cutoff = perc_cutoff)
  avg_sens_clustering <- compute_avg_sens_clustering(g, nodes_df = nodes)
  c(
    path_stats,
    avg_sens_clustering = avg_sens_clustering,
    density = igraph::edge_density(g, loops = FALSE),
    modularity = tryCatch({
      comm <- igraph::cluster_edge_betweenness(g)
      igraph::modularity(comm)
    }, error = function(e) NA),
    assortativity = tryCatch(igraph::assortativity_degree(g, directed = FALSE), error = function(e) NA)
  )
}

compute_network_metrics_kinase <- function(nodes, edges, missing_nodes, relative_to,
                               res.path, condition, perc_cutoff){
  
  g <- graph_from_data_frame(edges, directed = F, vertices = nodes)
  
  make_shortestpath_all <- function(g, missing_nodes, relative_to = "n_nodes", 
                                    inverse_score = FALSE, add_missing = T, a = 1){
    nodes_sd <- V(g)[V(g)$type != "Kinase"]
    all_shortestpaths <- as.vector(distances(g, mode = "all"))
    all_shortestpaths[is.infinite(all_shortestpaths)] <- NA
    # Add missing_nodes nodes as maximum path length
    if (!is.null(missing_nodes) && nrow(missing_nodes) > 0) {
      n_nodes <- vcount(g) + nrow(missing_nodes)
      max_path_all <- max(all_shortestpaths, na.rm = TRUE)
      all_shortestpaths <- c(all_shortestpaths, rep(max_path_all, nrow(missing_nodes)))
    } else {
      n_nodes <- vcount(g)
    }
    if (relative_to == "n_sens_all"){
      relative_to <- length(nodes_sd)
    } else if (relative_to == "n_nodes"){
      relative_to <- n_nodes
    } else {
      stop("parameter defining what the network score is relative to is missing")
    }
    #  inverse shortest path so that higher score means shorter
    if (inverse_score) {
      all_shortestpaths <- 1/(a + all_shortestpaths)
    }
    
    return(median(all_shortestpaths, na.rm = T) / relative_to)
  }
  
  rel_med_path_all_nodes = make_shortestpath_all(g = g, missing_nodes = missing_nodes, relative_to = relative_to)
  
  list(
    rel_med_path_all_nodes = rel_med_path_all_nodes,
    rel_med_path_all_nodes_core = NA,
    rel_med_path_all_core_inv = NA,
    rel_med_path_KT = NA,
    rel_med_path_KT_core = NA,
    rel_med_path_KT_vs_all_nodes = NA,
    rel_med_path_KT_vs_all_nodes_core = NA,
    rel_med_path_KT_core_inv = NA,
    rel_med_path_KT_vs_all_nodes_core_inv = NA,
    avg_sens_clustering = NA,
    density = NA,
    modularity = NA,
    assortativity = NA
  )
  
}

# Helper: plot histogram of random metrics with observed value
plot_network_metric_hist <- function(metric_name, random_vals, obs_val, res.path, cell, perc_cutoff) {
  dir.create(file.path(res.path, "Network_metrics_plots"), showWarnings = FALSE, recursive = TRUE)
  df <- data.frame(value = random_vals)
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "#BDBDBD", color = "black") +
    geom_vline(xintercept = obs_val, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste0(metric_name, " (", cell, ", perc=", perc_cutoff, ")"),
         x = metric_name, y = "Frequency") +
    theme_bw()
  fname <- file.path(res.path, "Network_metrics_plots",
                     paste0("hist_", metric_name, "_", cell, "_p", perc_cutoff, ".png"))
  ggsave(fname, p, width = 6, height = 4)
}

make_ukadb_net_pws <- function(ukadb, ppi_network, slk, res.path, b = 2, highlight_degree = 10 ){
  # slk: signalink protein to pathway
  if (!dir.exists(res.path)) {
    dir.create(res.path, recursive = TRUE, showWarnings = FALSE)
  }
  
  ukadb_to_kinograte <- data.frame(name = unique(ukadb$Kinase_Name),
                                   prize = 1,
                                   LogFC = 1,
                                   type = "Kinase")
  
  kinograte_res <- kinograte_pg_pcsf(df = ukadb_to_kinograte, ppi_network = ppi_network, 
                                     perc_cutoff = 0,
                                     res.path = res.path, condition = "UKAdb",
                                     cluster = T, b = b, write = T)
  nodes2pw <- kinograte_res$nodes %>%
    left_join(slk, by = c("Protein" = "name")) %>% 
    dplyr::rename(pathway = pathways) %>%
    mutate(pathway = gsub("\\|", ",", pathway),
           pathway = gsub("\\/", "-", pathway)) %>%
    # group_by(Protein, type) %>%
    # summarize(pathway = paste0(pathway, collapse = ","))  %>% # this shows that hidden nodes are not from wnt mostly :(
    distinct(Protein, .keep_all = T) %>% mutate(pathway = "WNT")
  
  mynet <- visualize_network_pg(nodes = nodes2pw,
                                edges = kinograte_res$edges,
                                maintitle = kinograte_res$maintitle,
                                cluster_df = kinograte_res$wc_df, options_by = 'pathway',
                                highlight_degree = highlight_degree)
  mynet
  
  # Save with error handling to prevent HTML widget issues
  html_file <- paste0(res.path, "/WNT_ukadb_kinases_network_SignaLink_pws.html")
  visSave(mynet, html_file, selfcontained = TRUE, background = "white")
  
  
}

make_network_and_stats <- function(uka, sens, perc_cutoff, spec_cutoff, res.path, condition = NULL, 
                                   write = F, ppi_network, relative_to, b,
                                   wp_ontology_names = NULL, highlight_degree = 5) {
  # Makes a network from kinase and sensitivity data.
  # Then computes network statistics.
  # if write = T, enrich network with pathways and plot the graph with the pathways.
  # 1. combine uka and sens for the network 
  combined_df <- bind_rows(list(Sensitivity = sens, Kinase = uka), .id = "type") %>% 
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
  # 2. Make network
  kinograte_res <- kinograte_pg_pcsf(df = combined_df, ppi_network = ppi_network, 
                                perc_cutoff = perc_cutoff,
                                res.path = res.path, condition = condition,
                                cluster = T, b = b, write = write)
  
  # 3. Compute network stats 
  if (!is.null(kinograte_res$network)) {
    g <- kinograte_res$network
    metrics <- compute_network_metrics(g, kinograte_res$nodes, kinograte_res$edges, kinograte_res$missing_nodes, relative_to = relative_to,
                                       res.path = res.path, condition = condition, perc_cutoff = perc_cutoff)
    
    
    # 4. Enrich network with pathways, then network with pathways
    if (write == T){
      enrich_file <- paste0(res.path, "/pathways_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv")
      if (!file.exists(enrich_file)) {
        pws <- do_network_enrichment(kinograte_res$network, pval = 0.05,
                            perc_cutoff = perc_cutoff, folder = res.path, condition = condition,
                            database = c("WikiPathways_2024_Human"), wp_ontology_names = wp_ontology_names)
      } else {
        pws <- readr::read_csv(enrich_file, show_col_types = FALSE)
      }
      
      prot2pw <-  pws %>%
        dplyr::select(Pathway, Genes) %>%
        separate_longer_delim(Genes, delim = ";") %>%
        distinct() %>%
        dplyr::rename("Protein" = "Genes") %>%
        group_by(Protein) %>%
        summarize(pathway = paste0(Pathway, collapse = ","))
      
      write_csv(prot2pw, paste0(res.path, "/nodes_", condition, "_p", sub(".*0\\.", "", perc_cutoff), "_with_pathways.csv"))
      

      nodes2pw <- kinograte_res$nodes  %>%
        left_join(prot2pw, by = "Protein")

      mynet <- visualize_network_pg(nodes = nodes2pw,
                               edges = kinograte_res$edges,
                               maintitle = kinograte_res$maintitle,
                               cluster_df = kinograte_res$wc_df, options_by = 'pathway',
                               highlight_degree = highlight_degree)
      mynet

      # Save with error handling to prevent HTML widget issues
      html_file <- paste0(res.path, "/", condition, "_network_p", perc_cutoff, "_pws.html")
      tryCatch({
        visSave(mynet, html_file, selfcontained = TRUE, background = "white")
      }, error = function(e) {
        warning(paste("Failed to save HTML file:", html_file, "Error:", e$message))
        # Try saving without selfcontained as fallback
        tryCatch({
          visSave(mynet, html_file, selfcontained = FALSE, background = "white")
        }, error = function(e2) {
          warning(paste("Fallback save also failed:", e2$message))
        })
      })
    }
    return(metrics)
  }
  return(NULL)
}




# compute_med_path_kin_tar <- function(nodes, edges, perc_cutoff, res.path, condition = NULL, write = FALSE){
#   # Compute Sensitivity node interactors ratio and network stats, return as data.frame
#   # mean_n_interactor_ratio <- get_sensitivity_interactors(
#   #   kinograte_nodes = nodes, 
#   #   ppi_network = ppi_networkv12, 
#   #   res.path = res.path, 
#   #   perc_cutoff = perc_cutoff,
#   #   condition = condition
#   # )
#   # Nodes that were input to the network but were dropped are also needed

#   path_stats <- compute_path_stats_kin_tar(nodes, edges, missing_nodes, perc_cutoff, res.path, condition)
#   rel_med_path_kinase_target <- path_stats$rel_med_path_kinase_target
#   rel_med_path_all_nodes <- path_stats$rel_med_path_all_nodes
#   rel_med_path_kinase_target_vs_all_nodes <- path_stats$rel_med_path_kinase_target_vs_all_nodes
#   n_nodes <- path_stats$n_nodes

#   stats <- data.frame(
#     condition = condition,
#     Percentile.cutoff = perc_cutoff,
#     n_nodes = n_nodes,
#     n_kin = nrow(nodes %>% filter(type == "Kinase")),
#     n_target = nrow(nodes %>% filter(type == "Sensitivity")),
#     n_double = nrow(nodes %>% filter(type == "Sensitivity-Kinase")),
#     Metric = c("Relative med path kinase-target", "Relative med path all nodes", "Kinase-target / all nodes path ratio"), 
#     value = c(rel_med_path_kinase_target, rel_med_path_all_nodes, rel_med_path_kinase_target_vs_all_nodes)
#   )
#   if (write){
#     write_csv(stats, paste0(res.path, "/stats_", condition, ".csv")) 
#   }
#   # Metric = c("Relative med path kinase-target", "Relative med path all nodes",
#   #            "n overlaps / n targets", "mean of n interactors / refnet interactors"), 
#   # value = c(rel_med_path_kinase_target, rel_med_path_all_nodes, 
#   #           rel_overlap_kin_sens, mean_n_interactor_ratio)
#   return(list(
#     rel_med_path_kinase_target = rel_med_path_kinase_target,
#     rel_med_path_all_nodes = rel_med_path_all_nodes,
#     rel_med_path_kinase_target_vs_all_nodes = rel_med_path_kinase_target_vs_all_nodes
#   ))
# }







get_sensitivity_interactors <- function(kinograte_nodes, ppi_network, res.path, perc_cutoff, condition) {
  # Filter for Sensitivity or Sensitivity-Kinase nodes
  sens_nodes <- kinograte_nodes %>%
    dplyr::filter(type %in% c("Sensitivity", "Sensitivity-Kinase")) %>%
    dplyr::pull(Protein)
  
  if (length(sens_nodes) == 0) {
    message("No Sensitivity or Sensitivity-Kinase nodes found.")
    return(invisible(NULL))
  }
  
  # Interactors in kinograte result network
  # Build edge list from kinograte_nodes (use all nodes in kinograte_nodes)
  # But the kinograte result edges are not available here, so use degree column
  kinograte_degrees <- kinograte_nodes %>%
    dplyr::filter(Protein %in% sens_nodes) %>%
    dplyr::select(Protein, n_interactors_net = degree)
  
  # Interactors in reference network
  # ppi_network is assumed to have columns head, tail
  ppi_edges <- ppi_network
  if (!all(c("head", "tail") %in% colnames(ppi_edges))) {
    stop("ppi_network must have columns 'head' and 'tail'")
  }
  # For each sens node, count unique interactors in ref net
  n_interactors_refnet <- sapply(sens_nodes, function(prot) {
    sum(ppi_edges$head == prot | ppi_edges$tail == prot)
  })
  refnet_df <- data.frame(Protein = sens_nodes, n_interactors_refnet = n_interactors_refnet)
  
  # Merge with kinograte degrees
  result <- dplyr::left_join(refnet_df, kinograte_degrees, by = "Protein")
  result$n_interactors_net[is.na(result$n_interactors_net)] <- 0
  result$ratio <- round((result$n_interactors_net / result$n_interactors_refnet), 5)
  result <- distinct(result)
  
  result_med <- round(mean(result$ratio), 3)
  result_med_row <- data.frame(Protein = "Mean", n_interactors_refnet = NA, 
                              n_interactors_net = NA, ratio = result_med)
  
  result <- rbind(result, result_med_row)
  
  # Save
  write_csv(
    result,
    paste0(res.path, "/stats_interactors_of_sens_nodes_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv")
  )
  
  return(result_med)
}




