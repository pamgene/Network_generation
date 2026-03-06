library(kinograte)
library(visNetwork)
library(igraph)
library(httr)
library(ggplot2)
library(PCSF)
source("R/kinograte_PG_core.R")
library(data.table)

augment_by_threshold_steps <- function(pcsf_graph, full_graph, cost_threshold = 0.1, steps = 1) {
  # Get the core nodes from the PCSF result
  core_nodes <- V(pcsf_graph)$name
  print(paste0("number of input nodes: ", length(unique(core_nodes))))
  
  # Start with core nodes
  current_nodes <- core_nodes0
  
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






compute_network_metrics_kinase <- function(nodes, edges, missing_nodes, 
                               res.path, condition, perc_cutoff){
  
  g <- graph_from_data_frame(edges, directed = F, vertices = nodes)
  
  make_shortestpath_all <- function(g, missing_nodes, relative_to = "n_nodes", 
                                    inverse_score = FALSE, add_missing = T, a = 1){
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

    #  inverse shortest path so that higher score means shorter
    if (inverse_score) {
      all_shortestpaths <- 1/(a + all_shortestpaths)
    }
    
    return(median(all_shortestpaths, na.rm = T) / n_nodes)
  }
  
  rel_med_path = make_shortestpath_all(g = g, missing_nodes = missing_nodes)
  rel_med_path_inv = make_shortestpath_all(g = g, missing_nodes = missing_nodes, inverse_score = TRUE)
  
  list(
    rel_med_path = rel_med_path,
    rel_med_path_inv = rel_med_path_inv,
    density = igraph::edge_density(g, loops = FALSE),
    clustering = mean(igraph::transitivity(g, type = "local"))
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



make_network_and_stats <- function(uka, art_nodes = NULL, art_lfc = NULL, perc_cutoff, fscore_cutoff, res.path, condition = NULL, write = F, 
ppi_network,  b, wp_ontology_names = NULL, highlight_degree = 5) {
  # Makes a network then computes network statistics.
  # if write = T, enrich network with pathways and plot the graph with the pathways.
  
  # 0. Artificial nodes
  if (!is.null(art_nodes)){
    art_df <- data.frame(name = art_nodes, prize = 1, type = "Artificial", LogFC = art_lfc)
    uka <- bind_rows(uka, art_df)
  }
  
  # 1. Make network
  kinograte_res <- kinograte_pg(df = uka, ppi_network = ppi_network, 
                                perc_cutoff = perc_cutoff,
                                res.path = res.path, condition = condition,
                                cluster = T, b = b, write = write)
  
  
  g <- kinograte_res$network
  metrics <- compute_network_metrics_kinase(kinograte_res$nodes, kinograte_res$edges, kinograte_res$missing_nodes, 
                                            res.path = res.path, condition = condition, perc_cutoff = perc_cutoff)
  
  
  if (write == T){
    
    # save node_cluster data
    
    nodes_clusters <- kinograte_res$wc_df %>% 
      group_by(cluster) %>%
      summarize(nodes = paste0(id, collapse = "; "))
    write_csv(nodes_clusters, paste0(res.path, "/nodes_clusters_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv"))
    
    # 4. Enrich network with pathways, then network with pathways
    enrich_file <- paste0(res.path, "/enrich_results_per_pathway_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv")
    if (!file.exists(enrich_file)) {
      pws <- do_network_enrichment(kinograte_res$network, pval = 0.05,
                                   perc_cutoff = perc_cutoff, folder = res.path, condition = condition,
                                   database = c("Reactome_Pathways_2024"), wp_ontology_names = wp_ontology_names)
      # dbs: WikiPathways_2024_Human, KEGG_2021_Human, Reactome_Pathways_2024
    } else {
      pws <- readr::read_csv(enrich_file, show_col_types = FALSE)
    }
    prot2pw <-  pws %>%
      dplyr::select(Group_Pw, Genes) %>%
      separate_longer_delim(Genes, delim = ";") %>%
      distinct() %>%
      dplyr::rename("Protein" = "Genes") %>%
      group_by(Protein) %>%
      summarize(pathway = paste0(Group_Pw, collapse = ","))
    
    write_csv(prot2pw, paste0(res.path, "/nodes_", condition, "_p", sub(".*0\\.", "", perc_cutoff), "_with_pathways.csv"))
    
    
    nodes2pw <- kinograte_res$nodes  %>%
      left_join(prot2pw, by = "Protein")
    
    # Check if we have valid data before creating visualization
    if (nrow(nodes2pw) == 0 || is.null(kinograte_res$edges) || nrow(kinograte_res$edges) == 0) {
      warning(paste("No nodes or edges available for pathway visualization of", condition, "- skipping"))
      mynet <- NULL
    } else {
      mynet <- visualize_network_pg(nodes = nodes2pw,
                                    edges = kinograte_res$edges,
                                    maintitle = kinograte_res$maintitle,
                                    cluster_df = kinograte_res$wc_df, options_by = 'pathway',
                                    highlight_degree = highlight_degree)
    }
    mynet
    
    # Check if network object is valid before saving
    if (is.null(mynet) || !inherits(mynet, "visNetwork")) {
      warning(paste("Invalid network object for", condition, "_pws - skipping HTML save"))
    } else {
      # Ensure output directory exists and use normalized path
      output_dir <- normalizePath(res.path, mustWork = FALSE)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      condition <- gsub(" ", "_", condition)
      html_file <- file.path(output_dir, paste0(condition, "_p", perc_cutoff, ".html"))
      
      
      # Use htmlwidgets::saveWidget instead of visSave for more robust saving
      tryCatch({
        htmlwidgets::saveWidget(mynet, 
                                file = html_file,
                                selfcontained = TRUE,
                                background = "white")
        cat("Successfully saved:", html_file, "\n")
      }, error = function(e) {
        warning(paste("Failed to save HTML for", condition, "_pws:", e$message))
        # Try alternative method if htmlwidgets fails
        tryCatch({
          visSave(mynet, html_file, selfcontained = F, background = "white")
          cat("Saved with visSave (not self-contained):", html_file, "\n")
        }, error = function(e2) {
          warning(paste("Both saveWidget and visSave failed for", condition, "_pws:", e2$message))
        })
      })
    }
  }
  return(metrics)
  
}




