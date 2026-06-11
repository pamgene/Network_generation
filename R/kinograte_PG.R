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
                               res.path, condition, spec_cutoff){
  
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
plot_network_metric_hist <- function(metric_name, random_vals, obs_val, res.path, cell, spec_cutoff) {
  dir.create(file.path(res.path, "Network_metrics_plots"), showWarnings = FALSE, recursive = TRUE)
  df <- data.frame(value = random_vals)
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "#BDBDBD", color = "black") +
    geom_vline(xintercept = obs_val, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = paste0(metric_name, " (", cell, ", spec=", spec_cutoff, ")"),
         x = metric_name, y = "Frequency") +
    theme_bw()
  fname <- file.path(res.path, "Network_metrics_plots",
                     paste0("hist_", metric_name, "_", cell, "_spec", spec_cutoff, ".png"))
  ggsave(fname, p, width = 6, height = 4)
}


make_network_and_stats <- function(uka, sens, art_nodes = NULL, art_lfc = NULL, spec_cutoff, res.path, condition = NULL, 
                                   write = F, ppi_network, relative_to, b,
                                   highlight_degree = 5) {
  # Makes a network from kinase and sensitivity data.
  # Then computes network statistics.
  # if write = T, enrich network with pathways and plot the graph with the pathways.
  
  # 0. Artificial nodes
  if (!is.null(art_nodes)){
    art_df <- data.frame(name = art_nodes, prize = 1, type = "Artificial", LogFC = art_lfc)
    uka <- bind_rows(uka, art_df)
  }
  # 1. combine uka and sens for the network.
  # LogFC_all is used in network visualization as a label that pops up when hovering mouse. LogFC is used to color the node.
  combined_df <- bind_rows(list(Sensitivity = sens, Kinase = uka), .id = "type") %>% 
    dplyr::group_by(name) %>%
    summarize(
      type = paste0(type, collapse = "-"),
      prize = mean(prize),
      LogFC_all = paste0(round(LogFC, 2), collapse = ", ")
    ) %>%
    distinct() %>%
    left_join(uka[,c('name', 'LogFC')], by = "name") %>%
    left_join(sens[,c('name', 'LogFC')] %>% dplyr::rename(LogFC_s = LogFC), by = "name") %>%
    mutate(LogFC = round(coalesce(LogFC, LogFC_s), 2)) %>%
    select(-LogFC_s)
  

  # combined_df$direction[grep(",", combined_df$LogFC)] <- 0
  sumdf <- combined_df %>% group_by(type) %>% summarize(n_nodes = n())
  print(sumdf) 
  # 2. Make network
  kinograte_res <- kinograte_pg_pcsf(df = combined_df, ppi_network = ppi_network, 
                                     spec_cutoff = spec_cutoff,
                                     res.path = res.path, condition = condition,
                                     cluster = T, b = b, write = write)
  
  # 3. Compute network stats 
  if (!is.null(kinograte_res$network)) {
    # g <- kinograte_res$network
    # metrics <- compute_network_metrics(g, kinograte_res$nodes, kinograte_res$edges, kinograte_res$missing_nodes, relative_to = relative_to,
    #                                    res.path = res.path, condition = condition, spec_cutoff = spec_cutoff)
    
    
    # 4. Enrich network with pathways
    if (write == T){
      
      nodes_clusters <- kinograte_res$wc_df %>% 
        group_by(cluster) %>%
        summarize(nodes = paste0(id, collapse = "; "))
      write_csv(nodes_clusters, paste0(res.path, "/nodes_clusters_", condition, "_spec", spec_cutoff, ".csv"))
      
      enrich_file <- paste0(res.path, "/pathways_", condition, "_spec", spec_cutoff, ".csv")
      if (!file.exists(enrich_file)) {
        pws <- do_network_enrichment(kinograte_res$network, pval = 0.05,
                                     spec_cutoff = spec_cutoff, folder = res.path, condition = condition,
                                     database = c("Reactome_Pathways_2024"), min_n_hits = min_n_hits,
                                     nodes = kinograte_res$nodes)
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
      
      write_csv(prot2pw, paste0(res.path, "/nodes_", condition, "_spec", spec_cutoff, "_with_pathways.csv"))
      
      
      nodes2pw <- kinograte_res$nodes  %>%
        left_join(prot2pw, by = "Protein")

      # # this is to select proteins for AACR abstract, WSUDCL2 figure
      # filtered_proteins <- c("PIK3CA", "MTOR", "AKT1", "RPS6KB1", "EIF2AK3", "EIF2AK1", "IFIT1B",
      #                        "CSNK2A1", "PRKACA", "CSNK1E", "PLK1", "TTK",
      #                        "FBXO43", "CAMK2A", "CAMK2KD", "PNCK", "ATF1", "RPS6KA5")
      # # this is to select proteins for AACR abstract, NUDUL figure
      # filtered_proteins <- c("CHEK1", "ATR", "ATM", "MAPK8", "MAPK9", "MAPK13", "MAPK3",
      #                        "CDK2", "CDK4", "CDK6", "CDK3", "CDK1", "AURKA", "PLK3", "CDK7",
      #                        "ORC2", "MAK", "MASTL", "VRK1")
      # nodes2pw <- nodes2pw %>% filter(Protein %in% filtered_proteins)
      
      
      mynet <- visualize_network_pg(nodes = nodes2pw,
                                    edges = kinograte_res$edges,
                                    maintitle = kinograte_res$maintitle,
                                    cluster_df = kinograte_res$wc_df, options_by = 'pathway',
                                    highlight_degree = 10)
      mynet
      
      # Save with error handling to prevent HTML widget issues
      html_file <- paste0(res.path, "/", condition, "_network_spec", spec_cutoff, "_pws.html")
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
  }
}

make_network_and_stats_kinase <- function(uka, art_nodes = NULL, art_lfc = NULL, spec_cutoff,
                                          res.path, condition = NULL, write = F,
ppi_network,  b, highlight_degree = 5, min_n_hits = 2) {
  # Makes a network then computes network statistics.
  # if write = T, enrich network with pathways and plot the graph with the pathways.
  
  # 0. Artificial nodes
  if (!is.null(art_nodes)){
    art_df <- data.frame(name = art_nodes, prize = 1, type = "Artificial", LogFC = art_lfc)
    uka <- bind_rows(uka, art_df)
  }
  
  # 1. Make network, or load from disk if already done
  nodes_file <- paste0(res.path, "/nodes_", condition, "_spec", spec_cutoff, ".csv")

  edges_file <- paste0(res.path, "/edges_", condition, "_spec", spec_cutoff, ".csv")
  wc_file    <- paste0(res.path, "/wc_",    condition, "_spec", spec_cutoff, ".csv")

  if (file.exists(nodes_file) && file.exists(edges_file)) {
    cat("\nNetwork already exists for ", condition, " spec ", spec_cutoff, ", loading from disk.\n")
    nodes <- readr::read_csv(nodes_file, show_col_types = FALSE)
    edges <- readr::read_csv(edges_file, show_col_types = FALSE)
    wc_df <- if (file.exists(wc_file)) readr::read_csv(wc_file, show_col_types = FALSE) else NULL
    missing_df <- {
      mf <- paste0(res.path, "/missing_nodes_", condition, "_spec", spec_cutoff, ".csv")
      if (file.exists(mf)) readr::read_csv(mf, show_col_types = FALSE) else NULL
    }
    maintitle <- paste0(condition, " - Network with specificity cutoff = ", spec_cutoff,
                        ", Number of nodes = ", nrow(nodes))
    subnet <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = nodes$Protein))
    kinograte_res <- list(network = subnet, nodes = nodes, edges = edges,
                          missing = missing_df, wc_df = wc_df, maintitle = maintitle)
  } else {
    kinograte_res <- tryCatch(
      {
        kinograte_pg(df = uka, ppi_network = ppi_network,
                     spec_cutoff = spec_cutoff,
                     res.path = res.path, condition = condition,
                     cluster = T, b = b, write = write)
      }, error = function(e) {
        cat("\nCan't generate network for ", condition, " with spec cutoff ", spec_cutoff, ": ", e$message, "\n")
        return(NULL)
      }
    )
    if (is.null(kinograte_res)) return(NULL)
  }

  g <- kinograte_res$network
  metrics <- compute_network_metrics_kinase(kinograte_res$nodes, kinograte_res$edges, kinograte_res$missing_nodes,
                                            res.path = res.path, condition = condition, spec_cutoff = spec_cutoff)
  
  
  if (write == T){
    
    # save node_cluster data
    
    if (!is.null(kinograte_res$wc_df)) {
      nodes_clusters <- kinograte_res$wc_df %>%
        group_by(cluster) %>%
        summarize(nodes = paste0(id, collapse = "; "))
      write_csv(nodes_clusters, paste0(res.path, "/nodes_clusters_", condition, "_spec", spec_cutoff, ".csv"))
    }
    
    # 4. Enrich network with pathways, then network with pathways
    enrich_file <- paste0(res.path, "/pathways_", condition, "_spec", spec_cutoff, ".csv")
    if (!file.exists(enrich_file)) {
      pws <- tryCatch(
        {
          do_network_enrichment(kinograte_res$network, pval = 0.05,
                                spec_cutoff = spec_cutoff, folder = res.path, condition = condition,
                                database = c("Reactome_Pathways_2024"), min_n_hits = min_n_hits,
                                nodes = kinograte_res$nodes)
        }, error = function(e) {
            # Handle the error
            cat("\nCan't generate pathway enrichment for condition: ", condition, " and spec cutoff: ", spec_cutoff, "\n")
        }
        
      )
      # dbs: WikiPathways_2024_Human, KEGG_2021_Human, Reactome_Pathways_2024
    } else {
      pws <- readr::read_csv(enrich_file, show_col_types = FALSE)
    }
    
    if (!is.null(pws)){
      prot2pw <-  pws %>%
        dplyr::select(Group_Pw, Genes) %>%
        separate_longer_delim(Genes, delim = ";") %>%
        distinct() %>%
        dplyr::rename("Protein" = "Genes") %>%
        group_by(Protein) %>%
        summarize(pathway = paste0(Group_Pw, collapse = ","))
      
      write_csv(prot2pw, paste0(res.path, "/nodes_", condition, "_spec", spec_cutoff, "_with_pathways.csv"))
      
      nodes2pw <- kinograte_res$nodes  %>%
        left_join(prot2pw, by = "Protein")
    }
    
    # Skip visualization if HTML already exists
    output_dir <- normalizePath(res.path, mustWork = FALSE)
    condition_clean <- gsub(" ", "_", condition)
    html_file <- file.path(output_dir, paste0(condition_clean, "_spec", spec_cutoff, ".html"))

    if (!file.exists(html_file)) {
      if (exists("nodes2pw") && nrow(nodes2pw) > 0){
        # all data available
        mynet <- visualize_network_pg(nodes = nodes2pw,
                                      edges = kinograte_res$edges,
                                      maintitle = kinograte_res$maintitle,
                                      cluster_df = kinograte_res$wc_df, options_by = 'pathway',
                                      highlight_degree = highlight_degree)
      } else {
        # network available, pathways not
        mynet <- visualize_network_pg(nodes = kinograte_res$nodes,
                                      edges = kinograte_res$edges,
                                      maintitle = kinograte_res$maintitle,
                                      cluster_df = kinograte_res$wc_df, options_by = 'group',
                                      highlight_degree = highlight_degree)
      }

      # Check if network object is valid before saving
      if (is.null(mynet) || !inherits(mynet, "visNetwork")) {
        warning(paste("Invalid network object for", condition, "_pws - skipping HTML save"))
      } else {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
      
      
        # Use htmlwidgets::saveWidget instead of visSave for more robust saving
        tryCatch({
          htmlwidgets::saveWidget(mynet,
                                  file = html_file,
                                  selfcontained = TRUE,
                                  background = "white")
          cat("\nSuccessfully saved:", html_file, "\n")
        }, error = function(e) {
          warning(paste("Failed to save HTML for", condition, "_pws:", e$message))
          tryCatch({
            visSave(mynet, html_file, selfcontained = F, background = "white")
            cat("Saved with visSave (not self-contained):", html_file, "\n")
          }, error = function(e2) {
            warning(paste("Both saveWidget and visSave failed for", condition, "_pws:", e2$message))
          })
        })
      } # end valid mynet
    } # end !file.exists(html_file)
  } # end write == T
  return(metrics)
  
}




