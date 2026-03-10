# Network enrichment and visualization

library(kinograte)
library(visNetwork)
library(igraph)
library(httr)
library(ggplot2)
library(rWikiPathways)
library(httr)
library(jsonlite)
library(ggraph)


visualize_network_pg <- function (
    nodes, edges, maintitle,
    cluster_df = NULL,
    layout = "layout_with_fr",
    seed = 123,
    options_by = "group",
    highlight_nodes = NULL,
    highlight_degree = 5
) {
  
  edges_vis <- edges
  # Set edge colors to black
  edges_vis$color <- "#000000"
  
  # Prepare nodes
  nodes_vis <- nodes %>%
    dplyr::rename(id = Protein, group = type, value = prize) %>%
    dplyr::distinct() %>%
    dplyr::mutate(label = id)
  # Add cluster info
  if (!is.null(cluster_df)) {
    nodes_vis <- dplyr::left_join(nodes_vis, cluster_df, by = "id")
  }
  # If options_by is pathway and pathway column exists, create pathway_multi for multi-select
  if (!is.null(options_by) && options_by == "pathway" && "pathway" %in% colnames(nodes_vis)) {
    nodes_vis <- nodes_vis %>%
      dplyr::mutate(
        pathway_multi = ifelse(
          is.na(pathway) | pathway == "",
          NA_character_,
          gsub("\\s*,\\s*", ";", gsub('["()]', '', trimws(as.character(pathway))))
        )
      )
  }
  
  
  
  # Node shape based on type
  type_shapes <- c(
    "Kinase" = "dot",
    "Artificial" = "star",
    "Peptide" = "ellipse",
    "Peptide-Kinase" = "box",
    "Sensitivity" = "square",
    "Sensitivity-Kinase" = 'hexagon',
    # "RNA" = "square",
    # "RNA-Peptide" = 'hexagon',
    "RNA-Kinase"= "diamond",
    'RNA-Peptide-Kinase' = 'star',
    "Protein" = "triangle",
    "Protein-Kinase" = "diamond",
    "Protein-Peptide" = "hexagon",
    "Protein-Peptide-Kinase" = "star",
    "Hidden" = "text"
  )
  
  nodes_vis$shape <- type_shapes[as.character(nodes_vis$group)]
  nodes_vis$shape[is.na(nodes_vis$shape)] <- "ellipse"
  
  # Adaptive diverging color scale: blue -> gray -> red
  logfc_vals <- nodes_vis$LogFC %>% as.numeric()
  logfc_vals[is.infinite(logfc_vals)] <- NA
  
  p25_neg <- stats::quantile(logfc_vals[logfc_vals < 0], probs = 0.25, na.rm = TRUE)
  p75_pos <- stats::quantile(logfc_vals[logfc_vals > 0], probs = 0.75, na.rm = TRUE)
  
  limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
  color_limits <- c(-limit_val, limit_val)
  
  diverge_pal <- colorRampPalette(c("#0000CC", "#BDBDBD", "#C62828"))  # blue -> gray -> red
  col_vals <- diverge_pal(201)  # 201 values, midpoint at 101
  
  # Clamp and scale LogFC
  logfc_clamped <- pmax(pmin(logfc_vals, color_limits[2]), color_limits[1])
  scaled_fc <- scales::rescale(logfc_clamped, to = c(-100, 100), from = color_limits)
  scaled_fc[is.na(scaled_fc)] <- 0
  color_index <- round(scaled_fc) + 101
  nodes_vis$color.background <- col_vals[color_index]
  
  # Highlight high-degree nodes
  if (!is.null(highlight_degree)) {
    nodes_vis <- nodes_vis %>%
      dplyr::mutate(
        color.border = ifelse(degree >= highlight_degree, "#4C9900", "#454545"),
        borderWidth = ifelse(degree >= highlight_degree, 2.5, 1)
      )
  } else {
    nodes_vis$color.border <- "#454545"
    nodes_vis$borderWidth <- 1
  }
  
  # Build network
  net <- visNetwork::visNetwork(nodes = distinct(nodes_vis, id, .keep_all = TRUE), edges = edges_vis, main = maintitle) %>%
    visNetwork::visPhysics() %>%
    visNetwork::visNodes(font = list(size = 10)) %>%
    visNetwork::visEdges(smooth = FALSE) %>%
    visNetwork::visIgraphLayout(layout = layout, randomSeed = seed) %>%
    visNetwork::visInteraction(multiselect = TRUE)
  
  # Legend for node shapes
  present_types <- unique(nodes_vis$group)
  present_shapes <- type_shapes[names(type_shapes) %in% present_types]
  
  type_legend <- data.frame(
    label = names(present_shapes),
    shape = unname(present_shapes),
    color = rep("#BDBDBD", length(present_shapes)),
    borderWidth = 1,
    font.size = 12
  )
  
  # Legend for color scale
  color_legend <- data.frame(
    label = c(
      paste0("≤ ", round(-limit_val, 2)),
      "0 (neutral)",
      paste0("≥ ", round(limit_val, 2))
    ),
    shape = "dot",
    color = c("#0000CC", "#BDBDBD", "#C62828"),
    borderWidth = 0,
    font.size = 12
  )
  
  legend_df <- rbind(type_legend, color_legend)
  
  net <- net %>%
    visNetwork::visLegend(
      addNodes = legend_df,
      useGroups = FALSE
    )
  
  
  # Interaction options
  if (!is.null(options_by)) {
    if (options_by == "pathway" && "pathway_multi" %in% colnames(nodes_vis)) {
      net <- net %>% visNetwork::visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = TRUE,
        selectedBy = list(
          variable = "pathway",
          multiple = TRUE,
          main = "Select Pathway:"
        )
      )
    } else {
      net <- net %>% visNetwork::visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = TRUE
      )
    }
  }
  
  return(net)
}


network_enrichment_pg <- function(network, database, wp_ontology_names, ...) {
  
  enrich_res <- suppressWarnings(enrichment_analysis_pg(subnet = network, database = database, wp_ontology_names = wp_ontology_names, ...))
  
  enrich_res_df <- enrich_res$enrichment %>% 
    dplyr::mutate_at(c("P.value", "Adjusted.P.value", "Combined.Score"), as.numeric)
  
  num_clusters <- base::unique(enrich_res_df$Cluster)
  enrich_res_df <- dplyr::select(enrich_res_df, Cluster,
                                 Term, Overlap, P.value, Adjusted.P.value, 
                                 Combined.Score, Genes) %>% 
    dplyr::mutate_if(is.numeric, signif, 3)
  
  
  return(
    list(
      pathways = enrich_res_df,
      cluster = num_clusters
      
    )
  )


}

filter_wps_by_ontology <- function(pw_df, 
                                   ontology_names = c("classic metabolic pathway", "signaling pathway", "regulatory pathway")) {

  pw_df_wikipw <- pw_df %>% filter(Database == "WikiPathways_2024_Human") %>%
    separate_wider_delim(Term, delim = " WP",  names = c("Term", "ID"))  %>%
    mutate(ID = paste0("WP", ID),
           Term = tolower(trimws(Term))) %>%
    mutate(Term = sub("-", " ", Term ))
  pw_df_not_wikipw <- pw_df %>% filter(Database != "WikiPathways_2024_Human") %>%
    mutate(ID = "")
  
  query_wikipathways_sparql <- function(sparql_query) {
    endpoint <- "https://sparql.wikipathways.org/sparql"
    
    response <- GET(
      endpoint,
      query = list(
        query = sparql_query,
        format = "json"
      )
    )
    
    if (status_code(response) == 200) {
      result <- fromJSON(content(response, type = "text", encoding = "UTF-8"))
      
      # Extract the bindings into a dataframe
      bindings <- result$results$bindings
      
      if (length(bindings) > 0) {
        # Convert list of lists to dataframe
        df <- as.data.frame(lapply(bindings, function(col) {
          sapply(col, function(x) if(is.list(x)) x$value else x)
        }), stringsAsFactors = FALSE)
        
        return(df)
      }
    }
    return(NULL)
  }
  
  # SPARQL query
  sparql_query <- '
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT DISTINCT ?wpid ?pathwayLabel ?ontologyLabel
WHERE {
  ?pathway a wp:Pathway ;
           dcterms:identifier ?wpid ;
           dc:title ?pathwayLabel ;
           wp:ontologyTag ?ontologyTerm ;
           wp:organismName "Homo sapiens" .
  ?ontologyTerm rdfs:label ?ontologyLabel .
}
'

  # Execute query
  pathway_ontology_data <- query_wikipathways_sparql(sparql_query)
  pathway_ontology_data <- pathway_ontology_data %>%
    mutate(ontologyLabel.value = tolower(trimws(ontologyLabel.value)),
           pathwayLabel.value = tolower(trimws(pathwayLabel.value)),
           pathwayLabel.value = sub("-", " ", pathwayLabel.value))
    
  
  # pathway_ontology_data <- pathway_ontology_data %>% 
  #   mutate(ontologyLabel.value = trimws(ontologyLabel.value),
  #          pathwayLabel.value= trimws(pathwayLabel.value)) 
  
  pw_ont_sum <- pathway_ontology_data %>% group_by(wpid.value, pathwayLabel.value) %>%
    summarize(ontologies = paste0(ontologyLabel.value, collapse = "&"))
  # antipattern <- "disease|cancer|drug|syndrome|disorder|development"
  antipattern <- "drug|syndrome|development|disorder|infectious disease|epilepsy|chromosomal duplication|parkinson|alzheimer|charcot|pallister|bubonic plague|impotence"
  
  print(paste0("Filtering out the following pathway ontologies: ", antipattern))
  
  pw_ont_sum_in <- pw_ont_sum %>% 
    filter(!grepl(antipattern, ontologies))
  pw_ont_sum_out <- pw_ont_sum %>% 
    filter(grepl(antipattern, ontologies))
  
  pathway_ontology_data_filt <- pathway_ontology_data %>% 
    filter(wpid.value %in% pw_ont_sum_in$wpid.value)
  
  
  if (!is.null(pathway_ontology_data_filt) && nrow(pathway_ontology_data_filt) > 0) {
    
    # Filter your dataframe
    pw_df_filt <- pw_df_wikipw %>%
      filter(ID %in% trimws(pathway_ontology_data_filt$wpid.value))
    
    print(paste("Found ", nrow(pw_df_filt), "pathways"))
    
    
    return(bind_rows(pw_df_filt, pw_df_not_wikipw))
  } else {
    warning("No ontology match found, returning original df")
    return(pw_df)
  }
  
}


do_network_enrichment <- function(network, pval = 0.05, perc_cutoff, folder, condition = NULL,
                                  database, wp_ontology_names){
  enrichr_res <- network_enrichment_pg(network = network, database = database, wp_ontology_names = wp_ontology_names)
  enrichr_res_pw <- enrichr_res$pathways
  
  enrichr_res_pw[,"P.value"] = signif(as.numeric(as.character(enrichr_res_pw[,"P.value"])), 3)
  enrichr_res_pw[,"Adjusted.P.value"] = signif(as.numeric(as.character(enrichr_res_pw[,"Adjusted.P.value"])), 3)
  
  # Filter pathways: require minimum 3 overlapping genes across all clusters
  # When a pathway appears in multiple clusters, sum the overlapping genes
  enrichr_res_pw_filt <- enrichr_res_pw %>% 
    group_by(Term) %>%
    summarize(
      n_hits = sum(as.numeric(str_extract(Overlap, "^[0-9]+")), na.rm = TRUE),
      pathway_size = unique(as.numeric(str_extract(Overlap, "(?<=/)[0-9]+")))[1],
      overlap = n_hits / pathway_size,
      n_clusters = n(),
      Genes = paste0(Genes, collapse = ";"),
      Adjusted.P.value.mean = mean(Adjusted.P.value),
      Combined.Score.mean = mean(Combined.Score),
      Clusters = paste0(Cluster, collapse = ", "),
      .groups = 'drop'
    ) %>%
    filter(n_hits > 2) %>%  # Remove pathways with ≤2 overlapping genes
    dplyr::rename(Pathway = Term)
  
  #filter for the top 75 percentile based on the Combined.Score column
  if ("Combined.Score" %in% colnames(enrichr_res_pw)) {
    combined_score_threshold <- quantile(as.numeric(enrichr_res_pw_filt$Combined.Score.mean), 0.25, na.rm = TRUE)
    enrichr_res_pw_filt <- enrichr_res_pw_filt %>% dplyr::filter(as.numeric(Combined.Score.mean) >= combined_score_threshold)
  }
  
  if ("Reactome_Pathways_2024" %in% database){
    enrichr_res_pw_filt <- add_reactome_hierarchy(pw_df = enrichr_res_pw_filt, rpws = rpws, rel = rpwrel)
  }
  write_csv(enrichr_res_pw_filt,
            paste0(folder, "/pathways_", condition, "_p", sub(".*0\\.", "", perc_cutoff), "_all.csv"))
  
  # For each Genes value, choose one pathway based on hierarchy_n

  enrichr_res_pw_filt_short <- enrichr_res_pw_filt %>%
    group_by(Genes) %>%
    slice_min(order_by = hierarchy_n, with_ties = FALSE) %>%
    ungroup() 
  write_csv(enrichr_res_pw_filt_short,
            paste0(folder, "/pathways_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv"))
  
  return(enrichr_res_pw_filt_short)
  
  # browser()
  # enrichr_res_pw_clusters <- enrichr_res_pw %>% filter(Term %in% enrichr_res_pw_filt$Pathway)
  # 
  # # # filter enrich_res_pw with p-value
  # # enrichr_res_pw_proc <- enrichr_res_pw %>% filter(P.value < pval) %>%
  # #   group_by(Cluster, Overlap, P.value, Adjusted.P.value, Combined.Score, Genes) %>% 
  # #   summarize(pathways = paste(Term, collapse = ";"))
  # # # browser()
  # 
  # 
  # # Make shorter version
  # enrichr_res_pw_short <- enrichr_res_pw_clusters %>% 
  #   group_by(Cluster, Genes) %>% 
  #   summarize(pathways = paste(Term, collapse = " | "))
  # 
  # # the following lines result in keeping the longest subsets 
  # remove_subset_rows <- function(df, column, delimiter = ";") {
  #   # Split the strings into lists of words
  #   word_lists <- strsplit(df[[column]], delimiter)
  #   
  #   # Function to check if a list of words is a subset of another list
  #   is_subset <- function(words, other_words) {
  #     all(words %in% other_words)
  #   }
  #   
  #   # Create a logical vector to mark rows to keep
  #   keep <- sapply(seq_along(word_lists), function(i) {
  #     # Compare the current row with all other rows
  #     !any(sapply(seq_along(word_lists), function(j) {
  #       i != j && is_subset(word_lists[[i]], word_lists[[j]])
  #     }))
  #   })
  #   
  #   # Return the filtered dataframe
  #   df[keep, ]
  # }
  # 
  # # Apply the function
  # enrichr_res_pw_cleaned <- remove_subset_rows(enrichr_res_pw_short, "Genes")
  # 
  # write_csv(enrichr_res_pw_cleaned, 
  #           paste0(folder, "/enrich_results_short_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv"))
  # 
  # pw_summary <- enrichr_res_pw_cleaned %>% separate_longer_delim(pathways, " | ") %>%
  #   group_by(pathways) %>% summarize(n_clusters = n(),
  #                                    clusters = paste0(Genes, collapse = " | "))
  # 
  # write_csv(pw_summary, 
  #           paste0(folder, "/enrich_results_per_pathway_", condition, "_p", sub(".*0\\.", "", perc_cutoff), ".csv"))
  # return(pw_summary)
 
  
  
  
}



id_of <- function(name, rpws) {
  match_idx <- match(gsub("/", " ", tolower(name)), gsub("/", " ", tolower(rpws$pathway_name)))
  rpws$pathway_id[match_idx]
}

get_ancestors <- function(stId, rel, rpws) {
  parents <- rel$parent[rel$child == stId]
  if (length(parents) == 0) return(character(0))
  parent_names <- rpws$pathway_name[match(parents, rpws$pathway_id)]
  # recursively include ancestors of each parent
  unique(c(parent_names, unlist(lapply(parents, get_ancestors, rel = rel, rpws = rpws))))
}

add_reactome_hierarchy <- function(pw_df, rpws, rel) {
  pw_df$pathway_id <- sapply(pw_df$Pathway, function(name) {
    st <- id_of(name, rpws)
    if (is.na(st)) return(NA_character_)
    st
  })
  
  pw_df$hierarchy_ids <- sapply(pw_df$Pathway, function(name) {
    st <- id_of(name, rpws)
    if (is.na(st)) return(NA)
    paste(get_ancestors(st, rel, rpws), collapse = ";")
  })
  
  pw_df <- pw_df %>%
    mutate(hierarchy_ids = ifelse(hierarchy_ids == "", NA_character_, hierarchy_ids),
           hierarchy_ids = coalesce(hierarchy_ids, Pathway))

  pw_df$group <- sub("^.*;", "", pw_df$hierarchy_ids)

  # get top level hierarchy to exclude from results
  top_level_ids <- setdiff(unique(rel$parent), unique(rel$child))
  top_level_names <- rpws$pathway_name[match(top_level_ids, rpws$pathway_id)]
  top_level_names <- top_level_names[!is.na(top_level_names)]

  pw_df <- pw_df %>%
    filter(!Pathway %in% top_level_names) 

  # redundant pws have the same Genes and upstream pw in the hierarchy. 
  # Delete redundant pws by grouping with their shared first upstream pw in the hierarchy
  pw_df_n_grouped <- pw_df %>%
    mutate(hierarchy_ids = paste0(hierarchy_ids, ";"),
      Immediate_Upstream = str_extract(hierarchy_ids, "^[^;]+")) %>%
    group_by(Clusters, Genes, Immediate_Upstream) %>%
    # Count the number of pathways that share this specific Upstream parent, cluster, and Genes
    mutate(N_Redundant = n()) %>%
    ungroup() %>%
    mutate(
      Final_Pathway = if_else(
        N_Redundant > 1, 
        Immediate_Upstream, # redundant pws are collapsed to their immediate upstream pw name
        Pathway # use the original Pathway name.
      )
    ) %>%
    relocate(Final_Pathway) %>%
    distinct(Clusters, Genes, Final_Pathway, .keep_all = T) %>%
    select(-Pathway, -N_Redundant, -Immediate_Upstream) %>%
    mutate(Group_Pw = ifelse(
      !is.na(pathway_id),
      paste0(group, "_", Final_Pathway),
      Final_Pathway
    )) %>%
    relocate(Group_Pw, Final_Pathway) %>% rename("Pathway" = "Final_Pathway")

  first_hierarchy_v <- sub(";.*", "\\1", pw_df_n_grouped$hierarchy_ids)
  pw_v <- pw_df_n_grouped$Pathway
  
  # Collapse rows where Pathway matches first member of hierarchy_ids in another row and Genes also match ---
  # For each row, get first member of hierarchy_ids
  pw_df_n_grouped1 <- pw_df_n_grouped %>%
    mutate(
      first_hierarchy = sub(";.*", "\\1", hierarchy_ids),
      is_upstream = if_else(
        Pathway %in% first_hierarchy_v,
        1, 0
      ),
      is_redundant = if_else(
        first_hierarchy %in% pw_v,
        1, 0
      )
    ) %>%
    relocate(Pathway, first_hierarchy, hierarchy_ids, Genes, is_upstream, is_redundant) %>%
    arrange(Genes, desc(is_upstream))
 
  pw_df_n_grouped2 <- pw_df_n_grouped1 %>%
    filter(is_redundant != 1) %>%
    select(-first_hierarchy, -is_upstream, -is_redundant) %>%
    relocate(Group_Pw, Pathway, Genes, n_hits, pathway_size, overlap, Adjusted.P.value.mean, Combined.Score.mean, Clusters)
  pw_df_n_grouped2$hierarchy_n = lengths(strsplit(pw_df_n_grouped2$hierarchy_ids, ";"))
  pw_df_n_grouped2 <- pw_df_n_grouped2[,c("Group_Pw", "Pathway", "Genes", "hierarchy_n", "n_hits", "pathway_size", "overlap",  "Adjusted.P.value.mean", 
                                          "Combined.Score.mean", "Clusters", "hierarchy_ids", "pathway_id", "group")]
  pw_df_n_grouped2 %>% filter(hierarchy_n > 1)
  
  
  
  # reactome pw hierarchy network
  # reactome_pw_hierarchy_vis(pw_df = pw_df_n_grouped2, result_path = "results")
  
  
}


reactome_pw_hierarchy_vis <- function(pw_df, result_path){
  wrap_long_names <- function(name) {
    # Inserts a newline character \n after the last space found before the 25th character
    # This tries to break the line on whole words.
    str_replace_all(name, "(.{1,25})(\\s)", "\\1\n")
  }
  
  network_edges <- pw_df %>%
    filter(group == "Signal Transduction") %>%
    mutate(H = paste0(Pathway, ";", hierarchy_ids)) %>%
    mutate(H = sub(";$", "", H)) %>%
    select(H) %>% 
    mutate(Hierarchy_list = str_split(H, ";")) %>%
    # Unnest the list to create one row for every pathway name in the chain
    # This makes it easier to compare adjacent nodes
    unnest(Hierarchy_list) %>%
    group_by(H) %>%
    mutate(
      to = Hierarchy_list,
      from = lead(Hierarchy_list)
    ) %>%
    ungroup() %>%
    filter(!is.na(from)) %>%
    select(from, to) %>%
    distinct()
  
  network_nodes <- data.frame(id = union(network_edges$from, network_edges$to)) 
  
  my_graph <- graph_from_data_frame(network_edges, directed = TRUE)
  
  V(my_graph)$name <- sapply(V(my_graph)$name, wrap_long_names)
  
  p <- ggraph(my_graph, layout = 'tree') + 
    geom_edge_fan(
      arrow = arrow(length = unit(3, 'mm')), 
      end_cap = circle(5, 'mm') # Adjust this value if your boxes are too big/small
    ) + 
    geom_node_label(
      aes(label = name), 
      fill = "#ADD8E6", # Light blue background for the box
      color = "black",  # Text color
      size = 4,         # Font size
      alpha = 0.9,      # Slightly transparent box
      label.padding = unit(0.3, "lines"), # Padding inside the box
      label.r = unit(0.2, "lines")      # Slightly rounded corners
    ) + 
    theme_graph()
  browser()
  ggsave(paste0(result_path, "/Signal_Transduction_reactome_pw_tree.png"), p, w = 30, h = 20, unit = "cm")
  
}


enrichment_analysis_pg <-function(subnet, mode=NULL, gene_universe, database, wp_ontology_names){
  # Checking function arguments
  if (missing(subnet))
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  # if (class(subnet)[1] != "PCSF" || class(subnet)[2] != "igraph")
  #   stop("The subnetwork must be a \"PCSF\" object derived from an \"igraph\" class.")
  if (!is.null(mode)){
    if(mode==1 && missing(gene_universe))
      stop("Need to specify a list of genes (vector of gene symbols) used as background in enrichment analysis by topGO package")
  }
  
  cat("  Performing enrichment analysis...\n\n")
  
  # Obtain clusters in the subnet using edge betweenness clustering algorithm from igraph package.
  clusters = cluster_edge_betweenness(subnet)
  
  # Perform ebrichment analysis for each cluster using EnrichR through its API or topGO.
  
  havingInternet <- function() {
    if (.Platform$OS.type == "windows") {
      ipmessage <- system("ipconfig", intern = TRUE)
    } else {
      ipmessage <- system("ifconfig", intern = TRUE)
    }
    validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
    any(grep(validIP, ipmessage))
  }
  
  internet_connection <- havingInternet()
  
  if(!is.null(mode)){
    if(mode==0){
      if(internet_connection){
        cat("  Enrichment is being performed by EnrichR (http://amp.pharm.mssm.edu/Enrichr) API ...\n")
        enrich = call_enr_pg(clusters = clusters, mode = 0, gene_universe = gene_universe, database = database)
      }
      else{
        stop("There is no working Internet connection, perform your enrichment with topGO package with mode=1 by providing background gene list ...\n")
      }
    }
    else{
      cat("  Enrichment is being performed by topGO package ...\n")
      enrich = call_enr_pg(clusters = clusters, mode = mode, gene_universe = gene_universe)
    }
  } else {
    if(internet_connection){
      cat("  Enrichment is being performed by EnrichR (http://amp.pharm.mssm.edu/Enrichr) API ...\n")
      enrich = call_enr_pg(clusters = clusters, mode = 0, gene_universe = gene_universe, database = database)
    }
    else{
      stop("There is no working Internet connection, perform your enrichment with topGO package with mode=1 by providing background gene list ...\n")
    }
  }
  
  if('Compound'%in% V(subnet)$type){##then we have drugs!
    require(dplyr)
    comps=data.frame(Drug=V(subnet)$name[which(V(subnet)$type=='Compound')],
                     Cluster=clusters$membership[which(V(subnet)$type=='Compound')])%>%
      dplyr::group_by(Cluster)%>%
      dplyr::summarise(DrugsByBetweenness=paste(Drug,collapse=';'))
    
  }
  else{
    comps <-NULL
  }
  enrichment = enrich[[1]]
  enrichment_complete = enrich[[2]]
  
  novals<-which(unlist(sapply(enrich[[2]],function(x) is.null(dim(x)))))
  if(length(novals)>0)
    enrichment_complete <- enrichment_complete[-novals]
  enrichment_tab = do.call(rbind,lapply(c(1:length(enrichment_complete)),function(x) data.frame(Cluster=x,enrichment_complete[[x]])))
  
  if ("WikiPathways_2024_Human" %in% database) {
    enrichment_tab <- filter_wps_by_ontology(pw_df = enrichment_tab, ontology_names = wp_ontology_names)
  }

  more.than.two=which(sapply(enrichment_tab$Genes,function(x) length(unlist(strsplit(x,split=';')))>2))
  if(length(more.than.two)>0)
    enrichment_tab=enrichment_tab[more.than.two,]
  if(!is.null(comps))
    enrichment_tab = enrichment_tab%>%dplyr::left_join(comps,by='Cluster')
  
  # Add 'group" and 'title' attributes to subnet
  V(subnet)$group = clusters$membership
  V(subnet)$title = paste0("Cluster ",clusters$membership,": Enrichment analysis")
  for( i in 1:length(V(subnet))){
    V(subnet)$title[i] = paste0( V(subnet)$title[i], enrichment[[V(subnet)$group[i]]])
  }
  
  # Derive a "PCSFe" object from an "igraph" class.
  class(subnet) <- c("PCSFe", "igraph")
  # Combine the subnetwork and colplete enrichment analysis tables.
  
  output = list(subnet, enrichment_tab)
  names(output) = c("subnet", "enrichment")
  
  return (output)
}

call_enr_pg <- function(clusters, mode=0, gene_universe, 
                        database = c("KEGG_2021_HUMAN", "Reactome_Pathways_2024", 
                                     "WikiPathways_2024_Human")){
  
  if( mode == 0){
    
    # Enrichment analysis is performed with ENRICHR
    ENRICHR_ADDLIST = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_EXPORT = 'http://amp.pharm.mssm.edu/Enrichr/export'
    
    # Enrichment results
    enrichment_result = as.list(1:length(clusters))
    enrichment_result_complete = as.list(1:length(clusters))

    # Perform Enrichment Analysis for each cluster in the forest
    for( a in 1:length(clusters)){
      
      # List of genes to be regusted for enrichment via ENRICHR API
      genes = clusters[[a]]
      request = list(list = paste(genes, collapse = "\n"))
      complete_request = httr::POST(ENRICHR_ADDLIST, body = request)
      output =httr::content(complete_request, type = "text", encoding = "ISO-8859-1")
      userListID  = strsplit(strsplit(output, "\n")[[1]][3], ": ")[[1]][2]
      
      response_collection=NULL
      
      # Request enrichment for each database and combine them all
      for( b in 1:length(database)){
        
        # Gather an EXPORT URL and the Response
        url = paste0(ENRICHR_EXPORT, "?userListId=",userListID, "&backgroundType=", database[b])
        response = GET(url)
        response = httr::content(response, type = "text",  encoding = "ISO-8859-1")
        response = strsplit(response, "\n")[[1]]
        #  if(length(unlist(response  ))==1)
        #       next
        response = lapply(response, function(x){sp = strsplit(x, "\t")[[1]]; return (sp)})
        
        # If the response contains some elements then combine it
        # Defensive: skip if response is header only or malformed
        if(length(response) > 1 && length(response[[1]]) > 0) {
          x = length(response) - 1
          m_resp = as.data.frame(matrix(0, nrow = x, ncol = length(response[[1]])))
          colnames(m_resp) = response[[1]]
          for(i in 1:x){
            # Defensive: only assign if the row has the right length
            if(length(response[[i+1]]) == length(response[[1]])) {
              m_resp[i,] = response[[i+1]]
            }
          }
          m_resp$Database = database[b]
          response_collection = rbind(response_collection, m_resp)
        }
      }
      
      if(is.null(response_collection))
        next
      # Reorder the enrichment according to the "Adjusted P-value"
      ordered_resp = data.frame(response_collection$`Term`, response_collection$`Adjusted P-value`, response_collection$`Combined Score`, response_collection$Database)
      ordered_resp = ordered_resp[order(ordered_resp[,2]),]
      ordered_resp[,2] = signif(as.numeric(as.character(ordered_resp[,2])), 3)
      ordered_resp[,3] = signif(as.numeric(as.character(ordered_resp[,3])), 3)
      # Convert the enrichment table into HTML format in order to display it
      enrich = "<!DOCTYPE html> <html> <head> <style>
      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,
      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}
      tr:nth-child(even) {background-color: #dddddd;}
      </style> </head> <body>
      <table> <tr>  <th>Term</th> <th>Adjusted P-value</th> <th>Combined Score</th> <th>Database</th> </tr>";
      for(i in 1:nrow(ordered_resp)){
        enrich = paste0(enrich, " <tr>")
        for(j in 1:ncol(ordered_resp)){
          enrich = paste0(enrich, "<td>",ordered_resp[i,j], "</td>")
        }
        enrich = paste0(enrich, "</tr> ")
      }
      
      enrich = paste0(enrich, "</table> </body> </html>")
      
      # Attach the Enrichment Analysis for the current cluster

      enrichment_result[[a]] = enrich
      enrichment_result_complete[[a]] = response_collection
    }
    
  } else{
    
    # Enrichment analysis is performed by topGO
    
    # Enrichment results
    enrichment_result = as.list(1:length(clusters))
    enrichment_result_complete = as.list(1:length(clusters))
    
    # Perform Enrichment Analysis for each cluster in the forest
    for( a in 1:length(clusters)){
      
      # List of genes to be regusted for enrichment by topGO
      genes = clusters[[a]]
      fg <- factor(as.integer(gene_universe %in% genes))
      names(fg) <- gene_universe
      tgData <- new("topGOdata", description = "simple_session", ontology = "BP",
                    allGenes=fg, nodeSize=15, annot=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
      
      resultFisher <- runTest(tgData, algorithm = "classic", statistic = "fisher")
      resultKS <- runTest(tgData, algorithm = "classic", statistic = "ks")
      
      res_table_top15 <- GenTable(tgData, classicFisher = resultFisher,
                                  classicKS = resultKS,
                                  orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 15)
      res_table_top1000 <- GenTable(tgData, classicFisher = resultFisher,
                                    classicKS = resultKS,
                                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000)
      
      
      # Reorder the enrichment according to the "Adjusted P-value" and select the top 15 enrichments
      #ordered_resp[,2] = signif(as.numeric(as.character(ordered_resp[,2])), 3)
      #ordered_resp[,3] = signif(as.numeric(as.character(ordered_resp[,3])), 3)
      
      # Convert the enrichment table into HTML format in order to display it
      enrich = "<!DOCTYPE html> <html> <head> <style>
      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,
      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}
      tr:nth-child(even) {background-color: #dddddd;}
      </style> </head> <body>
      <table> <tr> <th>GO.ID</th> <th>Term</th> <th>Annotated</th> <th>Significant</th> <th>Expected</th> <th>classicFisher</th> <th>classicKS</th> </tr>";
      for(i in 1:nrow(res_table_top15)){
        enrich = paste0(enrich, " <tr>")
        for(j in 1:ncol(res_table_top15)){
          enrich = paste0(enrich, "<td>",res_table_top15[i,j], "</td>")
        }
        enrich = paste0(enrich, "</tr> ")
      }
      enrich = paste0(enrich, "</table> </body> </html>")
      
      # Attach the Enrichment Analysis for the current cluster
      enrichment_result[[a]] = enrich
      enrichment_result_complete[[a]] = res_table_top1000
    }
    
  }
  
  
  return (list(enrichment_result, enrichment_result_complete))
  
}


plot_cell_networks_kinase <- function(uka, art_nodes, art_lfc, spec_cutoff, respath, perc_cutoffs,
                               rank_uka_abs = T, ppi_network = ppi_networkv12, b, cs,
                               wp_ontology_names = c("metabolic pathway", "signaling pathway", "signaling", "regulatory pathway"), 
                               highlight_degree) {
  # Ensure base output directory exists
  if (!dir.exists(respath)) {
    dir.create(respath, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Clean up any existing temporary files to avoid conflicts
  temp_pattern <- file.path(respath, "*_files")
  existing_temp_dirs <- list.dirs(respath, recursive = FALSE)[grepl("_files$", list.dirs(respath, recursive = FALSE))]
  if (length(existing_temp_dirs) > 0) {
    sapply(existing_temp_dirs, function(d) {
      if (dir.exists(d)) {
        tryCatch(unlink(d, recursive = TRUE, force = TRUE), error = function(e) {})
      }
    })
  }
  
  # Parse input data, filter for comparisons vs control
  uka_parsed <- clean_uka_to_kinograte(uka, cs = cs)
  conditions <- unique(uka_parsed$Sgroup_contrast)
  
  for (perc_cutoff in perc_cutoffs) {
    # Sequential loop over cells for this percentile cutoff
    temp_files <- c()
    for (condition in conditions) {
      print(paste0("Constructing network for perc cutoff: ", perc_cutoff, " and condition: ", condition, "..."))
      # Take top hits for uka for observed score
      uka_filt <- uka_parsed %>%
        filter(Sgroup_contrast == condition) %>%
        uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff, cs = cs)
      make_network_and_stats(
        uka = uka_filt, art_nodes = art_nodes, art_lfc = art_lfc, perc_cutoff = perc_cutoff,
        spec_cutoff = spec_cutoff, res.path = respath, condition = condition,
        write = T, ppi_network = ppi_network,
        b = b, wp_ontology_names = wp_ontology_names, highlight_degree = highlight_degree
      )
    }
  }
  
}


plot_cell_networks <- function(uka, sens, del_cells = NULL, 
                               control, spec_cutoff, zscore, best_drug_per_target = NULL, respath, perc_cutoffs, 
                               rank_uka_abs = T, balance = F, ppi_network = ppi_networkv12, only_kinase = F, b, cs,
                               wp_ontology_names, highlight_degree) {
  
  # Ensure base output directory exists
  if (!dir.exists(respath)) {
    dir.create(respath, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Clean up any existing temporary files to avoid conflicts
  temp_pattern <- file.path(respath, "*_files")
  existing_temp_dirs <- list.dirs(respath, recursive = FALSE)[grepl("_files$", list.dirs(respath, recursive = FALSE))]
  if (length(existing_temp_dirs) > 0) {
    sapply(existing_temp_dirs, function(d) {
      if (dir.exists(d)) {
        tryCatch(unlink(d, recursive = TRUE, force = TRUE), error = function(e) {})
      }
    })
  }
  
  # Parse input data, filter for comparisons vs control
  sens_parsed <- clean_sens_to_kinograte(sens, control = control, zscore = zscore,
                                         best_drug_per_target = best_drug_per_target)
  
  uka_parsed <- clean_uka_to_kinograte(uka, control = control, cs = cs)
  # Find common cell lines in uka and sens
  common_cells <- intersect(unique(uka_parsed$cell_line), unique(sens_parsed$cell_line))
  if (!is.null(del_cells)){
    common_cells <- common_cells[!common_cells %in% del_cells]   
  }
  # Sequential loop over percentile cutoffs
  loop_start <- Sys.time()
  logs <- character() # Collect log messages
  for (perc_cutoff in perc_cutoffs) {
    # Sequential loop over cells for this percentile cutoff
    temp_files <- c()
    for (cell in common_cells) {
      # Take top hits for uka and sens for observed score
      uka_filt <- uka_parsed %>% filter(cell_line == cell) %>% uka_top(spec_cutoff = spec_cutoff, rank_uka_abs = rank_uka_abs, perc_cutoff = perc_cutoff)
      sens_filt <- sens_parsed %>% filter(cell_line == cell) %>% sens_top(perc_cutoff, balance)
      # For random sampling, filter for cell
      uka_cell_all <- uka_parsed %>% filter(cell_line == cell)
      
      obs_res <- make_network_and_stats(
        uka = uka_filt, sens = sens_filt,
        perc_cutoff = perc_cutoff, spec_cutoff = spec_cutoff,
        res.path = respath, condition = cell, write = T,
        ppi_network = ppi_network, relative_to = "n_nodes", only_kinase = only_kinase, b = b, 
        wp_ontology_names = wp_ontology_names, highlight_degree = highlight_degree
      )
    }
  }
}

