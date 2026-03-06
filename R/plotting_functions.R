library(magrittr)
library(ggalt)
library(ggfortify)
library(ggtext)
# library(ggVennDiagram)
library(ComplexHeatmap)
library(RColorBrewer)

plot_var_explained <- function (object, x = "view", y = "factor", split_by = NA, plot_total = FALSE, 
          factors = "all", min_r2 = 0, max_r2 = NULL, legend = TRUE, 
          use_cache = TRUE, ...) 
{

  if (length(unique(c(x, y, split_by))) != 3) {
    stop(paste0("Please ensure x, y, and split_by arguments are different.\n", 
                "  Possible values are `view`, `group`, and `factor`."))
  }
  if (is.na(split_by)) 
    split_by <- setdiff(c("view", "factor", "group"), c(x, 
                                                        y, split_by))
  if ((use_cache) & .hasSlot(object, "cache") && ("variance_explained" %in% 
                                                  names(object@cache))) {
    r2_list <- object@cache$variance_explained
  }
  else {
    r2_list <- calculate_variance_explained(object, factors = factors, 
                                            ...)
  }
  r2_mk <- r2_list$r2_per_factor
  r2_mk_df <- melt(lapply(r2_mk, function(x) melt(as.matrix(x), 
                                                  varnames = c("factor", "view"))), id.vars = c("factor", 
                                                                                                "view", "value"))
  colnames(r2_mk_df)[ncol(r2_mk_df)] <- "group"
  if ((length(factors) == 1) && (factors[1] == "all")) {
    factors <- factors_names(object)
  }
  else {
    if (is.numeric(factors)) {
      factors <- factors_names(object)[factors]
    }
    else {
      stopifnot(all(factors %in% factors_names(object)))
    }
    r2_mk_df <- r2_mk_df[r2_mk_df$factor %in% factors, ]
  }
  r2_mk_df$factor <- factor(r2_mk_df$factor, levels = factors)
  r2_mk_df$group <- factor(r2_mk_df$group, levels = groups_names(object))
  r2_mk_df$view <- factor(r2_mk_df$view, levels = views_names(object))
  groups <- names(r2_list$r2_total)
  views <- colnames(r2_list$r2_per_factor[[1]])
  if (!is.null(min_r2)) 
    r2_mk_df$value[r2_mk_df$value < min_r2] <- 0.001
  min_r2 = 0
  if (!is.null(max_r2)) {
    r2_mk_df$value[r2_mk_df$value > max_r2] <- max_r2
  }
  else {
    max_r2 = max(r2_mk_df$value)
  }
  
  title <- round(head(model@cache$variance_explained$r2_total[[1]]), 1)
  title <- paste0(paste0(c(names(object@data[1]), names(object@data[2])), c(": ", ": "), 
                         as.character(title), c("%", "%")), collapse = ', ')
  print(title)
  
  p1 <- ggplot(r2_mk_df, aes(x = .data[[x]], y = .data[[y]])) + 
    geom_tile(aes(fill = .data$value), color = "black") + 
    geom_text(aes(.data[[x]], .data[[y]], label=round(.data$value, 1)), colour = "red", check_overlap = TRUE) +
    facet_wrap(as.formula(sprintf("~%s", split_by)), nrow = 1) + 
    labs(x = "", y = "", title = "") + scale_fill_gradientn(colors = c("gray97", 
                                                                       "darkblue"), guide = "colorbar", limits = c(min_r2, max_r2)) + 
    guides(fill = guide_colorbar("Var. (%)")) + theme(axis.text.x = element_text(size = rel(1), 
                                                                                 color = "black"), axis.text.y = element_text(size = rel(1.1), 
                                                                                                                              color = "black"), axis.line = element_blank(), axis.ticks = element_blank(), 
                                                      panel.background = element_blank(), strip.background = element_blank(), 
                                                      strip.text = element_text(size = rel(1))) +
    ggtitle(title)
  
  if (isFALSE(legend)) 
    p1 <- p1 + theme(legend.position = "none")
  if (length(unique(r2_mk_df[, split_by])) == 1) 
    p1 <- p1 + theme(strip.text = element_blank())
  if (plot_total) {
    r2_m_df <- melt(lapply(r2_list$r2_total, function(x) lapply(x, 
                                                                function(z) z)), varnames = c("view", "group"), value.name = "R2")
    colnames(r2_m_df)[(ncol(r2_m_df) - 1):ncol(r2_m_df)] <- c("view", 
                                                              "group")
    r2_m_df$group <- factor(r2_m_df$group, levels = MOFA2::groups_names(object))
    r2_m_df$view <- factor(r2_m_df$view, levels = views_names(object))
    min_lim_bplt <- min(0, r2_m_df$R2)
    max_lim_bplt <- max(r2_m_df$R2)
    p2 <- ggplot(r2_m_df, aes(x = .data[[x]], y = .data$R2)) + 
      geom_bar(stat = "identity", fill = "deepskyblue4", 
               color = "black", width = 0.9) + facet_wrap(as.formula(sprintf("~%s", 
                                                                             split_by)), nrow = 1) + xlab("") + ylab("Variance explained (%)") + 
      scale_y_continuous(limits = c(min_lim_bplt, max_lim_bplt), 
                         expand = c(0.005, 0.005)) + theme(axis.ticks.x = element_blank(), 
                                                           axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"), 
                                                           axis.title.y = element_text(color = "black"), axis.line = element_line(color = "black"), 
                                                           panel.background = element_blank(), strip.background = element_blank(), 
                                                           strip.text = element_text())
    if (length(unique(r2_m_df[, split_by])) == 1) 
      p2 <- p2 + theme(strip.text = element_blank())
    plot_list <- list(p1, p2)
  }
  else {
    plot_list <- p1
  }
  return(plot_list)
}



make_pathway_pcas <- function(pw_list, bestf_list, idmap, w = 2500, kegg.gs, e, mycolor, myshape, axis1 = 1, axis2 = 2){ #f
  for (i in 1:length(pw_list)){
    for (j in 1:3){
      pw <- kegg.gs[grepl(pw_list[[i]], names(kegg.gs))]	
      name <- paste0(names(bestf_list[j]), " - ", names(pw_list[i]), " pathway")
      print(name)
      mydf <- bestf_list[[j]]  %>%
        filter(type == "pep") %>% 
        left_join(idmap, by = c("ID" = "anonID")) %>%
        filter(!is.na(entrezid), entrezid %in% pw[[1]]) %>%
        select(matches("inspire")) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Patientcode")
      
      p <- plot_pca(df_topca = mydf, e = e, mycolor, myshape, mytitle = name, axis1 = axis1, axis2 = axis2)
      
      ggsave(paste0("plots/pca/pca_per_pathway/", name, ".png"), p, width = 11, height = 8, units = "cm")
      

    }
  }
  
}

filter_data_for_pw <- function(fc_df, idmap, pw){
  fc_df$entrezid <- as.character(fc_df$entrezid)
  idmap$entrezid <- as.character(idmap$entrezid)
  fc_df %>%
    filter(entrezid %in% pw[[1]]) %>% 
    left_join(idmap, by = "entrezid")  %>% select(-entrezid)

}


parse_dfs_to_heatmapformat <- function(prot_df, pep_df, idmap, pw){
  prot_df_pw <- prot_df %>%
    filter(!is.na(ID), ID %in% pw[[1]]) %>% 
    left_join(idmap, by = c("ID" = "entrezid")) %>% mutate(type = "Pr")
  pep_df_pw <- pep_df %>%
    filter(!is.na(ID), ID %in% pw[[1]]) %>% 
    left_join(idmap, by = c("ID" = "entrezid")) %>% mutate(type = "")
  
  alldf <- rbind(prot_df_pw, pep_df_pw) %>% select(-ID)
  return(alldf)
  
}


# make_heatmaps <- function(heatmap_mat, name, ha, idmap, w = 2500, e, rows_clustered = T, 
#                           cols_clustered = F, coloring = "normal", row_annotation = F, res.path){ 
#   
#   
#   
#   plot_heatmap(df = heatmap_mat, ha = ha, name = name, e = e, w = w, coloring = coloring, 
#                rows_clustered = rows_clustered, cols_clustered = cols_clustered, row_annotation = row_annotation,
#                res.path = res.path)
# }

make_pathway_tilemaps <- function(pw_list, fc_df, idmap, manual_size = F, w = NULL, h = NULL, kegg.gs, rows_clustered = T, 
                          cols_clustered = F, coloring = "normal", res.path = "results/pathway_tilemaps/"){ 
  # Makes an FC tilemap for each pathway.
  # Input:
  # coloring: "normal" or "limited_scale"
  # fc_df: columns should be: comparison (values e.g. "T vs C"; arranged factor), entrezid, data_from (values e.g. "Pep", 'Ki', 'MS'; arranged factor), and value (LFC).
  # idmap: entrezid to uniprot name, contains all data
  
  
  
  if (!dir.exists(res.path)){
    dir.create(res.path, recursive = T)
  }
  
  # for each pathway:
  for (i in 1:length(pw_list)){
    # get pw elements

    pw <- kegg.gs[grepl(pw_list[[i]], names(kegg.gs))]	
    name <- paste0(names(pw_list[i]))
    print(name)
    #browser()
    pws_not_found <- c()
    if (length(pw) == 0) {
      print(paste("Pathway not found: ", name))
      pws_not_found <- c(pws_not_found, name)
      next  # Skip to the next iteration
    } else {
      # filter data for pw elements
      mymat <- filter_data_for_pw(fc_df = fc_df, idmap = idmap, pw = pw)
      plot_pathway_tilemap(df = mymat, name = name, manual_size = F,  w = NULL, h = NULL, coloring = coloring, 
                           rows_clustered = rows_clustered, cols_clustered = cols_clustered,
                           res.path = res.path)
    }
    print("Pathways not found")
    print(pws_not_found)
  }
}



parse_bestf_to_pcaformat <- function(bestf, idmap, pw){ 
  bestf %>%
    filter(type == "pep") %>% 
    left_join(idmap, by = c("ID" = "anonID")) %>%
    filter(!is.na(entrezid), entrezid %in% pw[[1]]) %>%
    select(-ID) %>%
    rename("ID" = "ID.y") %>%
    select(-FactorNum, -type, -coef, -PepProtein_UniprotID, -PepProtein_SeqSimilarity, -Sequence, -entrezid, -PepProtein_UniprotName) %>%
    distinct() %>%
    arrange(ID) %>%
    column_to_rownames("ID")%>%
    as.matrix()
}



plot_heatmap <- function(df, ha = NULL, value_limits = NULL, name,legend_name, w = 3000, h = 2000, rows_clustered = F, cols_clustered = F, 
                         show_row_names = F, coloring = "normal", res.path){
  # df: matrix with optional row and colnames
  # ha: heatmap annotation df with first column containing the colnames of X, 2nd col is annotation
  # value_limits: if not null, clips the df according to min and max value. example: c(0,1)
  
  parse_ha <- function(ha, X){
    # parses the ha input to a vector with names of the same colnames as X, values as the annotation. 
    # Orders the vector in the same way as X
    first_col <- colnames(ha)[1]
    second_col <- colnames(ha)[2]
    col_names <- colnames(X)
    ha_vec <- setNames(ha[[second_col]][match(col_names, ha[[first_col]])], col_names)
    return(ha_vec)
  }
  
  
  
  if(!is.null(value_limits)){
    df <- as.data.frame(df) %>%
      rownames_to_column("row_id") %>%
      pivot_longer(-row_id, names_to = "col_id", values_to = "value") %>%
      mutate(value = pmax(value_limits[1], pmin(value, value_limits[2]))) %>%
      pivot_wider(names_from = col_id, values_from = value) %>%
      column_to_rownames("row_id") %>% as.matrix()
  }
  
  if (!dir.exists(res.path)){
    dir.create(res.path, recursive = T)
  }

  # coloring
  if (coloring == "diff_normal"){
    colorfunc <- circlize::colorRamp2(seq(min(df), max(df), length = 3), c("green", "black", "red"))
  } else if (coloring == 'diff_limited'){
    colorfunc <- circlize::colorRamp2(seq(-1, 1, length = 3), c("green", "black", "red"))
  } else if (coloring == "viridis"){
    colorfunc <- viridis::viridis(100)
  } else if (coloring == "onlyneg"){
    colorfunc <- circlize::colorRamp2(seq(-7, 0, length = 2), c("green", "black"))
  } else if (coloring == "onlypos"){
    colorfunc <- circlize::colorRamp2(seq(0, (max(df)), length = 2), c("blue", "red"))
  } else if (coloring == 'discrete'){
    colorfunc <- c("DOWN" = "green", "0" = "black", "UP" = "red") 
  }

  # Kinase group annotation and colors
  ha <- NULL
  if (!is.null(ha)) {
    ha_vec <- parse_ha(ha, X)
    kinase_groups_unique <- unique(ha_vec)
    kinase_group_colors <- setNames(RColorBrewer::brewer.pal(min(length(kinase_groups_unique), 12), "Paired")[seq_along(kinase_groups_unique)], kinase_groups_unique)
    ha <- ComplexHeatmap::HeatmapAnnotation(Kinase_Group = ha_vec, 
                                            col = list(Kinase_Group = kinase_group_colors))
  }

  png(filename = paste0(res.path, "/Heatmap_", name, ".png"), 
      width = w, height = h, units = "px", res = 150)

  heatmap_args <- list(
    df,
    name = legend_name, 
    cluster_columns = cols_clustered,
    cluster_rows = rows_clustered,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    column_names_side = "top",
    row_km = 0,
    column_km = 0,
    column_title_rot = 0, 
    column_title = name,
    col = colorfunc,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8)
  )

  if (!is.null(ha)) {
    heatmap_args$top_annotation <- ha
    p <- do.call(Heatmap, heatmap_args)
  } else {
    p <- do.call(Heatmap, heatmap_args)
  }

  p <- draw(p)
  dev.off()
}


plot_pathway_tilemap <- function(df, name, manual_size, w, h, rows_clustered = F, 
                                 cols_clustered = F, coloring = "normal",
                                 res.path){
  

  if (coloring == "limited_scale"){
    p25_neg <- stats::quantile(df$value[df$value < 0], probs = 0.25, na.rm = TRUE)
    p75_pos <- stats::quantile(df$value[df$value > 0], probs = 0.75, na.rm = TRUE)
    
    limit_val <- max(abs(p25_neg), abs(p75_pos), na.rm = TRUE)
    color_limits <- c(-limit_val, limit_val)
    
    df$value <- pmax(pmin(df$value, color_limits[2]), color_limits[1])
  }  
  
  

  p <- ggplot(df, aes(x = data_from, y = uniprotname, fill = value)) + 
    geom_tile() +
    facet_wrap(df$comparison, nrow = 1) + 
    labs(x = "Comparison", y = "", title = name) + 
    theme(axis.text.x = element_text(size = 15, angle = 35, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 20),
          strip.text = element_text(size = 15),
          title = element_text(size = 20)) +
    scale_fill_gradient2(low = "green",
                         mid = "black",
                         high = "red")
  
  if (manual_size == F){
    # set width
    n_data_from <- length(unique(df$data_from))
    
    if(n_data_from == 1){
      w = 30
    } else if (n_data_from == 2){
      w = 40
    } else if (n_data_from == 3){
      w = 50
    } else{
      w = 60
    }
    
    
    if (nrow(df) < 400){
      h = 60
    } else if (nrow(df) >= 400  & nrow(df) <= 500){
      h = 70
    } else if (nrow(df) >500 & nrow(df <= 750)){
      h = 80
    } else if (nrow(df) > 750){
      h = 90
    }
  }
  
  ggsave(paste0(res.path, "/", name, ".png"), p, width = w, height = h, units = "cm" )
  
}


plot_rna_pep_hist <- function(rnav, pepv, titles, data_from){
  ncol_fig <- length(rnav)
  
  jpeg(paste0("results/plots/hist_prot_pep_normalizations_", data_from, ".jpg"), width = 1200, height = 600, units = "px")
  
  par(mfrow= c(2,ncol_fig))
  #par(mar = c(2.5, 2.5, 2.5, 2.5))
  for (i in 1:ncol_fig){
    hist(rnav[[i]], 20, main = titles[i], xlab = 'Protein')
  } 
  for (j in 1:ncol_fig){
    hist(pepv[[j]], 20, main = NULL, xlab = "pep")
  }
  
  dev.off()
}




############
### MOFA ###
############





plot_factor_pca <- function(model, mycolor, myshape, f1, f2, perfactor_var, w = 14, h = 14){
  
  a <- plot_factors_pg(model, 
                    factors = c(f1, f2),
                    color_by = mycolor,
                    shape_by = myshape,
                    #color_name = "sample_type",
                    dot_size = 2.5,
                    scale = F,
                    stroke = NA) + 
    xlab(paste0("Factor ", f1, " (R2 = ", round(sum(perfactor_var[paste0('Factor', f1),]), 1), "%)")) + 
    ylab(paste0("Factor ", f2, " (R2 = ", round(sum(perfactor_var[paste0('Factor', f2),]), 1), "%)")) +
    theme_gray() +
    theme(axis.text.x= element_text(colour = 'black'),
          axis.text.y = element_text(colour = 'black')) +
    ggtitle(paste0(modelname))
  
  
  
  myfilename = paste0("./results/plots/mofamodel_", modelname, 
                      "/factorplot_color_", mycolor, 
                      "_shape_", myshape,
                      "_f", f1, "-f", f2,
                      ".png")
  ggsave(myfilename, a, width = w, height = h, units = "cm")
  a
}

plot_factors_pg <- function (object, factors = c(1, 2), groups = "all", show_missing = TRUE, 
          scale = FALSE, color_by = NULL, shape_by = NULL, color_name = NULL, 
          shape_name = NULL, dot_size = 2, alpha = 1, legend = TRUE, 
          stroke = NULL, return_data = FALSE) {
  
  
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (length(unique(factors)) == 1) {
    .args <- as.list(match.call()[-1])
    .args <- .args[names(.args) != "factors"]
    return(do.call(plot_factor, c(.args, list(factors = unique(factors)))))
  }
  else if (length(factors) > 2) {
    .args <- as.list(match.call()[-1])
    p <- do.call(.plot_multiple_factors, .args)
    return(p)
  }
  
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name)) 
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name)) 
    shape_name <- shape_by
  factors <- .check_and_get_factors(object, factors)
  Z <- get_factors(object, factors = factors, groups = groups, 
                   as.data.frame = TRUE)

  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  Z <- Z[complete.cases(Z), ]
  df <- merge(Z, color_by, by = "sample")
  df <- merge(df, shape_by, by = "sample")
  df$shape_by <- as.character(df$shape_by)
  if (isFALSE(show_missing)) 
    df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  df <- spread(df, key = "factor", value = "value")
  df <- df[, c(colnames(df)[seq_len(4)], factors)]
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "x", "y"))
  
  if (scale) {
    df$x <- df$x/max(abs(df$x))
    df$y <- df$y/max(abs(df$y))
  }
  if (return_data) 
    return(df)
  if (is.null(stroke)) {
    stroke <- .select_stroke(N = length(unique(df$sample)))
  }

  p <- ggplot(df, aes(x = .data$x, y = .data$y, color = .data$color_by, 
                      shape = .data$shape_by)) + 
    geom_point(size = dot_size, alpha = alpha) +
    #geom_encircle(aes(group = .data$shape_by, fill = .data$shape_by), alpha = 0.1, expand = 0, s_shape = 1) +
    labs(x = factors[1], y = factors[2]) + 
    #scale_shape_manual(values=c(16, 17)) + 
    scale_shape_manual(values=1:nlevels(df %>% pull(shape_by) %>% as.factor())) +
    theme(axis.text = element_text(size = rel(0.8), color = "black"), 
          axis.title = element_text(size = rel(1.1), color = "black"), 
          axis.line = element_line(color = "black", linewidth = 0.5), 
          axis.ticks = element_line(color = "black", linewidth = 0.5))
  
  #browser()
  #p <- .add_legend(p, df, legend, color_name, shape_name)
  if (!is.null(color_name)) {
    p <- p + labs(color = color_name)
  }
  if (!is.null(shape_name)) {
    p <- p + labs(shape = shape_name)
  }
  return(p)
}

.plot_multiple_factors <- function(object, factors = "all", show_missing = TRUE, dot_size = 1,
                                   color_by = NULL, color_name = "", shape_by = NULL, shape_name = "") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Collect relevant data
  Z <- get_factors(object, factors=factors, as.data.frame=TRUE)
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing factor values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # Spread over factors
  df <- tidyr::spread(df, key="factor", value="value")
  
  # Prepare the legend
  p <- ggplot(df, aes_string(x=factors[1], y=factors[2], color="color_by", shape="shape_by")) +
    geom_point() +
    scale_shape_manual(values=1:nlevels(df %>% pull(shape_by) %>% as.factor()))
    theme(
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_text(size=rel(1.2))
    )
  if (length(unique(df$color))>1) { p <- p + labs(color=color_name) } else { p <- p + guides(color=FALSE) + scale_color_manual(values="black") }
  if (is.numeric(df$color)) p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n=5, name="RdYlBu")))(10)) 
  if (length(unique(df$shape))>1) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (length(unique(df$color))>1 || length(unique(df$shape))>1) { legend <- GGally::grab_legend(p) } else { legend <- NULL }
  
  
  # Generate plot
  p <- GGally::ggpairs(df, 
                       columns = factors,
                       lower = list(continuous=GGally::wrap("points", size=dot_size)), 
                       diag = list(continuous='densityDiag'), 
                       upper = list(continuous=GGally::wrap("points", size=dot_size)), 
                       mapping = aes_string(color="color_by", shape="shape_by"), 
                       title = "", 
                       legend = legend
  )
  p <- p + theme_bw() + theme(
    # axis.line = element_line(color="black", size=rel(1.0)),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
  
  return(p)
}

.set_shapeby <- function (object, shape_by) 
{
  if (is.null(shape_by)) {
    shape_by <- rep("1", sum(object@dimensions[["N"]]))
  }
  else if (shape_by[1] == "group") {
    shape_by <- samples_metadata(object)$group
  }
  else if ((length(shape_by) == 1) && is.character(shape_by) & 
           (shape_by %in% colnames(samples_metadata(object)))) {
    shape_by <- samples_metadata(object)[, shape_by]
  }
  else if ((length(shape_by) == 1) && is.character(shape_by) && 
           (shape_by[1] %in% unlist(features_names(object)))) {
    viewidx <- which(sapply(features_names(object), function(x) shape_by %in% 
                              x))
    foo <- list(shape_by)
    names(foo) <- names(viewidx)
    shape_by <- lapply(get_data(object, features = foo), 
                       function(l) Reduce(cbind, l))[[1]][1, ]
  }
  else if (is(shape_by, "data.frame")) {
    stopifnot(all(colnames(shape_by) %in% c("sample", "color")))
    stopifnot(all(unique(shape_by$sample) %in% unlist(samples_names(object))))
  }
  else if (length(shape_by) > 1) {
    stopifnot(length(shape_by) == sum(object@dimensions[["N"]]))
  }
  else {
    stop("'shape_by' was specified but it was not recognised, please read the documentation")
  }
  if (!is(shape_by, "data.frame")) {
    df = data.frame(sample = samples_metadata(object)$sample, 
                    shape_by = as.factor(shape_by), stringsAsFactors = FALSE)
  }
  return(df)
}

.set_colorby <- function (object, color_by) 
{
  if (is.null(color_by)) {
    color_by <- rep("1", sum(object@dimensions[["N"]]))
  }
  else if (color_by[1] == "group") {
    color_by <- samples_metadata(object)$group
  }
  else if ((length(color_by) == 1) && (is.character(color_by) | 
                                       is.factor(color_by)) & (color_by[1] %in% colnames(samples_metadata(object)))) {
    color_by <- samples_metadata(object)[, color_by]
  }
  else if ((length(color_by) == 1) && is.character(color_by) && 
           (color_by[1] %in% unlist(features_names(object)))) {
    viewidx <- which(sapply(features_names(object), function(x) color_by %in% 
                              x))
    foo <- list(color_by)
    names(foo) <- names(viewidx)
    color_by <- lapply(get_data(object, features = foo), 
                       function(l) Reduce(cbind, l))[[1]][1, ]
  }
  else if ((length(color_by) == 1) && is.character(color_by) && 
           (color_by[1] %in% colnames(get_factors(object)[[1]]))) {
    color_by <- do.call(rbind, get_factors(object))[, color_by]
  }
  else if (is(color_by, "data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("sample", "color")))
    stopifnot(all(unique(color_by$sample) %in% unlist(samples_names(object))))
  }
  else if (length(color_by) > 1) {
    stopifnot(length(color_by) == sum(get_dimensions(object)$N))
  }
  else {
    stop("'color_by' was specified but it was not recognised, please read the documentation")
  }
  if (!is(color_by, "data.frame")) {
    df <- data.frame(sample = samples_metadata(object)$sample, 
                     color_by = color_by, stringsAsFactors = FALSE)
  }
  if (length(unique(df$color_by)) < 5) 
    df$color_by <- as.factor(df$color_by)
  return(df)
}

.check_and_get_factors <- function (object, factors) 
{
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(factors)))
  if (is.numeric(factors)) {
    stopifnot(all(factors <= object@dimensions$K))
    factors_names(object)[factors]
  }
  else {
    if (paste0(factors, collapse = "") == "all") {
      factors_names(object)
    }
    else {
      stopifnot(all(factors %in% factors_names(object)))
      factors
    }
  }
}


.select_stroke <- function (N) 
{
  if (N <= 1000) {
    stroke <- 0.5
  }
  else if (N > 1000 & N <= 10000) {
    stroke <- 0.2
  }
  else {
    stroke <- 0.05
  }
}


.add_legend <- function (p, df, legend, color_name, shape_name) 
{
  if (is.numeric(df$color_by)) {
    p <- p + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 5, 
                                                                           name = "RdYlBu")))(10)) + labs(fill = color_name)
  }
  else {
    if (length(unique(df$color_by)) > 1) {
      p <- p + guides(fill = guide_legend(override.aes = list(shape = 21))) + 
        labs(fill = color_name)
    }
    else {
      p <- p + guides(fill = "none", color = "none") + 
        scale_color_manual(values = "black") + scale_fill_manual(values = "gray60")
    }
  }
  if (length(unique(df$shape_by)) > 1) {
    p <- p + scale_shape_manual(values = c(21, 23, 24, 25)[1:length(unique(df$shape_by))]) + 
      guides(shape = guide_legend(override.aes = list(fill = "black"))) + 
      labs(shape = shape_name)
  }
  else {
    p <- p + scale_shape_manual(values = c(21)) + guides(shape = "none")
  }
  if (legend) {
    p <- p + guides(color = guide_legend(override.aes = list(fill = "white"))) + 
      theme(legend.text = element_text(size = rel(0.8)), 
            legend.title = element_text(size = rel(0.8)), 
            legend.key = element_rect(fill = "white", color = "white"))
  }
  else {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

###########
### PCA ###
###########

plot_pca_continous_color <- function(df_topca, e, joinby, mycolor, myshape, 
                     mytitle, axis1, axis2){
  
  pca.res <- prcomp(df_topca[,-1], center = T)
  df_e <- left_join(df_topca, e, by = joinby) %>% select(all_of(c(mycolor, myshape)))
  
  pca_data <- as.data.frame(pca.res$x)
  pca_data_e <- cbind(pca_data, df_e)
  #browser()
  scoreplot <- ggplot(pca_data_e,
                      aes(x = .data[[axis1]], y = .data[[axis2]], color = .data[[mycolor]], shape = .data[[myshape]])) + 
    geom_point(size = 2.5) + 
    scale_color_gradient(low = "#32CD32", high = "red") +
    ggtitle(mytitle) +
    theme_classic() +
    theme(axis.text.x= element_text(colour = 'black'),
          axis.text.y = element_text(colour = 'black'),
          plot.title = element_textbox_simple())
  
  if (!is.null(myshape)){
    scoreplot <- scoreplot + scale_shape_manual(values=1:nlevels(pca_data_e %>% pull(myshape) %>% as.factor()))
  }

  scoreplot
  
}

plot_pca2 <- function(df_topca, e, joinby, mycolor, myshape = NULL,
                      mytitle, toscale = FALSE, sel_pcs = 1:5, mypalette) {

  # Set rownames
  df_topca <- df_topca %>% column_to_rownames(joinby)
  
  # PCA
  pca.res <- prcomp(df_topca, center = TRUE, scale. = toscale)
  
  # Get variation for the figure
  pca_var <- pca.res$sdev^2
  pca_var_exp <- round(100 * pca_var / sum(pca_var), 1)
  
  new_names <- paste0("PC", 1:length(pca_var_exp), " (", pca_var_exp, "%)")
  colnames(pca.res$x) <- new_names
  pc_labels <- new_names[sel_pcs] 
  
  # Join scores with metadata
  pcaScores <- pca.res$x %>% 
    as.data.frame() %>% 
    mutate(join_id = rownames(.)) %>% 
    left_join(e, by = c('join_id' = joinby))
  #browser()
  # make shape/color columns factors if needed
  if (!is.null(myshape)) {
    pcaScores[[myshape]] <- as.factor(pcaScores[[myshape]])
  }

  pcaScores[[mycolor]] <- as.factor(pcaScores[[mycolor]])
  
  # Build aesthetics string for both color and shape
  mapping <- if (!is.null(myshape)) {
    ggplot2::aes_string(color = mycolor, shape = myshape)
  } else {
    ggplot2::aes_string(color = mycolor)
  }

  # Create plot
  scoreplot <- GGally::ggpairs(
    data = pcaScores,
    columns = pc_labels,
    mapping = mapping,
    upper = list(continuous = "blank")
  ) + 
    scale_color_brewer(palette = mypalette) + 
    ggtitle(mytitle) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # Dummy plot for legend
  if (!is.null(myshape)) {
    dummy <- ggplot(pcaScores, aes(
      x = !!sym(pc_labels[1]),
      y = !!sym(pc_labels[2]),
      color = !!sym(mycolor),
      shape = !!sym(myshape)
    )) +
      geom_point() +
      scale_color_brewer(palette = mypalette) + # ensure legend matches plot
      theme_void() +
      theme(legend.position = "right")
  } else {
    dummy <- ggplot(pcaScores, aes(
      x = !!sym(pc_labels[1]),
      y = !!sym(pc_labels[2]),
      color = !!sym(mycolor)
    )) +
      geom_point() +
      scale_color_brewer(palette = mypalette) + # ensure legend matches plot
      theme_void() +
      theme(legend.position = "right")
  }

  legend <- cowplot::get_legend(dummy)
  plot_grob <- GGally::ggmatrix_gtable(scoreplot)

  # Add a new column to the right for the legend
  combined <- gtable::gtable_add_cols(plot_grob, widths = unit(1.5, "in"))

  # Add legend to all rows, centered vertically
  combined <- gtable::gtable_add_grob(
    combined,
    grobs = legend,
    t = 1, l = ncol(combined), b = nrow(combined)
  )

  return(combined)
  
}



plot_pcas2 <- function(df_list, colors, shapes, e, joinby, 
                       namesfrom, valuesfrom,
                       savetodir = "results/plots/PCA",
                       w = 12, h = 10, discretecolor = T,
                       oneplot = F, toscale = F, sel_pcs = 1:5,
                       axis1 = 1, axis2 = 2, mypalette = "Paired"){
  # mycolor should be a vector of colors
  # oneplot means to output only one pca plot with specified axes

  p_list <- as.list(rep(1, length(df_list)))


  lapply(1:length(colors), function(i){
    colorx <- colors[i]
    shapex <- shapes[i]
    # If shapex is NA or NULL, set to NULL for downstream logic
    if (is.na(shapex) || is.null(shapex)) shapex <- NULL
    for (i in 1:length(df_list)){
      topca <- df_list[[i]] %>% 
        pivot_wider(names_from = all_of(namesfrom), values_from = all_of(valuesfrom), values_fill = 0)
      print(paste0("Working on: ", names(df_list[i]), ", color: ", colorx, ", shape: ", shapex))
      
      if(discretecolor){
        if (oneplot){
          combs <- combn(sel_pcs, 2, simplify = FALSE)
          # browser()
          lapply(combs, function(x){
            # make dir
            axisdir <- paste0(savetodir, "/ax", x[1], "-", x[2])
            if (!dir.exists(axisdir)){
              dir.create(axisdir, recursive = T)
            }
            myfilename = paste0(savetodir, "/ax", x[1], "-", x[2], "/", names(df_list[i]), "_PCA_color-", colorx,
                                "-shape-", shapex, ".png")
            myfilename_bip = paste0(savetodir, "/ax", x[1], "-", x[2], "/", names(df_list[i]), "_PCA_color-", colorx,
                                "-shape-", shapex, "_biplot.png")
            
            ps <- plot_pca_individual(df_topca = topca, e = e, joinby = joinby, mycolor = colorx, myshape = shapex, 
                                     mytitle = names(df_list[i]), axis1 = x[1], axis2 = x[2], mypalette = mypalette)
            ggsave(myfilename, ps$s, width = w, height = h, units = "cm")
            ggsave(myfilename_bip, ps$b, width = w, height = h, units = "cm")
            p_list[[i]] <- p
          })
        } else {
          # multiple plots with GGally
          myfilename = paste0(savetodir, "/", names(df_list[i]), "_PCA_color-", colorx, "_pcs", paste(sel_pcs, collapse = "-"), 
                              "-shape-", shapex, "_scale", toscale, ".png")
          mytitle = paste0( names(df_list[i]), " - ", colorx, " - ", shapex)
          p <- plot_pca2(df_topca = topca, e = e, joinby = joinby, mycolor = colorx, myshape = shapex,
                         mytitle = mytitle, toscale = toscale, sel_pcs = sel_pcs, mypalette = mypalette)
          ggsave(myfilename, p, width = w, height = h, units = "cm")
          p_list[[i]] <- p
        }
      } else{
        # Continous color pca
        myfilename = paste0(savetodir, "/ax", axis1, "-", axis2, "/", names(df_list[i]), "_PCA_color-", colorx,
                            "-shape-", shapex, ".png")
        p <- plot_pca_continous_color(df_topca = topca, e = e, joinby = joinby, mycolor = colorx, myshape = shapex, 
                                      mytitle = names(df_list[i]), axis1 = axis1, axis2 = axis2)
        ggsave(myfilename, p, width = w, height = h, units = "cm")
        p_list[[i]] <- p
      }
      
      
      
      
      
      names(p_list) <- paste0(names(df_list), "_p")
      p_list
    }
  })
}



plot_pca_individual <- function(df_topca, e, joinby, mycolor, myshape,
                     mytitle, axis1, axis2, frame = F, mypalette){
  #browser()
  pca.res <- prcomp(df_topca %>% column_to_rownames(joinby), center = T)
  biplot <- factoextra::fviz_pca_biplot(pca.res, repel = TRUE,
                                        col.var = "#2E9FDF", col.ind = "#696969",
                                        select.var = list(contrib = 20),
                                        axes = c(axis1, axis2)) +
    theme_bw()
  
  
  df_e <- left_join(df_topca, e, by = joinby)
  
  scoreplot <- autoplot(pca.res, data = df_e,
                        color = mycolor, size = 2.5, shape = myshape,
                        x = axis1, y = axis2, scale = 0,
                        frame = frame, frame.colour = myshape, frame.alpha = 0.05) +
    #scale_x_continuous(limits = c(-0.3, 1)) +
    #xlab("Dimension 1 (variance = 12.4%)") +
    scale_color_brewer(palette = mypalette) + 
    ggtitle(mytitle) +
    theme_classic() +
    theme(axis.text.x= element_text(colour = 'black'),
          axis.text.y= element_text(colour = 'black'),
          plot.title = element_textbox_simple())

  # Only add shape scale if myshape is not NULL and is a valid column in e
  if (!is.null(myshape) && myshape %in% colnames(e)) {
    scoreplot <- scoreplot + scale_shape_manual(values=1:nlevels(df_e %>% pull(myshape) %>% as.factor()))
  }

  return(list(s = scoreplot, b = biplot))

}

plot_pcas <- function(df_list, mycolor, myshape, e, joinby,
                      namesfrom, valuesfrom,
                      axis1 = 1, axis2 = 2, savetodir = "results/plots/PCA",
                      w = 12, h = 10, discretecolor = T){
  p_list <- as.list(rep(1, length(df_list)))

  if (!dir.exists(savetodir)){
    dir.create(savetodir, recursive = T)
  }

  for (i in 1:length(df_list)){
    myfilename = paste0(savetodir, "/", names(df_list[i]), "_PCA_color-", mycolor,
                        "-shape-", myshape, "-ax-", axis1, "-", axis2, ".png")
    topca <- df_list[[i]] %>%
      pivot_wider(names_from = all_of(namesfrom), values_from = all_of(valuesfrom), values_fill = 0)
    print(paste0("Working on: ", names(df_list[i])))
    if(discretecolor){
      p <- plot_pca(df_topca = topca, e = e, joinby = joinby, mycolor = mycolor, myshape = myshape,
                    mytitle = names(df_list[i]), axis1 = axis1, axis2 = axis2)
    } else{
      p <- plot_pca_continous_color(df_topca = topca, e = e, joinby = joinby, mycolor = mycolor, myshape = myshape,
                    mytitle = names(df_list[i]), axis1 = axis1, axis2 = axis2)
    }

    p_list[[i]] <- p
    ggsave(myfilename, p, width = w, height = h, units = "cm")
  }
  names(p_list) <- paste0(names(df_list), "_p")
  p_list
}



# 
# do_venn <- function(list_of_sets, savepath, mytitle, w = 11, h = 8){
#   p <- ggVennDiagram(list_of_sets,
#                      label = "count",
#                      label_alpha = 0 ) +
#     scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
#     scale_x_continuous(expand = expansion(mult = .5)) +
#     ggtitle(mytitle)
# 
#   filename = paste0(savepath, "/Venn ", paste0(names(list_of_sets), collapse = "-"), "-", mytitle, '.jpg')
#   
#   ggsave(filename, p, width = 11, height = 8, units = "cm")
# }
#   ggsave(filename, p, width = 11, height = 8, units = "cm")
# }
#   ggsave(filename, p, width = 11, height = 8, units = "cm")
# }
# }
#   ggsave(filename, p, width = 11, height = 8, units = "cm")
# }
#   ggsave(filename, p, width = 11, height = 8, units = "cm")
# }
